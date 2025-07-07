package orf_finder

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"strconv"

	"lab_buddy_go/utils"
)

type ORF struct {
	SeqID  string
	Start    int
	End      int
	Strand   string
	Frame    int
	Length_nt   int
	Length_aa	int
	StartCodon string
}

func findORFs(seq_id string, seq string, frame []int, strand string, startCodons map[string]bool) []ORF {
	var orfs []ORF															// Slice to store identified ORFs

	stopCodons := map[string]bool{"TAA": true, "TAG": true, "TGA": true}	// Define stop codons

	s := strings.ToLower(strand)											// Variable to save strand option

	if s == "positive" || s == "both" {										// For ORFs on positive strand
		for _, f := range frame {											// For each frame 
			for i := f - 1; i <= len(seq)-3; i += 3 {						// Start at the 0 for specific frame and end at the last viable codon
				codon := seq[i : i+3] 										// Grab the whole 3-letter codon
				if startCodons[codon] {										// If the codon is in the start codons map:
					orfFound := false										// Track if stop codon was found
					for j := i + 3; j <= len(seq)-3; j += 3 {				// Scan plus 3 each iteration 
						stop := seq[j : j+3]								// Declare and update stop variable
						if stopCodons[stop] {								// If the current codon is in the stop codons map:
							start := i										// Save 'i' index as the start
							end := j + 3									// Save the last index of the stop codon to the end
							orfLength := end - start						// Calculate the length of the ORF

							orfs = append(orfs, ORF{						// Append data to ORF struct
								SeqID:     seq_id,							// Sequence name
								Start:     start,							// ORF start index
								End:       end,								// ORF end index 
								Strand:    "+",								// ORF strand
								Length_nt: orfLength,						// ORF length in nucleotides
								Length_aa: orfLength / 3,					// ORF length in amino acids
								Frame:     f,								// ORF frame
								StartCodon: codon,
							})
							orfFound = true
							break											// Move onto next start codon
						}
					}
					if !orfFound {
						start := i											// Save 'i' index as the start
						end := len(seq)										// Extend to end of sequence
						orfLength := end - start							// Calculate the length

						orfs = append(orfs, ORF{							// Append incomplete ORF
							SeqID:     seq_id,							
							Start:     start,
							End:       -5,									// Use -1 or placeholder for ">0"
							Strand:    "+",
							Length_nt: orfLength,
							Length_aa: orfLength / 3,
							Frame:     f,
							StartCodon: codon,
						})
					}
				}
			}
		}
	}

	if s == "negative" || s == "both" {										// For ORFs on negative strand
		rcSeq := common.ReverseComplement(seq)								// Compute reverse complement of sequence
		for _, f := range frame {
			for i := f - 1; i <= len(rcSeq)-3; i += 3 {
				codon := rcSeq[i : i+3]
				if startCodons[codon] {
					orfFound := false
					for j := i + 3; j <= len(rcSeq)-3; j += 3 {
						stop := rcSeq[j : j+3]
						if stopCodons[stop] {
							start := len(seq) - (j + 3)						// Convert reverse coords to original sequence
							end := len(seq) - i
							orfLength := end - start

							orfs = append(orfs, ORF{
								SeqID:     seq_id,
								Start:     start,
								End:       end,
								Strand:    "-",
								Length_nt: orfLength,
								Length_aa: orfLength / 3,
								Frame:     -f,
								StartCodon: codon,
							})
							orfFound = true
							break											// Move to next start codon
						}
					}
					if !orfFound {
						start := -5										// No stop codon found
						end := len(seq) - i								// Position of start codon in original strand
						orfLength := end - start							// Calculate approximate ORF length

						orfs = append(orfs, ORF{
							SeqID:     seq_id,
							Start:     start,
							End:       end,
							Strand:    "-",
							Length_nt: orfLength,
							Length_aa: (orfLength / 3) - 1,		// Minus one for stop codon
							Frame:     -f,
							StartCodon: codon,
						})
					}
				}
			}
		}
	}

	return orfs
}

func orfHandler(id string, seq string, opts map[string]interface{}) error {
	frames := opts["frames"].([]int)							// List of frames to check
	strand := opts["strand"].(string)							// Strand option
	minLen := opts["minLen"].(int)								// Minimum ORF length
	startCodons := opts["start_codons"].(map[string]bool)
	orfs := findORFs(id, seq, frames, strand, startCodons)					// Run ORF finder

	offset := 0
	if val, ok := opts["chunk_start"].(int); ok {
		offset = val
	}

	suppInc := false											// Set default: include incomplete ORFs
	if val, ok := opts["supp_inc"].(bool); ok {
		suppInc = val											// Update if flag provided
	}

	writer := opts["writer"].(*bufio.Writer)					// Output writer (stdout or file)

	for i, orf := range orfs {
		if suppInc && (orf.Start == -5 || orf.End == -5) {
			continue											// Skip incomplete ORFs if user requests suppression
		}
		if orf.Length_nt >= minLen {

			// GFF3 uses 1-based start coordinates
			start := orf.Start
			end := orf.End
			
			if start != -5 {
				start += offset
			}
			if end != -5 {
				end += offset
			}

			// Set phase (0-based codon offset)
			absFrame := orf.Frame
			if absFrame < 0 {
				absFrame = -absFrame
			}
			phase := (absFrame - 1) % 3
			if phase < 0 {
				phase += 3
			}		

			// Build attribute string
			attrs := fmt.Sprintf(
				"ID=orf%d;Length_nt=%d;Length_aa=%d;Frame=%d;StartCodon=%s",
				i+1, orf.Length_nt, orf.Length_aa, orf.Frame, orf.StartCodon,
			)			

			if orf.Start == -5 || orf.End == -5 {
				attrs += ";Partial=Yes"							// Add Partial flag for incomplete ORFs
			}

			// Construct GFF3 line
			gffLine := fmt.Sprintf(
				"%s\tLabBuddy\tORF\t%d\t%d\t.\t%s\t%d\t%s\n",
				orf.SeqID,
				start+1,											// Convert to 1-based
				end,
				orf.Strand,
				phase,
				attrs,
			)
			writer.WriteString(gffLine)
		}
	}

	return nil
}


func parseFrames(frameStr string) []int {
	var frames []int
	for _, s := range strings.Split(frameStr, ",") {
		s = strings.TrimSpace(s)
		if f, err := strconv.Atoi(s); err == nil {
			frames = append(frames, f)
		}
	}
	return frames
}


func Run(args []string) {
	fs := flag.NewFlagSet("orf_finder", flag.ExitOnError)

	inputFile := fs.String("in_file", "", "Input FASTA file")
	minLen := fs.Int("minlen", 100, "Minimum ORF length")
	frameFlag := fs.String("frame", "1,2,3", "Comma-separated frame(s): 1,2,3")
	strand := fs.String("strand", "both", "DNA directionality for analysis (both/positive/negative)")
	outFile := fs.String("out_file", "", "Output file (default is stdout)")
	suppInc := fs.Bool("supp_inc", false, "Suppress incomplete ORFs (those without stop codons)")
	startCodonsFlag := fs.String("start", "ATG", "Comma-separated list of start codons (e.g., ATG,GTG,TTG)")

	err := fs.Parse(args)
	if err != nil {
		fmt.Println("Error parsing flags:", err)
		os.Exit(1)
	}
	if len(fs.Args()) > 0 {
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}
	if *inputFile == "" {
		log.Fatal("Error: -in_file is required")
	}

	validBases := map[rune]bool{'A': true, 'T': true, 'G': true, 'C': true}

	codonSet := make(map[string]bool)
	for _, codon := range strings.Split(strings.ToUpper(*startCodonsFlag), ",") {
		codon = strings.TrimSpace(codon)
		if len(codon) != 3 {
			fmt.Fprintf(os.Stderr, "Warning: start codon %q ignored (not 3 letters)\n", codon)
			continue
		}
		valid := true
		for _, base := range codon {
			if !validBases[base] {
				valid = false
				break
			}
		}
		if !valid {
			fmt.Fprintf(os.Stderr, "Warning: start codon %q ignored (contains non-ATCG letters)\n", codon)
			continue
		}
		codonSet[codon] = true
	}
	if len(codonSet) == 0 {
		fmt.Fprintln(os.Stderr, "No valid start codons provided. Defaulting to ATG.")
		codonSet["ATG"] = true
	}
	

	frames := parseFrames(*frameFlag)
	for _, f := range frames {
		if f <= 0 || f > 3 {
			log.Fatalf("Invalid frame: %d. Only 1, 2, 3 are allowed.", f)
		}
	}

	s := strings.ToLower(*strand)
	acceptableStrand := map[string]bool{"positive": true, "negative": true, "both": true}
	if !acceptableStrand[s] {
		log.Fatalf("Invalid strand: %s. Allowed values are 'positive', 'negative', or 'both'.", *strand)
	}

	var writer *bufio.Writer

	if *outFile == "" {
		// Default to stdout
		writer = bufio.NewWriter(os.Stdout)
	} else {
		// Write to specified file
		file, err := os.Create(*outFile)
		if err != nil {
			log.Fatalf("Failed to create output file: %v", err)
		}
		writer = bufio.NewWriter(file)
		defer file.Close() // Close after the run is complete
	}	

	opts := map[string]interface{}{
		"frames": frames,
		"strand": *strand,
		"minLen": *minLen,
		"writer": writer,
		"supp_inc": *suppInc,
		"start_codons": codonSet,
	}

	writer.WriteString("##gff-version 3\n")

	err = common.StreamFastaWithOpts(*inputFile, orfHandler, opts)
	if err != nil {
		log.Fatalf("error running ORF finder: %v", err)
	}

	writer.Flush()
}
