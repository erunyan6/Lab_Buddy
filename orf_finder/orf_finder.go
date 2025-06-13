package orf_finder

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"io"
)

func findORFs(sequence string, strand string, minLength int, showSeq bool, seqID string,
	originalLength int, strictMode bool, startCodons map[string]bool, frameSet map[string]bool,
	outFmt string, orfCounter *int, output io.Writer, summary *SummaryStats) {
	stopCodons := map[string]bool{"TAA": true, "TAG": true, "TGA": true}
	usedStops := make(map[int]bool)

	for i := 0; i <= len(sequence)-3; i++ {
		codon := sequence[i : i+3]
		if startCodons[codon] {
			for j := i + 3; j <= len(sequence)-3; j += 3 {
				stop := sequence[j : j+3]
				if stopCodons[stop] {
					if strictMode && usedStops[j] {
						break
					}
					usedStops[j] = true

					orfLength := j + 3 - i
					if orfLength >= minLength {
						var start, end, frame int
						if strand == "-" {
							start = originalLength - (j + 3) + 1
							end = originalLength - i
							frame = -((i % 3) + 1)
						} else {
							start = i + 1
							end = j + 3
							frame = (i % 3) + 1
						}

						summary.Total++
						summary.TotalLength += orfLength
						
						if strand == "+" {
							summary.Forward++
						} else {
							summary.Reverse++
						}
						
						frameString := fmt.Sprintf("%+d", frame)
						summary.FrameCounts[frameString]++
						
						if orfLength > summary.LongestLength {
							summary.LongestLength = orfLength
							summary.LongestSeqID = seqID
							summary.LongestStart = start
							summary.LongestEnd = end
						}						

						orfID := fmt.Sprintf("orf%d", *orfCounter)
						(*orfCounter)++

						if outFmt == "gff" {
							phase := (3 - (start-1)%3) % 3
							fmt.Fprintf(output, "%s\torf_finder\tORF\t%d\t%d\t.\t%s\t%d\tID=%s;Length=%d\n",
								seqID, start, end, strand, phase, orfID, orfLength)
						} else {
							if showSeq {
								orfSeq := sequence[i : j+3]
								fmt.Fprintf(output, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n", seqID, strand, start, end, orfLength, frame, orfSeq)
							} else {
								fmt.Fprintf(output, "%s\t%s\t%d\t%d\t%d\t%d\n", seqID, strand, start, end, orfLength, frame)
							}
						}
					}
					break
				}
			}
		}
	}
}

func reverseComplement(seq string) string {
	var rc strings.Builder
	for i := len(seq) - 1; i >= 0; i-- {
		switch seq[i] {
		case 'A':
			rc.WriteByte('T')
		case 'T':
			rc.WriteByte('A')
		case 'C':
			rc.WriteByte('G')
		case 'G':
			rc.WriteByte('C')
		default:
			rc.WriteByte('N')
		}
	}
	return rc.String()
}

func readMultiFasta(file string) (map[string]string, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	sequences := make(map[string]string)
	var currentID string
	var seqBuilder strings.Builder

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(line, ">") {
			if currentID != "" {
				sequences[currentID] = seqBuilder.String()
				seqBuilder.Reset()
			}
			currentID = strings.TrimPrefix(line, ">")
		} else {
			seqBuilder.WriteString(strings.ToUpper(line))
		}
	}
	if currentID != "" {
		sequences[currentID] = seqBuilder.String()
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return sequences, nil
}

type SummaryStats struct {
	Total         int
	Forward       int
	Reverse       int
	LongestLength int
	LongestSeqID  string
	LongestStart  int
	LongestEnd    int
	FrameCounts   map[string]int
	TotalLength   int
}


func Run(args []string) {

	fs := flag.NewFlagSet("orf_finder", flag.ExitOnError)

	// Define flags
	inputFile := fs.String("in_file", "", "Input FASTA file")
	minLen := fs.Int("minlen", 0, "Minimum ORF length")
	showSeq := fs.Bool("showseq", false, "Show ORF sequence. Suppressed on GFF3 output mode")
	mode := fs.String("mode", "thorough", "Mode: 'strict' or 'thorough'")
	startCodonsFlag := fs.String("start", "ATG", "Comma-separated list of start codons (e.g., ATG,GTG,TTG)")
	strandFlag := fs.String("strand", "both", "Strand to search: +, -, or both")
	frameFlag := fs.String("frame", "all", "Comma-separated frame(s): +1,+2,+3,-1,-2,-3 or all")
	outFmt := fs.String("outfmt", "tsv", "Output format: 'tsv' or 'gff'")
	outFile := fs.String("out", "", "Write output to file (optional)")
	summaryFlag := fs.Bool("summary", false, "Print ORF summary to stdout")
	fs.Parse(args)

	// Validate strand
	validStrands := map[string]bool{"+": true, "-": true, "both": true}
	if !validStrands[*strandFlag] {
		log.Fatalf("Invalid strand: %s (choose +, -, or both)", *strandFlag)
	}

	// Validate and parse frames
	validFrames := map[string]bool{
		"+1": true, "+2": true, "+3": true,
		"-1": true, "-2": true, "-3": true,
		"all": true,
	}
	frameSet := make(map[string]bool)
	if *frameFlag == "all" {
		frameSet["all"] = true
	} else {
		for _, f := range strings.Split(*frameFlag, ",") {
			f = strings.TrimSpace(f)
			if !validFrames[f] {
				log.Fatalf("Invalid frame: %s (choose from +1 to -3, or all)", f)
			}
			frameSet[f] = true
		}
	}

	summary := &SummaryStats{FrameCounts: make(map[string]int)}

	// Validate output format
	if *outFmt != "tsv" && *outFmt != "gff" {
		log.Fatalf("Invalid --outfmt: %s (choose 'tsv' or 'gff')", *outFmt)
	}

	// Parse start codons
	startCodons := make(map[string]bool)
	for _, codon := range strings.Split(*startCodonsFlag, ",") {
		codon = strings.ToUpper(strings.TrimSpace(codon))
		if len(codon) == 3 {
			startCodons[codon] = true
		}
	}

	// Validate input
	if *inputFile == "" {
		log.Fatal("Error: --in_file is required")
	}

	// Mode check
	strictMode := false
	if *mode == "strict" {
		strictMode = true
	} else if *mode != "thorough" {
		log.Fatalf("Unknown mode: %s. Use 'strict' or 'thorough'", *mode)
	}

	// Prepare output destination
	var output io.Writer = os.Stdout
	if *outFile != "" {
		f, err := os.Create(*outFile)
		if err != nil {
			log.Fatalf("Failed to create output file: %v", err)
		}
		defer f.Close()
		output = f
	}
	
	seqMap, err := readMultiFasta(*inputFile)
	if err != nil {
		log.Fatalf("Failed to read FASTA: %v", err)
	}
	
	if *outFmt == "gff" {
		fmt.Fprintln(output, "##gff-version 3")  // use output writer
	}
	
	orfCounter := 1
	for seqID, sequence := range seqMap {
		originalLength := len(sequence)
		if *strandFlag == "+" || *strandFlag == "both" {
			findORFs(sequence, "+", *minLen, *showSeq, seqID, originalLength, strictMode, startCodons, frameSet, *outFmt, &orfCounter, output, summary)
		}
		if *strandFlag == "-" || *strandFlag == "both" {
			rc := reverseComplement(sequence)
			findORFs(rc, "-", *minLen, *showSeq, seqID, originalLength, strictMode, startCodons, frameSet, *outFmt, &orfCounter, output, summary)
		}
	}

	if *summaryFlag {
		avg := 0.0
		if summary.Total > 0 {
			avg = float64(summary.TotalLength) / float64(summary.Total)
		}
	
		fmt.Fprintln(os.Stdout, "\n=== ORF Summary ===")
		fmt.Fprintf(os.Stdout, "Total ORFs: %d\n", summary.Total)
		fmt.Fprintf(os.Stdout, "  Forward strand: %d\n", summary.Forward)
		fmt.Fprintf(os.Stdout, "  Reverse strand: %d\n", summary.Reverse)
		fmt.Fprintf(os.Stdout, "Longest ORF: %d bp (%s:%d-%d)\n",
			summary.LongestLength, summary.LongestSeqID, summary.LongestStart, summary.LongestEnd)
		fmt.Fprintf(os.Stdout, "Average ORF length: %.1f bp\n", avg)
		fmt.Fprintln(os.Stdout, "Frame usage:")
		for frame, count := range summary.FrameCounts {
			fmt.Fprintf(os.Stdout, "  %s: %d\n", frame, count)
		}
	}
	
}
