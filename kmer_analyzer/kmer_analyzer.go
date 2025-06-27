package kmer_analyzer

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"
)

// define_mer_pairs returns all possible k-length nucleotide strings (k-mers)
// using A, C, G, and T. If includeN is true, 'N' is also included.
func define_mer_pairs(k_value int, includeN bool) []string {
	nucleotides := []rune{'A', 'C', 'G', 'T'}		// Standard nucleotide options 
	if includeN {
		nucleotides = append(nucleotides, 'N')		// optionally include single ambigious base 
	}
	var kmers []string								// String slice to hold all possible kmers
	
	// Recursive function to build k-mers one base at a time.
	// Prefix: partial k-mer (built so far)
	// Depth: how many positions remain to reach full k-mer length
	var build func(prefix string, depth int)		
	build = func(prefix string, depth int) {		
		if depth == 0 {								// If k-mer is complete:
			kmers = append(kmers, prefix)			// Add it to the result slice
			return									// Move onto the next kmer
		}
		for _, base := range nucleotides {			// Loop through each nucleotide:
			build(prefix+string(base), depth-1)		// Recurse with new base and reduced depth
		}
	}
	build("", k_value)								// Start recursion with empty string and full depth
	return kmers									// When finished, return all generated k-mers
}


// reverseComplement returns the complimentary values for a provided DNA sequence in the negative direction.
func reverseComplement(seq string) string {
	var rc strings.Builder							// Initializes string builder
	for i := len(seq) - 1; i >= 0; i-- {			// Starts at the end and works backwards
		switch seq[i] {								// Processess by base
		case 'A': rc.WriteByte('T')
		case 'C': rc.WriteByte('G')
		case 'G': rc.WriteByte('C')
		case 'T': rc.WriteByte('A')
		default:  rc.WriteByte('N') 				// Fallback for ambigous characters 
		}
	}
	return rc.String()								// Return inverted, complimentary sequence
}


// countKmers returns k-mer frequencies in a FASTA file, along with the total number of valid k-mers found.
// It processes the FASTA file line-by-line and uses a rolling window to avoid loading the sequence into memory.
// If ignoreNs is true, k-mers containing 'N' are excluded.
func countKmers(filename string, k int, ignoreNs bool, strand string, frame int) (map[string]int, int, error) {
	file, err := os.Open(filename)					// Attempt to open the file
	if err != nil {
		return nil, 0, err							// Return error if file cannot be opened
	}
	defer file.Close()								// Ensure the file is closed when the function exits

	kmerCounts := make(map[string]int)				// Map to store k-mer (string) counts (int)
	total := 0										// Count of total valid k-mers
	var buffer []rune								// Rolling window of current sequence
	position := 0									// Tracks base position in sequence for frame tracking

	invalidBases := make(map[rune]int)				// Map of invalid bases detected (e.g., 'R')

	scanner := bufio.NewScanner(file)				// Read input line-by-line
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())	// Remove whitespace
		if strings.HasPrefix(line, ">") {			// If header is detected:
			continue 								// skip headers
		}

		// Only allow valid characters (including ambiguous base N if ignoreNs == false)
		for _, base := range strings.ToUpper(line) {	// Parses each base (uppercased) from the current sequence line
			if !strings.ContainsRune("ACGTN", base) {	// Check for invalid bases
				invalidBases[base]++
				continue 							// skip invalid characters
			}

			buffer = append(buffer, base)			// Appends the next base to the buffer
			if len(buffer) > k {					// If the buffer becomes longer than 'k':
				buffer = buffer[1:]					// Discard the leftmost (oldest) base to maintain k-length window
			}

			if len(buffer) == k {					// If the buffer is full:
				if frame == 0 || (position%3) == (frame-1) {		// Only count the kmer if it is in the appropriate frame
					kmer := string(buffer)				// Extract the kmer string
					if ignoreNs && strings.Contains(kmer, "N") {	// If ambigous base is detected:
						position++					// Move past the position without noting the kmer
						continue
					}

					switch strand {					// Handle strand-specific counting
					case "pos":						// If user specifies positive strand
						kmerCounts[kmer]++			// Add the kmer directly
					case "neg":						// If user specifies negative strand
						kmerCounts[reverseComplement(kmer)]++	// Reverse compliment the kmer, then add it
					default:						// Return error if invalid strand argument is provided
						return nil, 0, fmt.Errorf("invalid strand: %s", strand)
					}
					total++							// Increase total kmer count
				}
			}
			position++								// Move to the next position
		}
	}

	if len(invalidBases) > 0 {						// Display warning if invalid bases were detected
		fmt.Fprintln(os.Stderr, "Warning: FASTA input contains non-standard bases:")
		for base, count := range invalidBases {
			fmt.Fprintf(os.Stderr, "  %c: %d occurrences\n", base, count)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, 0, err
	}

	return kmerCounts, total, nil
}


// Run executes the kmer_analyzer command. 
// It expects a FASTA file and a k-mer size via command-line arguments,
// and prints the frequency of all k-mers found in the input sequence.
// Optionally reports all possible k-mers without quantity.
func Run_kmer_analyzer(args []string) {

	fs := flag.NewFlagSet("kmer_analyzer", flag.ExitOnError) 	// Isolated flag set specifically for "kmer_analyzer" subcommand 

	k_value := fs.Int("k_mer", 3, "K-mer value")	// Size of K-mer. 
	in_file := fs.String("in_file", "", "FASTA file input")		// Input file (FASTA)
	report_kmers := fs.Bool("report_kmer", false, "List all possible k-mers only")	// Option to generate and report all possible k-mers without frequency
	rel_freq := fs.Bool("rel_freq", true, "Output relative frequency (%)")			// Output relative frequency (%) if true (default: true)
	sort_by := fs.String("sort_by", "alpha", "Sort output by 'alpha' or 'freq'")	// Output sorting option for by alphabetical or by frequency 
	ignoreNs := fs.Bool("ignore_ns", false, "Ignore k-mers containing N")			// Option to ignore results with ambigous nucleotides
	frame := fs.Int("frame", 0, "Reading frame (0 = all (default), 1, 2, 3)")		// Optional frame-specific behavior (default '0' - All frames)
	strand := fs.String("strand", "pos", "Strand direction: pos, neg")				// Strand-specific directionality
	outFile := fs.String("out_file", "", "Optional: path to save output instead of printing to terminal") 	// Optional output file

	fs.Parse(args)					// Parse the flag set arguments

	allKmers := define_mer_pairs(*k_value, !*ignoreNs)	// Generate all possible kmers (+/- 'N')

	if *report_kmers {								// If user requests raw kmers:
		fmt.Println("All possible k-mers:")				
		fmt.Println(allKmers)						// Report all possible kmers without frequencies
		return
	}

	if *frame < 0 || *frame > 3 {					// Frame validation
		fmt.Println("Error: -frame must be 0 (all), 1, 2, or 3")
		os.Exit(1)
	}

	if *strand != "pos" && *strand != "neg" {		// Strand validation
		fmt.Println("Error: -strand must be 'pos' or 'neg'")
		os.Exit(1)
	}

	if *in_file == "" {								// User needs to provide a FASTA input or request raw kmers
		fmt.Println("Error: -in_file is required when not using -report_kmer")
		os.Exit(1)
	}

	kmerCounts, total, err := countKmers(*in_file, *k_value, *ignoreNs, *strand, *frame)		// Detects and counts relevant kmers
	if err != nil {
		fmt.Println("Error:", err)
		return
	}

	type kmerData struct {		// Declares struct for data organization
		Kmer   string			// Kmers 
		Count  int				// Kmer counts
		RelPct float64			// Percentage of kmers out of total
	}

	var result []kmerData		// Slice to hold merged k-mer results with count and percentage
	for _, kmer := range allKmers {		// For all generated kmers:
		count := kmerCounts[kmer]		// Get observed count; defaults to 0 if k-mer was not found
		pct := 0.0						// Initialize percentages
		if total > 0 {					// If kmers were detected:
			pct = float64(count) / float64(total) * 100		// Update percentage 
		}
		result = append(result, kmerData{kmer, count, pct})	// Prepares results
	}

	switch *sort_by {								// Output reporting option (sorting)				
	case "freq":											// Option to sort by kmer prevalence 
		sort.Slice(result, func(i, j int) bool {
			return result[i].Count > result[j].Count
		})
	default:												// Option to sort alphabetically by kmer
		sort.Slice(result, func(i, j int) bool {
			return result[i].Kmer < result[j].Kmer
		})
	}
	
	var out *os.File										// Define output file (if needed)
	if *outFile != "" {										// If outFile is not empty:
		var err error
		out, err = os.Create(*outFile)						// Create output file
		if err != nil {										// Return error creating output file if needed
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer out.Close()									// Ensure file is properly closed
	} else {
		out = os.Stdout										// If no outfile is provided, print to terminal
	}
	
	if *rel_freq {
		fmt.Fprintln(out, "K-mer\tCount\tRelative_Freq(%)")	// Print appropriate header
	} else {
		fmt.Fprintln(out, "K-mer\tCount")
	}
	
	for _, item := range result {							// Print Kmer results
		if *rel_freq {
			fmt.Fprintf(out, "%s\t%d\t%.2f\n", item.Kmer, item.Count, item.RelPct)
		} else {
			fmt.Fprintf(out, "%s\t%d\n", item.Kmer, item.Count)
		}
	}
	
}
