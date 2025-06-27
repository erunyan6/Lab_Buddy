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


// countKmers returns k-mer frequencies in a FASTA file, along with the total number of valid k-mers found.
// If ignoreNs is true, k-mers containing 'N' are excluded.
func countKmers(filename string, k int, ignoreNs bool) (map[string]int, int, error) {
	file, err := os.Open(filename)					// Attempt to open the file
	if err != nil {									
		return nil, 0, err							// Check for file readability
	}
	defer file.Close()								// Close to prevent file-handle leaks

	var sequenceBuilder strings.Builder				// Efficiently build the full DNA sequence
	scanner := bufio.NewScanner(file)				// Define memory-efficient line-by-line scanner
	for scanner.Scan() {							// Scan by line
		line := strings.TrimSpace(scanner.Text())	// Returns the line without white space
		if !strings.HasPrefix(line, ">") {			// If line is NOT a header:
			sequenceBuilder.WriteString(strings.ToUpper(line))		// Append sequence in uppercase
		}
	}
	sequence := sequenceBuilder.String()			// Convert builder to a single string
	kmerCounts := make(map[string]int)				// Defines map. kmer is key (string), kmer count is value pair (int)
	total := 0										// Track total valid k-mers counted

	for i := 0; i <= len(sequence)-k; i++ {			// Slides a window of size k across the sequence 
		kmer := sequence[i : i+k]					// Extract k-mer of length k
		if ignoreNs && strings.Contains(kmer, "N") {
			continue								// Skip k-mers with 'N' if requested
		}
		kmerCounts[kmer]++							// Increment k-mer count
		total++										// Increment total valid k-mers seen
	}
	return kmerCounts, total, nil					// Return results
}


// Run executes the kmer_analyzer command. 
// It expects a FASTA file and a k-mer size via command-line arguments,
// and prints the frequency of all k-mers found in the input sequence.
// Optionally reports all possible k-mers without quantity.
func Run_kmer_analyzer(args []string) {

	fs := flag.NewFlagSet("kmer_analyzer", flag.ExitOnError) 	// Isolated flag set specifically for "kmer_analyzer" subcommand 

	k_value := fs.Int("k_mer", 3, "K-mer value (default: 3)")	// Size of K-mer. 
	in_file := fs.String("in_file", "", "FASTA file input")		// Input file (FASTA)
	report_kmers := fs.Bool("report_kmer", false, "List all possible k-mers only")	// Option to generate and report all possible k-mers without frequency
	rel_freq := fs.Bool("rel_freq", true, "Output relative frequency (%)")			// Option to toggle kmer percentage
	sort_by := fs.String("sort_by", "alpha", "Sort output by 'alpha' or 'freq'")	// Output sorting option for by alphabetical or by frequency 
	ignoreNs := fs.Bool("ignore_ns", false, "Ignore k-mers containing N")			// Option to ignore results with ambigous nucleotides
	
	fs.Parse(args)					// Parse the flag set arguments

	allKmers := define_mer_pairs(*k_value, !*ignoreNs)	// Generate all possible kmers (+/- 'N')

	if *report_kmers {								// If user requests raw kmers:
		fmt.Println("All possible k-mers:")				
		fmt.Println(allKmers)						// Report all possible kmers without frequencies
		return
	}

	if *in_file == "" {								// User needs to provide a FASTA input or request raw kmers
		fmt.Println("Error: -in_file is required when not using -report_kmer")
		os.Exit(1)
	}

	kmerCounts, total, err := countKmers(*in_file, *k_value, *ignoreNs)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}

	type kmerData struct {
		Kmer   string
		Count  int
		RelPct float64
	}

	var result []kmerData
	for _, kmer := range allKmers {
		count := kmerCounts[kmer]
		pct := 0.0
		if total > 0 {
			pct = float64(count) / float64(total) * 100
		}
		result = append(result, kmerData{kmer, count, pct})
	}

	switch *sort_by {
	case "freq":
		sort.Slice(result, func(i, j int) bool {
			return result[i].Count > result[j].Count
		})
	default: // alpha
		sort.Slice(result, func(i, j int) bool {
			return result[i].Kmer < result[j].Kmer
		})
	}

	fmt.Println("K-mer\tCount\tRelative_Freq(%)")
	for _, item := range result {
		if *rel_freq {
			fmt.Printf("%s\t%d\t%.2f\n", item.Kmer, item.Count, item.RelPct)
		} else {
			fmt.Printf("%s\t%d\n", item.Kmer, item.Count)
		}
	}
}
