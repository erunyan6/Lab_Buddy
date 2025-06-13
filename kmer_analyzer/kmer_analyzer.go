package kmer_analyzer

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"
)

func define_mer_pairs(k_value int, includeN bool) []string {
	nucleotides := []rune{'A', 'C', 'G', 'T'}
	if includeN {
		nucleotides = append(nucleotides, 'N')
	}
	var kmers []string
	var build func(prefix string, depth int)
	build = func(prefix string, depth int) {
		if depth == 0 {
			kmers = append(kmers, prefix)
			return
		}
		for _, base := range nucleotides {
			build(prefix+string(base), depth-1)
		}
	}
	build("", k_value)
	return kmers
}

func countKmers(filename string, k int, ignoreNs bool) (map[string]int, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, err
	}
	defer file.Close()

	var sequenceBuilder strings.Builder
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if !strings.HasPrefix(line, ">") {
			sequenceBuilder.WriteString(strings.ToUpper(line))
		}
	}
	sequence := sequenceBuilder.String()
	kmerCounts := make(map[string]int)
	total := 0

	for i := 0; i <= len(sequence)-k; i++ {
		kmer := sequence[i : i+k]
		if ignoreNs && strings.Contains(kmer, "N") {
			continue
		}
		kmerCounts[kmer]++
		total++
	}
	return kmerCounts, total, nil
}

func Run(args []string) {

	fs := flag.NewFlagSet("kmer_analyzer", flag.ExitOnError)

	k_value := fs.Int("k_mer", 3, "K-mer value (default: 3)")
	in_file := fs.String("in_file", "", "FASTA file input")
	report_kmers := fs.Bool("report_kmer", false, "List all possible k-mers only")
	rel_freq := fs.Bool("rel_freq", true, "Output relative frequency (%)")
	sort_by := fs.String("sort_by", "alpha", "Sort output by 'alpha' or 'freq'")
	ignoreNs := fs.Bool("ignore_ns", false, "Ignore k-mers containing N")
	
	fs.Parse(args)

	allKmers := define_mer_pairs(*k_value, !*ignoreNs)

	if *report_kmers {
		fmt.Println("All possible k-mers:")
		fmt.Println(allKmers)
		return
	}

	if *in_file == "" {
		fmt.Println("Error: --in_file is required when not using --report_kmer")
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
