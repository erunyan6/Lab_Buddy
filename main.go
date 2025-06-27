package main

import (
    "fmt"
    "os"
    "strings"

    "lab_buddy_go/sanity_check"
	"lab_buddy_go/ran_dna_gen"
    "lab_buddy_go/kmer_analyzer"
    "lab_buddy_go/orf_finder"
    "lab_buddy_go/fasta3bit"
    "lab_buddy_go/fasta_overview"
    "lab_buddy_go/benchmark"
    
)

// Main controller 
func main() {
    if len(os.Args) < 2 {
        fmt.Println("Usage: ./lab_buddy <tool> [flags]")
        os.Exit(1)
    }
    
    toolName := os.Args[1]
    toolArgs := os.Args[2:]

    // Check for global --benchmark flag
	benchmarking := false
	var cleanedArgs []string
	for _, arg := range toolArgs {
		if arg == "-benchmark" {
			benchmarking = true
		} else {
			cleanedArgs = append(cleanedArgs, arg)
		}
	}

	// Tool execution wrapper
	run := func() {
		switch toolName {
		case "ran_dna_gen":
			ran_dna_gen.Run(cleanedArgs)
		case "check":
			sanity_check.Run(cleanedArgs)
		case "kmer_analyzer":
			kmer_analyzer.Run_kmer_analyzer(cleanedArgs)
		case "orf_finder":
			orf_finder.Run(cleanedArgs)
		case "encoder":
			fasta3bit.Run(cleanedArgs)
		case "fasta_overview":
			fasta_overview.Run(cleanedArgs)
		default:
			fmt.Printf("Unknown tool: %s\n", toolName)
			os.Exit(1)
		}
	}

	if benchmarking {
		label := fmt.Sprintf("lab_buddy %s %s", toolName, strings.Join(cleanedArgs, " "))
		benchmark.Run(label, run)
	} else {
		run()
	}
}
