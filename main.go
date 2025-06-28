package main

import (
	"fmt"
	"os"
	"strings"

	"lab_buddy_go/benchmark"
	version_control "lab_buddy_go/config"
	"lab_buddy_go/fasta3bit"
	"lab_buddy_go/fasta_overview"
	"lab_buddy_go/kmer_analyzer"
	"lab_buddy_go/orf_finder"
	"lab_buddy_go/ran_dna_gen"
	"lab_buddy_go/sanity_check"
)

// printCustomHelp formats a custom help menu
func printCustomHelp() {
	fmt.Println(`Lab Buddy - Custom Help Menu
Usage:
  lab_buddy <tool> [options]

Tools:
  kmer_analyzer		Analyze k-mer frequencies
  orf_finder		Find open reading frames
  sanity_check		Run diagnostic test
  ran_dna_gen		Generate random DNA sequence
  fasta_overview	Summary statistics of FASTA file
  fasta3bit		(Beta) Encode DNA in binary for rapid analysis

Global Flags:
  -h, -help		Show this help message
  -v, -version		Show version information

Benchmarking:
  -benchmark		Must be used in associtation with a tool.
			Displays computational resource usage and 
			pertinent operating system information
  `,
)
	os.Exit(0)
}

func printVersion() {
	fmt.Println("Lab Buddy - Version Information Menu")
	fmt.Println("Central Executable:")
	fmt.Printf("\tLab Buddy:\t\t%s\n", version_control.Main_version)
	fmt.Printf("\nModular tools:\n")
	fmt.Printf("\tKmer Analyzer:\t\t%s\n", version_control.Kmer_Analyzer)
	fmt.Printf("\tORF Finder:\t\t%s\n", version_control.ORF_Finder)
	fmt.Printf("\tRandom DNA Generator:\t%s\n", version_control.Ran_DNE_Gen)
	fmt.Printf("\tSanity Check:\t\t%s\n", version_control.Sanity_check)
	fmt.Printf("\tFASTA Overview:\t\t%s\n", version_control.FASTA_Overview)
	fmt.Printf("\tFASTA 3 Bit Encoder:\t%s\n", version_control.FASTA_3_Bit)
	fmt.Printf("\tBenchmark:\t\t%s\n", version_control.Benchmark)
	
	fmt.Println("")

	os.Exit(0)
}

// Main controller 
func main() {

	// If no arguments are given, show help
	if len(os.Args) < 2 {
		printCustomHelp()
	}

	// Scan for executible-specific help flags
	for _, arg := range os.Args[1:] {
		if len(os.Args) < 3 {
			if arg == "-h" || arg == "-help" {
				printCustomHelp()
			}
		}
	}

	// Version request
	for _, arg := range os.Args[1:] {
		if arg == "-v" || arg == "-version" {
			printVersion()
		}
	}

	// 
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
