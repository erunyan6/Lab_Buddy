package main

import (
	"fmt"
	"os"
	"strings"

	"lab_buddy_go/tools/benchmark"
	"lab_buddy_go/config"
	"lab_buddy_go/tools/fasta_overview"
	"lab_buddy_go/tools/kmer_analyzer"
	"lab_buddy_go/tools/orf_finder"
	"lab_buddy_go/tools/seq_generator"
	"lab_buddy_go/tools/sanity_check"
	"lab_buddy_go/tools/lab_buddy_art"
	"lab_buddy_go/tools/fasta_indexer"
	"lab_buddy_go/tools/orf_to_faa"
	"lab_buddy_go/tools/seq_sim"
	"lab_buddy_go/tools/fastqc_mimic"
	"lab_buddy_go/tools/fasta_isolate"
)

// printCustomHelp formats a custom help menu
func printCustomHelp() {
	fmt.Println(`Lab Buddy - Help Menu

Usage:
  lab_buddy <tool> [options]

Tools:
  kmer_analyzer		Analyze k-mer frequencies
  orf_finder		Find open reading frames
  check			Run diagnostic test
  seq_gen		Generate random DNA/RNA/Protein sequence(s)
  fasta_overview	Summary statistics of FASTA file
  lab_buddy_art		Cute and Fun ASCII art of Lab Buddy himself with an encouraging quote
  index_fasta 		Index FASTA for easy sequence access
  orf_to_faa        	Translate ORFs from orf_finder into FAA format
  seq_sim		Lightweight sequencing simulator for simple reads
  fastqc_mimic		Lab_Buddy version of the popular FASTQC analyzer and report generator
  fasta_isolate		Rapidly extract specific entries / ranges from FASTA files

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
	fmt.Printf("  Lab Buddy:\t\t%s\n", version_control.Main_version)
	fmt.Printf("\nModular tools:\n")
	fmt.Printf("  Kmer Analyzer:\t%s\n", version_control.Kmer_Analyzer)
	fmt.Printf("  ORF Finder:\t\t%s\n", version_control.ORF_Finder)
	fmt.Printf("  Seq Generator:\t%s\n", version_control.Seq_Generator)
	fmt.Printf("  Sanity Check:\t\t%s\n", version_control.Sanity_check)
	fmt.Printf("  FASTA Overview:\t%s\n", version_control.FASTA_Overview)
	fmt.Printf("  Benchmark:\t\t%s\n", version_control.Benchmark)
	fmt.Printf("  Lab Buddy Art\t\t%s\n", version_control.Lab_Buddy_Art)
	fmt.Printf("  FASTA Indexer:\t%s\n", version_control.FASTA_Indexer)
	fmt.Printf("  ORF to FAA:\t\t%s\n", version_control.ORF_to_FAA)
	fmt.Printf("  Seq Simulator:\t%s\n", version_control.Seq_Sim)
	fmt.Printf("  FASTQC_Mimic:\t\t%s\n", version_control.FastQC_Mimic)
	fmt.Printf("  FASTA_Isolate:\t%s\n", version_control.FASTA_Isolate)
	
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
		case "seq_gen":
			seq_generator.Run(cleanedArgs)
		case "check":
			sanity_check.Run(cleanedArgs)
		case "kmer_analyzer":
			kmer_analyzer.Run_kmer_analyzer(cleanedArgs)
		case "orf_finder":
			orf_finder.Run(cleanedArgs)
		case "fasta_overview":
			fasta_overview.Run(cleanedArgs)
		case "lab_buddy_art":
			lab_buddy_art.PrintLabBuddyArt()
		case "index_fasta":
			fasta_indexer.FastaIndex_Run(cleanedArgs)
		case "orf_to_faa":
			orf_to_faa.Orf_to_faa_Run(cleanedArgs)
		case "seq_sim":
			seq_sim.SeqSimRun(cleanedArgs)
		case "fastqc_mimic":
			fastqc_mimic.FASTQCmimic_Run(cleanedArgs)
		case "fasta_isolate":
			fasta_isolate.FastaIsolate_Run(cleanedArgs)
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
