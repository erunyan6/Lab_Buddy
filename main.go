package main

import (
    "fmt"
    "os"

    "lab_buddy_go/sanity_check"
	"lab_buddy_go/ran_dna_gen"
    "lab_buddy_go/kmer_analyzer"
    "lab_buddy_go/orf_finder"
    "lab_buddy_go/fasta3bit"
    
)

// Main controller 
func main() {
    if len(os.Args) < 2 {
        fmt.Println("Usage: ./lab_buddy <tool> [flags]")
        os.Exit(1)
    }
    
    toolName := os.Args[1]
    toolArgs := os.Args[2:]

    switch toolName {
    case "ran_dna_gen":
        ran_dna_gen.Run(toolArgs)
    case "check":
        sanity_check.Run(toolArgs)
    case "kmer_analyzer":
        kmer_analyzer.Run(toolArgs)
    case "orf_finder":
        orf_finder.Run(toolArgs)
    case "encoder":
        fasta3bit.Run(toolArgs)
    default: 
        fmt.Printf("Unknown tool: %s\n", toolName)
    }
    
}
