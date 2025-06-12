package main

import (
    "fmt"
    "os"

	"lab_buddy_go/config"

    "lab_buddy_go/example"
	"lab_buddy_go/ran_dna_gen"
    
)

// Main controller 
func main() {
    opts := config.ParseArgs(os.Args[1:])	// Parse CLI options

    // Use selected tool
    switch opts.Tool {
    case "check":
        example.Run(opts.Params)  // Sanity check 
	case "ran_dna_gen":
		ran_dna_gen.Run(opts.Params)	// Generate random DNA (FASTA format)
		fmt.Print("Parameters: ", opts.Params, "\n")
    default:
        fmt.Printf("Unknown tool: %s\n", opts.Tool)		// Default error message
    }
}
