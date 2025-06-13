package main

import (
    "fmt"
    "os"

    "lab_buddy_go/sanity_check"
	"lab_buddy_go/ran_dna_gen"
    
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
    default: 
        fmt.Printf("Unknown tool: %s\n", toolName)
    }
    
}
