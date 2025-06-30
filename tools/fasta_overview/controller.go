package fasta_overview

import (
	"flag"
	"fmt"
	"os"
	"strings"
	"io"
	"compress/gzip"
)

// openFileOrGzip opens a plain or gzip-compressed FASTA file
func openFileOrGzip(path string) (io.Reader, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	if strings.HasSuffix(path, ".gz") {
		gzReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, err
		}
		return gzReader, nil
	}

	return file, nil
}

func Run(args []string) {
	fs := flag.NewFlagSet("fasta_overview", flag.ExitOnError)
	inFile := fs.String("in_file", "", "Input FASTA file")
	mode := fs.String("mode", "dna", "Input mode: 'dna' or 'protein'")
	idMotif := fs.String("id_motif", "", "Only analyze sequences whose headers contain this substring")
	err := fs.Parse(args)										// Parse inputs 
	if err != nil {
		fmt.Println("Error parsing flags:", err)				// Check for outright input failures
		os.Exit(1)												// E.g., expected int by recieved str
	}

	if len(fs.Args()) > 0 {										// If unparsed arguments remain:
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())	// Flag the error and report it
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}

	if *inFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -in_file is required")
		fs.Usage()
		os.Exit(1)
	}

	switch *mode {
	case "dna":
		reader, err := openFileOrGzip(*inFile)
		if err != nil {
			fmt.Fprintln(os.Stderr, "Failed to open file:", err)
			os.Exit(1)
		}
		report := CheckFastaDNA(reader, *inFile, *idMotif)
		PrintDNAReport(report)
	case "protein":
		reader, err := openFileOrGzip(*inFile)
		if err != nil {
			fmt.Fprintln(os.Stderr, "Failed to open file:", err)
			os.Exit(1)
		}
		report := CheckFastaProtein(reader, *inFile, *idMotif)
		PrintProteinReport(report)
	
	default:
		fmt.Fprintf(os.Stderr, "Unsupported mode: %s\n", *mode)
		os.Exit(1)
	}
}
