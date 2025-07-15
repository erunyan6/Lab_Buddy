package fastqc_mimic

import (
	"flag"
	"fmt"
	"os"
)

func FASTQCmimic_Run(args []string) {

	fs := flag.NewFlagSet("fastqc_mimic", flag.ExitOnError) 	// Isolated flag set specifically for "fastqc_mimic" subcommand 
 
	inFile := fs.String("in_file", "", "FASTQ file input")		// Input file (FASTA)
	outFile := fs.String("out_file", "fastq_report", "Prefix for HTML report")
	csvOut := fs.Bool("csv_out", false, "Output FASTQ file statistics in csv form")
	perReadOut := fs.Bool("per_read", false, "Output per-read stats to CSV")
	htmlOut := fs.Bool("html", false, "Output FASTQ statistics and graphs to HTML file")

	err := fs.Parse(args)										// Parse inputs 
	if err != nil {
		fmt.Println("Error parsing flags:", err)				// Check for outright input failures
		os.Exit(1)												// e.g., expected int by recieved str
	}

	if len(fs.Args()) > 0 {										// If unparsed arguments remain:
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())	// Flag the error and report it
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}

	if *inFile == "" {
		fmt.Println("Error: in_file is required")
		fs.Usage()
		os.Exit(1)
	}

	if !*csvOut && !*perReadOut && !*htmlOut {
		fmt.Println("Error: No output format is selected")
		os.Exit(1)
	}

	// Functions to run inFile through
	records, err := ParseFastq(*inFile)
	if err != nil {
		fmt.Println("Failed to parse FASTQ:", err)
		os.Exit(1)
	}

	stats := ExtendedStats(records)
	lengths := make([]float64, len(records))
	for i, r := range records {
		lengths[i] = float64(len(r.Sequence))
	}
	var gcValues []float64
	for _, rec := range records {
		gc := calcGCContent(rec.Sequence)
		gcValues = append(gcValues, gc)
	}
	maxLength := 0
	for _, r := range records {
		if len(r.Quality) > maxLength {
			maxLength = len(r.Quality)
		}
	}
	// Initialize a slice of slices to collect quality per position
	perBaseQuals := make([][]float64, maxLength)
	for i := range perBaseQuals {
		perBaseQuals[i] = []float64{}
	}
	// Populate
	for _, r := range records {
		for i, qChar := range r.Quality {
			score := float64(qChar - 33)
			perBaseQuals[i] = append(perBaseQuals[i], score)
		}
	}


	if *csvOut {
		err := WriteCSVReport(*outFile, stats)
		if err != nil {
			fmt.Println("Failed to write CSV:", err)
		} else {
			fmt.Printf("Wrote FASTQ statistics to CSV file: %s.csv\n", *outFile)
		}
	}	

	if *perReadOut {
		err := WritePerReadCSV(*outFile, records)
		if err != nil {
			fmt.Println("Failed to write per-read CSV:", err)
		} else {
			fmt.Printf("Wrote FASTQ per-read statisitcs to CSV file: %s_per_read.csv\n", *outFile)
		}
	}
	
	var svgLength string
	var svgGC string
	var svgPQual string
	if *htmlOut {
		lengths := make([]float64, len(records))
		for i, r := range records {
			lengths[i] = float64(len(r.Sequence))
		}
		svgLength, err = GenerateLengthHistogramSVG(lengths)
		if err != nil {
			fmt.Println("Failed to generate SVG:", err)
			svgLength = "<p>Graph unavailable</p>"
		}
		svgGC, err = GenerateGCHistogramSVG(gcValues)
		if err != nil {
			fmt.Println("Failed to generate GC plot:", err)
			svgGC = "<p>Graph unavailable</p>"
		}
		svgPQual, err = GeneratePerBaseQualityBoxPlot(records)
		if err != nil {
			fmt.Println("Failed to generate Per-Base Quality plot:", err)
			svgGC = "<p>Graph unavailable</p>"
		}
	
		err = WriteHTMLReport(*outFile, stats, svgLength, svgGC, svgPQual)
		if err != nil {
			fmt.Println("Failed to write HTML:", err)
			os.Exit(1)
		} else {
			fmt.Printf("Wrote HTML file: %s.html\n", *outFile)
		}
	}
	
}



func calcGCContent(seq string) float64 {
	if len(seq) == 0 {
		return 0.0
	}
	var gc int
	for _, base := range seq {
		switch base {
		case 'G', 'g', 'C', 'c':
			gc++
		}
	}
	return float64(gc) / float64(len(seq)) * 100
}
