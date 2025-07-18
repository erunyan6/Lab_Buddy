package fastqc_mimic

import (
	"flag"
	"fmt"
	"os"
	"sync"
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
		err := WritePerReadCSVConcurrent(*outFile, records)
		if err != nil {
			fmt.Println("Failed to write per-read CSV:", err)
		} else {
			fmt.Printf("Wrote FASTQ per-read statisitcs to CSV file: %s_per_read.csv\n", *outFile)
		}
	}

	if *htmlOut {
		// Gather sample, instead of the entire freaking data 
		sampled := SampleReads(records, 100000)
		
		var (
			svgLength, svgGC, svgPQual, svgRQuality, svgGCBase, svgBaseContent, svgDuplication, svgKmerEnrichment string
		)
		
		var wg sync.WaitGroup
		wg.Add(8) // Number of concurrent graphs
		
		go func() {
			defer wg.Done()
			lengths := make([]float64, len(sampled))
			for i, r := range sampled {
				lengths[i] = float64(len(r.Sequence))
			}
			if s, err := GenerateLengthLinePlotSVG(lengths); err == nil {
				svgLength = s
			} else {
				fmt.Println("Failed to generate Read Length plot:", err)
				svgLength = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			maxLen := stats.MaxLength
			perBaseGC := ComputePerBaseGCContent(sampled, maxLen)
			if s, err := GeneratePerBaseGCPlot(perBaseGC); err == nil {
				svgGCBase = s
			} else {
				fmt.Println("Failed to generate Per Base GC plot:", err)
				svgGCBase = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			if s, err := GenerateGCContentLinePlot(gcValues); err == nil {
				svgGC = s
			} else {
				fmt.Println("Failed to generate GC plot:", err)
				svgGC = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			if s, err := GeneratePerBaseQualityLinePlot(sampled); err == nil {
				svgPQual = s
			} else {
				fmt.Println("Failed to generate Per-Base Quality plot:", err)
				svgPQual = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			means := computeMeanQuals(sampled)
			if s, err := GeneratePerReadQualityLinePlot(means); err == nil {
				svgRQuality = s
			} else {
				fmt.Println("Failed to generate Per-Read Quality plot:", err)
				svgRQuality = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			var maxLen1 int
			if stats.MaxLength > 100 {
				maxLen1 = 100
			} else {
				maxLen1 = stats.MaxLength
			}
			baseContent := ComputePerBaseSequenceContent(sampled, maxLen1)
			if s, err := GeneratePerBaseSeqContentPlot(baseContent, maxLen1); err == nil {
				svgBaseContent = s
			} else {
				fmt.Println("Failed to generate Per Base Sequence Content plot:", err)
				svgBaseContent = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			dupBuckets := ComputeDuplicationLevels(sampled, 200000)
			dupValues := DuplicationBucketsToPlotData(dupBuckets, len(sampled))
			if s, err := GenerateDuplicationLinePlot(dupValues); err == nil {
				svgDuplication = s
			} else {
				fmt.Println("Failed to generate duplication plot:", err)
				svgDuplication = "<p>Graph unavailable</p>"
			}
		}()
		
		go func() {
			defer wg.Done()
			k := 5
			maxReads := 100000
			trueMaxLen := GetMaxReadLength(sampled, maxReads)
			posCov := CountReadsPerPosition(sampled, trueMaxLen)
			kmerCounts, _ := CountKmerPositions(sampled, k, maxReads, trueMaxLen)
			topKmers := GetTopPositionalKmers(kmerCounts, 6)
			kmerTotals := make(map[string]int)
			for k, v := range kmerCounts {
				for _, c := range v {
					kmerTotals[k] += c
				}
			}
			enrich := ComputeKmerEnrichment(kmerCounts, kmerTotals, posCov, topKmers, trueMaxLen)
			if s, err := GenerateKmerEnrichmentPlot(enrich, topKmers); err == nil {
				svgKmerEnrichment = s
			} else {
				fmt.Println("Failed to generate k-mer enrichment plot:", err)
				svgKmerEnrichment = "<p>Graph unavailable</p>"
			}
		}()
		
		wg.Wait()
		
			

			err = WriteHTMLReport(*outFile, stats, svgLength, svgGC, svgPQual, svgRQuality, svgBaseContent, svgDuplication, svgKmerEnrichment, svgGCBase)
			if err != nil {
				fmt.Println("Failed to write HTML:", err)
				os.Exit(1)
			} else {
				fmt.Printf("Wrote HTML file: %s.html\n", *outFile)
			}
		}
	}

