package fasta_overview

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"bytes"
	"bufio"
	"unicode"
	"io"
	"compress/gzip"
	"sort"
)


// Define report structure
type FastaCheckReport struct {
	FileName            string
	CanOpen             bool
	HeaderCount         int
	DuplicateHeaders    int
	EmptyHeaders        int
	ShortSequences      int
	SequenceWithNoData  int
	InvalidBaseCounts   map[rune]int
	TotalBases          int
	TotalSequences      int
	UniqueHeaders       map[string]bool
	SequenceLengths     []int
	SequenceIDs         []string
	SequenceIDLengths   map[string]int
	EmptyLineWarnings   int
	SequenceBeforeHeader int
	Warnings            []string
	GCContent          map[string]float64
	NPercentage        map[string]float64
	MeanGCContent   float64
	MeanNPercentage float64
	WrappedSequenceLines    int
	UnwrappedSequenceCount  int
	SequenceLineLengthStats map[int]int
}


// Report Generator
func PrintReport(report FastaCheckReport) {
	fmt.Printf("FASTA Format Check Report: %s\n", report.FileName)
	fmt.Println("------------------------------------------")

	if report.CanOpen {
		fmt.Println("File opened and read successfully")
	} else {
		fmt.Println("!!! Failed to read the file !!!")
	}

	fmt.Printf("Headers found: %d\n", report.HeaderCount)
	fmt.Printf("Total sequences: %d\n", report.TotalSequences)

	if report.DuplicateHeaders > 0 {
		fmt.Printf("Duplicate headers found: %d\n", report.DuplicateHeaders)
	} else {
		fmt.Println("No duplicate headers found")
	}

	if report.EmptyHeaders > 0 {
		fmt.Printf("Empty headers found: %d\n", report.EmptyHeaders)
	} else {
		fmt.Println("All headers contain names")
	}

	if report.SequenceWithNoData > 0 {
		fmt.Printf("Headers with no sequence: %d\n", report.SequenceWithNoData)
	} else {
		fmt.Println("All headers have associated sequences")
	}

	if report.ShortSequences > 0 {
		fmt.Printf("Sequences under 10 bp: %d\n", report.ShortSequences)
	} else {
		fmt.Println("All sequences are at least 10 bp long")
	}

	if report.SequenceBeforeHeader > 0 {
		fmt.Printf("Sequence lines before first header: %d\n", report.SequenceBeforeHeader)
	}

	if report.EmptyLineWarnings > 0 {
		fmt.Printf("Empty lines within sequences: %d\n", report.EmptyLineWarnings)
	}

	if len(report.InvalidBaseCounts) > 0 {
		totalInvalid := 0
		fmt.Println("Invalid bases found:")
		for base, count := range report.InvalidBaseCounts {
			fmt.Printf("  %c: %d\n", base, count)
			totalInvalid += count
		}
		fmt.Printf("Total invalid (non-ATCGN) bases: %d\n", totalInvalid)
	} else {
		fmt.Println("No invalid bases found (all A, T, C, G, N)")
	}

	fmt.Printf("Total bases in all sequences: %d\n", report.TotalBases)

	if len(report.SequenceLengths) > 0 {
		minLen, maxLen, totalLen := report.SequenceLengths[0], report.SequenceLengths[0], 0
		for _, l := range report.SequenceLengths {
			totalLen += l
			if l < minLen {
				minLen = l
			}
			if l > maxLen {
				maxLen = l
			}
		}
		avgLen := float64(totalLen) / float64(len(report.SequenceLengths))
		fmt.Printf("\nSequence length statistics:\n")
		fmt.Printf("  Shortest: %d bp\n", minLen)
		fmt.Printf("  Longest:  %d bp\n", maxLen)
		fmt.Printf("  Average:  %.2f bp\n", avgLen)
	}

	fmt.Printf("\nPer-sequence lengths:\n")
	for _, id := range report.SequenceIDs {
		fmt.Printf("  %s: %d bp\n", id, report.SequenceIDLengths[id])
	}

	fmt.Printf("\nPer-sequence content statistics:\n")
	for _, id := range report.SequenceIDs {
		gc := report.GCContent[id]
		np := report.NPercentage[id]
		fmt.Printf("  %s: GC = %.2f%%, N = %.2f%%\n", id, gc, np)
	}

	fmt.Printf("\nAverage content across all sequences:\n")
	fmt.Printf("  Mean GC content: %.2f%%\n", report.MeanGCContent)
	fmt.Printf("  Mean N content:  %.2f%%\n", report.MeanNPercentage)

	if len(report.GCContent) > 0 {
		minGC, maxGC := 100.0, 0.0
		for _, gc := range report.GCContent {
			if gc < minGC {
				minGC = gc
			}
			if gc > maxGC {
				maxGC = gc
			}
		}
		fmt.Printf("GC content range: %.2f%% - %.2f%%\n", minGC, maxGC)
	}


	fmt.Println("\nLine wrapping:")
	if report.WrappedSequenceLines > 0 {
		fmt.Printf("  %d sequences appear to use line wrapping (multiple lines per sequence)\n", report.WrappedSequenceLines)
	}
	if report.UnwrappedSequenceCount > 0 {
		fmt.Printf("  %d sequences appear to be unwrapped (single line only)\n", report.UnwrappedSequenceCount)
	}

	if len(report.SequenceLineLengthStats) > 0 {
		fmt.Println("  Sequence line lengths observed (excluding headers):")
		keys := make([]int, 0, len(report.SequenceLineLengthStats))
		for k := range report.SequenceLineLengthStats {
			keys = append(keys, k)
		}
		sort.Ints(keys)
		for _, k := range keys {
			fmt.Printf("    %d bp: %d line(s)\n", k, report.SequenceLineLengthStats[k])
		}
	}


}


// Function to check the format of GZipped and plain FASTA files
func checkFastaFormatFromReader(r io.Reader, fileName string) FastaCheckReport {
	scanner := bufio.NewScanner(r)
	report := FastaCheckReport{
		FileName:                 fileName,
		CanOpen:                  true,
		InvalidBaseCounts:        make(map[rune]int),
		UniqueHeaders:            make(map[string]bool),
		SequenceIDLengths:        make(map[string]int),
		GCContent:                make(map[string]float64),
		NPercentage:              make(map[string]float64),
		SequenceLineLengthStats:  make(map[int]int),
	}

	inSequence := false
	lineNum := 0
	sequenceBuffer := strings.Builder{}
	var currentHeader string
	linesInCurrentSequence := 0 // NEW

	validBases := map[rune]bool{
		'A': true, 'T': true, 'C': true, 'G': true, 'N': true,
	}

	for scanner.Scan() {
		lineNum++
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			if inSequence {
				report.EmptyLineWarnings++
			}
			continue
		}

		if strings.HasPrefix(line, ">") {
			if currentHeader != "" {
				length := sequenceBuffer.Len()
				report.SequenceLengths = append(report.SequenceLengths, length)
				report.SequenceIDLengths[currentHeader] = length
				report.SequenceIDs = append(report.SequenceIDs, currentHeader)
				if length == 0 {
					report.SequenceWithNoData++
				} else if length < 10 {
					report.ShortSequences++
				}

				// Line wrapping tracking
				if linesInCurrentSequence == 1 {
					report.UnwrappedSequenceCount++
				} else if linesInCurrentSequence > 1 {
					report.WrappedSequenceLines++
				}

				// GC and N content
				seq := sequenceBuffer.String()
				if length > 0 {
					var gcCount, nCount int
					for _, base := range seq {
						switch unicode.ToUpper(base) {
						case 'G', 'C':
							gcCount++
						case 'N':
							nCount++
						}
					}
					gcPercent := float64(gcCount) / float64(length) * 100
					nPercent := float64(nCount) / float64(length) * 100
					report.GCContent[currentHeader] = gcPercent
					report.NPercentage[currentHeader] = nPercent
				}
			}

			report.HeaderCount++
			sequenceBuffer.Reset()
			inSequence = true
			linesInCurrentSequence = 0 // reset per sequence

			headerParts := strings.Fields(line[1:])
			if len(headerParts) == 0 {
				report.EmptyHeaders++
				currentHeader = fmt.Sprintf("unnamed_%d", lineNum)
			} else {
				currentHeader = headerParts[0]
			}

		// Make headers unique internally by appending a counter if duplicate
		originalHeader := currentHeader
			counter := 1
			for report.UniqueHeaders[currentHeader] {
				currentHeader = fmt.Sprintf("%s_dup%d", originalHeader, counter)
				counter++
			}
			report.UniqueHeaders[currentHeader] = true

			// Track that it was originally a duplicate
			if counter > 1 {
				report.DuplicateHeaders++
			}

		} else {
			if !inSequence {
				report.SequenceBeforeHeader++
			}

			linesInCurrentSequence++ // NEW
			lineLen := len(line)     // NEW
			report.SequenceLineLengthStats[lineLen]++ // NEW

			for _, base := range line {
				r := unicode.ToUpper(base)
				report.TotalBases++
				if !validBases[r] {
					report.InvalidBaseCounts[r]++
				}
			}
			sequenceBuffer.WriteString(line)
		}
	}

	if currentHeader != "" {
		length := sequenceBuffer.Len()
		report.SequenceLengths = append(report.SequenceLengths, length)
		report.SequenceIDLengths[currentHeader] = length
		report.SequenceIDs = append(report.SequenceIDs, currentHeader)
		if length == 0 {
			report.SequenceWithNoData++
		} else if length < 10 {
			report.ShortSequences++
		}

		// Line wrapping for final sequence
		if linesInCurrentSequence == 1 {
			report.UnwrappedSequenceCount++
		} else if linesInCurrentSequence > 1 {
			report.WrappedSequenceLines++
		}

		// GC and N content
		seq := sequenceBuffer.String()
		if length > 0 {
			var gcCount, nCount int
			for _, base := range seq {
				switch unicode.ToUpper(base) {
				case 'G', 'C':
					gcCount++
				case 'N':
					nCount++
				}
			}
			gcPercent := float64(gcCount) / float64(length) * 100
			nPercent := float64(nCount) / float64(length) * 100
			report.GCContent[currentHeader] = gcPercent
			report.NPercentage[currentHeader] = nPercent
		}
	}

	// Mean GC and N percentage
	var totalGC, totalN float64
	count := float64(len(report.SequenceIDs))
	for _, id := range report.SequenceIDs {
		totalGC += report.GCContent[id]
		totalN += report.NPercentage[id]
	}
	if count > 0 {
		report.MeanGCContent = totalGC / count
		report.MeanNPercentage = totalN / count
	}

	if err := scanner.Err(); err != nil {
		report.CanOpen = false
		report.Warnings = append(report.Warnings, "Error reading file: "+err.Error())
	}
	report.TotalSequences = report.HeaderCount

	return report
}

// Function to check FASTA format and return warnings
func check_format(in_file string) FastaCheckReport {
	// Open the file
	file, err := os.Open(in_file)
	if err != nil {
		return FastaCheckReport{
			FileName: in_file,
			CanOpen:  false,
			Warnings: []string{"Failed to open file: " + err.Error()},
		}
	}
	defer file.Close()

	// Check compression level
	file_ext := strings.ToLower(filepath.Ext(in_file))
	switch file_ext {

	case "fasta.gz", "fa.gz":
		expected := []byte{0x1F, 0x8B} // GZIP magic number
		buffer := make([]byte, len(expected))
		_, err := file.Read(buffer)
		if err != nil {
			return FastaCheckReport{
				FileName: in_file,
				CanOpen:  false,
				Warnings: []string{"Error reading file: " + err.Error()},
			}
		}
		if !bytes.Equal(buffer, expected) {
			return FastaCheckReport{
				FileName: in_file,
				CanOpen:  false,
				Warnings: []string{"File extension is .gz but does not have correct gzip header (expected 1F 8B)"},
			}
		}

		_, err = file.Seek(0, io.SeekStart)
		if err != nil {
			return FastaCheckReport{
				FileName: in_file,
				CanOpen:  false,
				Warnings: []string{"Failed to rewind file for gzip reader: " + err.Error()},
			}
		}

		gzReader, err := gzip.NewReader(file)
		if err != nil {
			return FastaCheckReport{
				FileName: in_file,
				CanOpen:  false,
				Warnings: []string{"Failed to create gzip reader: " + err.Error()},
			}
		}
		defer gzReader.Close()
		return checkFastaFormatFromReader(gzReader, in_file)

	case ".3bit":
		return FastaCheckReport{
			FileName: in_file,
			CanOpen:  false,
			Warnings: []string{"3bit encoded files are still under development"},
		}

	case ".fasta", ".fa":
		return checkFastaFormatFromReader(file, in_file)

	default:
		return FastaCheckReport{
			FileName: in_file,
			CanOpen:  false,
			Warnings: []string{"Unsupported file extension: " + file_ext},
		}
	}
}

// Run function to be called from Main
func Run(args []string) {

	fs := flag.NewFlagSet("fasta_overview", flag.ExitOnError)
	in_file := fs.String("in_file", "", "FASTA input file")
	err := fs.Parse(args)
	if err != nil {
		fmt.Println("Error parsing flags:", err)
		os.Exit(1)
	}

	if len(fs.Args()) > 0 {
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}

	if *in_file == "" {
		fmt.Fprintln(os.Stderr, "Error: -in_file is required")
		fs.Usage()
		os.Exit(1)
	}

	// Function to check FASTA format (compression level, unique headers, unexpected characters, line wrapping, etc.)
	report := check_format(*in_file)
	PrintReport(report)

}
