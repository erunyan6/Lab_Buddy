package fasta_overview

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strings"
	"unicode"
)

// Define report structure — eventually move this to common.go if shared with protein_checker.go
type FastaCheckReport struct {
	FileName                 string
	CanOpen                  bool
	HeaderCount              int
	DuplicateHeaders         int
	EmptyHeaders             int
	ShortSequences           int
	SequenceWithNoData       int
	InvalidBaseCounts        map[rune]int
	TotalBases               int
	TotalSequences           int
	UniqueHeaders            map[string]bool
	SequenceLengths          []int
	SequenceIDs              []string
	SequenceIDLengths        map[string]int
	EmptyLineWarnings        int
	SequenceBeforeHeader     int
	Warnings                 []string
	GCContent                map[string]float64
	NPercentage              map[string]float64
	MeanGCContent            float64
	MeanNPercentage          float64
	WrappedSequenceLines     int
	UnwrappedSequenceCount   int
	SequenceLineLengthStats  map[int]int
	FilteredByMotif  string
	SkippedSequences int

}

// Main DNA analysis function
func CheckFastaDNA(r io.Reader, fileName string, idMotif string) FastaCheckReport {
	scanner := bufio.NewScanner(r)
	report := FastaCheckReport{
		FileName:                fileName,
		CanOpen:                 true,
		InvalidBaseCounts:       make(map[rune]int),
		UniqueHeaders:           make(map[string]bool),
		SequenceIDLengths:       make(map[string]int),
		GCContent:               make(map[string]float64),
		NPercentage:             make(map[string]float64),
		SequenceLineLengthStats: make(map[int]int),
	}

	inSequence := false
	lineNum := 0
	sequenceBuffer := strings.Builder{}
	var currentHeader string
	linesInCurrentSequence := 0
	var lineLengths []int

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
				if idMotif == "" || strings.Contains(strings.ToLower(currentHeader), strings.ToLower(idMotif)) {
					finalizeSequence(&report, currentHeader, sequenceBuffer.String(), linesInCurrentSequence, lineLengths)
				} else {
					report.SkippedSequences++
				}
			}			

			report.HeaderCount++
			sequenceBuffer.Reset()
			lineLengths = lineLengths[:0]
			inSequence = true
			linesInCurrentSequence = 0

			headerParts := strings.Fields(line[1:])
			if len(headerParts) == 0 {
				report.EmptyHeaders++
				currentHeader = fmt.Sprintf("unnamed_%d", lineNum)
			} else {
				currentHeader = headerParts[0]
			}

			originalHeader := currentHeader
			counter := 1
			for report.UniqueHeaders[currentHeader] {
				currentHeader = fmt.Sprintf("%s_dup%d", originalHeader, counter)
				counter++
			}
			report.UniqueHeaders[currentHeader] = true
			if counter > 1 {
				report.DuplicateHeaders++
			}
		} else {
			if !inSequence {
				report.SequenceBeforeHeader++
			}
			
			linesInCurrentSequence++
			lineLen := len(line)
			lineLengths = append(lineLengths, lineLen)
			sequenceBuffer.WriteString(line)			
		}
	}

	if currentHeader != "" {
		if idMotif == "" || strings.Contains(strings.ToLower(currentHeader), strings.ToLower(idMotif)) {
			finalizeSequence(&report, currentHeader, sequenceBuffer.String(), linesInCurrentSequence, lineLengths)
		} else {
			report.SkippedSequences++
		}
	}
	

	if err := scanner.Err(); err != nil {
		report.CanOpen = false
		report.Warnings = append(report.Warnings, "Error reading file: "+err.Error())
	}
	report.FilteredByMotif = idMotif
	report.TotalSequences = len(report.SequenceIDs)

	return report
}

func finalizeSequence(report *FastaCheckReport, header, sequence string, lines int, lineLengths []int) {
	length := len(sequence)
	report.SequenceLengths = append(report.SequenceLengths, length)
	report.SequenceIDLengths[header] = length
	report.SequenceIDs = append(report.SequenceIDs, header)

	report.TotalBases += length

	if length == 0 {
		report.SequenceWithNoData++
	} else if length < 10 {
		report.ShortSequences++
	}

	if lines == 1 {
		report.UnwrappedSequenceCount++
	} else if lines > 1 {
		report.WrappedSequenceLines++
	}

	validBases := map[rune]bool{'A': true, 'T': true, 'C': true, 'G': true, 'N': true}

	for _, l := range lineLengths {
		report.SequenceLineLengthStats[l]++
	}	

	// GC and N content
	var gcCount, nCount int
	for _, base := range sequence {
		upper := unicode.ToUpper(base)
		switch upper {
		case 'G', 'C':
			gcCount++
		case 'N':
			nCount++
		}
		if !validBases[upper] {
			report.InvalidBaseCounts[upper]++
		}
	}
	if length > 0 {
		report.GCContent[header] = float64(gcCount) / float64(length) * 100
		report.NPercentage[header] = float64(nCount) / float64(length) * 100
	}

	// Update means
	var totalGC, totalN float64
	for _, id := range report.SequenceIDs {
		totalGC += report.GCContent[id]
		totalN += report.NPercentage[id]
	}
	count := float64(len(report.SequenceIDs))
	if count > 0 {
		report.MeanGCContent = totalGC / count
		report.MeanNPercentage = totalN / count
	}
}

// Report Printer
func PrintDNAReport(report FastaCheckReport) {
	fmt.Printf("FASTA Format Check Report: %s\n", report.FileName)
	fmt.Println(strings.Repeat("-", 40))

	if !report.CanOpen {
		fmt.Println("!!! Failed to read the file !!!")
		for _, w := range report.Warnings {
			fmt.Println("  -", w)
		}
		return
	}

	if report.FilteredByMotif != "" {
		fmt.Printf("Motif filter applied: only analyzing sequences containing \"%s\"\n", report.FilteredByMotif)
		fmt.Printf("Sequences skipped due to filter: %d\n", report.SkippedSequences)
	} else {
		fmt.Println("No motif filter applied; all sequences analyzed")
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
