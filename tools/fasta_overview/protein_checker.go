package fasta_overview

import (
	"bufio"
	"fmt"
	"io"
	"strings"
	"unicode"
	"sort"
)

// ProteinCheckReport defines structure for protein FASTA statistics
type ProteinCheckReport struct {
	FileName             string
	CanOpen              bool
	HeaderCount          int
	TotalSequences       int
	DuplicateHeaders     int
	EmptyHeaders         int
	SequenceBeforeHeader int
	SequenceLengths      []int
	SequenceIDs          []string
	SequenceIDLengths    map[string]int
	InvalidAminoAcids    map[rune]int
	Warnings             []string
	FilteredByMotif      string
	SkippedSequences     int
	AminoAcidCounts map[rune]int
	TotalResidues   int
	HydrophobicCount int
	HydrophilicCount int
	OtherCount       int
	ChargedPositive  int
	ChargedNegative  int
	MostCommonAA     rune
	LeastCommonAA    rune
	AmbiguousResidues map[rune]int
	MolecularWeights map[string]float64 // per sequence
	MinMolWeight     float64
	MaxMolWeight     float64
	MeanMolWeight    float64
}


// List of valid 1-letter amino acid codes (excluding B, J, O, U, X, Z)
var validAminoAcids = map[rune]bool{
	'A': true, 'C': true, 'D': true, 'E': true, 'F': true,
	'G': true, 'H': true, 'I': true, 'K': true, 'L': true,
	'M': true, 'N': true, 'P': true, 'Q': true, 'R': true,
	'S': true, 'T': true, 'V': true, 'W': true, 'Y': true,
}

var hydrophobic = map[rune]bool{'A': true, 'V': true, 'I': true, 'L': true, 'M': true, 'F': true, 'Y': true, 'W': true}
var hydrophilic = map[rune]bool{'R': true, 'N': true, 'D': true, 'Q': true, 'E': true, 'K': true, 'S': true, 'T': true, 'H': true}
var positiveCharged = map[rune]bool{'R': true, 'H': true, 'K': true}
var negativeCharged = map[rune]bool{'D': true, 'E': true}

var aaWeights = map[rune]float64{
	'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13,
	'F': 165.19, 'G': 75.07,  'H': 155.16, 'I': 131.17,
	'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
	'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
	'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19,
}


// CheckFastaProtein parses and analyzes a protein FASTA file
func CheckFastaProtein(r io.Reader, fileName string, idMotif string) ProteinCheckReport {
	scanner := bufio.NewScanner(r)
	report := ProteinCheckReport{
		FileName:          fileName,
		CanOpen:           true,
		InvalidAminoAcids: make(map[rune]int),
		SequenceIDLengths: make(map[string]int),
		AminoAcidCounts: make(map[rune]int),
		AmbiguousResidues: make(map[rune]int),
		MolecularWeights: make(map[string]float64),

	}

	inSequence := false
	lineNum := 0
	sequenceBuffer := strings.Builder{}
	var currentHeader string
	headerMap := make(map[string]bool)

	for scanner.Scan() {
		lineNum++
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, ">") {
			if currentHeader != "" {
				if idMotif == "" || strings.Contains(strings.ToLower(currentHeader), strings.ToLower(idMotif)) {
					finalizeProteinSequence(&report, currentHeader, sequenceBuffer.String())
				} else {
					report.SkippedSequences++
				}
			}		

			report.HeaderCount++
			sequenceBuffer.Reset()
			inSequence = true

			headerParts := strings.Fields(line[1:])
			if len(headerParts) == 0 {
				report.EmptyHeaders++
				currentHeader = fmt.Sprintf("unnamed_%d", lineNum)
			} else {
				currentHeader = headerParts[0]
			}

			// Make headers unique internally
			original := currentHeader
			counter := 1
			for headerMap[currentHeader] {
				currentHeader = fmt.Sprintf("%s_dup%d", original, counter)
				counter++
			}
			headerMap[currentHeader] = true
			if counter > 1 {
				report.DuplicateHeaders++
			}

		} else {
			if !inSequence {
				report.SequenceBeforeHeader++
			}
			for _, aa := range line {
				upper := unicode.ToUpper(aa)
				if !validAminoAcids[upper] {
					report.InvalidAminoAcids[upper]++
				}
			}		
			sequenceBuffer.WriteString(line)
		}
	}

	if currentHeader != "" {
		if idMotif == "" || strings.Contains(strings.ToLower(currentHeader), strings.ToLower(idMotif)) {
			finalizeProteinSequence(&report, currentHeader, sequenceBuffer.String())
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

func finalizeProteinSequence(report *ProteinCheckReport, header, sequence string) {
	length := len(sequence)
	report.SequenceIDs = append(report.SequenceIDs, header)
	report.SequenceLengths = append(report.SequenceLengths, length)
	report.SequenceIDLengths[header] = length

	ambiguousSet := map[rune]bool{
		'X': true, 'B': true, 'Z': true, 'J': true, 'U': true, 'O': true,
	}

	for _, aa := range sequence {
		upper := unicode.ToUpper(aa)
		if validAminoAcids[upper] {
			report.AminoAcidCounts[upper]++
			report.TotalResidues++
	
			if hydrophobic[upper] {
				report.HydrophobicCount++
			} else if hydrophilic[upper] {
				report.HydrophilicCount++
			} else {
				report.OtherCount++
			}
	
			if positiveCharged[upper] {
				report.ChargedPositive++
			}
			if negativeCharged[upper] {
				report.ChargedNegative++
			}
	
		} else {
			report.InvalidAminoAcids[upper]++
			if ambiguousSet[upper] {
				report.AmbiguousResidues[upper]++
			}
		}
	}

	var weight float64
	for _, aa := range sequence {
		upper := unicode.ToUpper(aa)
		if val, ok := aaWeights[upper]; ok {
			weight += val
		}
	}
	report.MolecularWeights[header] = weight

	min, max, total := 1e9, 0.0, 0.0
	for _, w := range report.MolecularWeights {
		if w < min {
			min = w
		}
		if w > max {
			max = w
		}
		total += w
	}
	report.MinMolWeight = min
	report.MaxMolWeight = max
	report.MeanMolWeight = total / float64(len(report.MolecularWeights))

	var mostCount = -1
	var leastCount = 1<<31 - 1 // max int
	for aa, count := range report.AminoAcidCounts {
		if count > mostCount {
			report.MostCommonAA = aa
			mostCount = count
		}
		if count < leastCount {
			report.LeastCommonAA = aa
			leastCount = count
		}
	}
}

// PrintProteinReport displays protein FASTA results
func PrintProteinReport(report ProteinCheckReport, mode string) {
	fmt.Printf("Protein FASTA Format Report: %s\n", report.FileName)
	fmt.Println(strings.Repeat("-", 40))

	if !report.CanOpen {
		fmt.Println("!!! Failed to read the file !!!")
		for _, w := range report.Warnings {
			fmt.Println("  -", w)
		}
		return
	}

	fmt.Printf("Input mode: %s\n", strings.ToUpper(mode))

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

	if report.SequenceBeforeHeader > 0 {
		fmt.Printf("Sequence lines before first header: %d\n", report.SequenceBeforeHeader)
	} else {
		fmt.Println("No sequence lines appear before the first header")
	}

	fmt.Printf("\nPer-sequence summary:\n")
	for _, id := range report.SequenceIDs {
		length := report.SequenceIDLengths[id]
		weight := report.MolecularWeights[id]
		fmt.Printf("  %s: %d aa\t\t%.2f Da\n", id, length, weight)
	}	

	if report.TotalResidues > 0 {
		fmt.Println("\nAmino acid composition:")
		keys := make([]rune, 0, len(report.AminoAcidCounts))
		for aa := range report.AminoAcidCounts {
			keys = append(keys, aa)
		}
		sort.Slice(keys, func(i, j int) bool {
			return keys[i] < keys[j]
		})
		for _, aa := range keys {
			count := report.AminoAcidCounts[aa]
			percent := float64(count) / float64(report.TotalResidues) * 100
			fmt.Printf("  %c: %4d (%.2f%%)\n", aa, count, percent)
		}
	}

	if report.TotalResidues > 0 {
		fmt.Printf("\nResidue class composition:\n")
		fmt.Printf("  Hydrophobic: %.2f%%\n", float64(report.HydrophobicCount)/float64(report.TotalResidues)*100)
		fmt.Printf("  Hydrophilic: %.2f%%\n", float64(report.HydrophilicCount)/float64(report.TotalResidues)*100)
		fmt.Printf("  Other:       %.2f%%\n", float64(report.OtherCount)/float64(report.TotalResidues)*100)
	
		fmt.Printf("\nCharged residues:\n")
		fmt.Printf("  Basic - Positive (R, H, K): %d\n", report.ChargedPositive)
		fmt.Printf("  Acidic - Negative (D, E):   %d\n", report.ChargedNegative)
	
		fmt.Printf("\nResidue abundance:\n")
		fmt.Printf("  Most common: %c (%.2f%%)\n", report.MostCommonAA,
			float64(report.AminoAcidCounts[report.MostCommonAA])/float64(report.TotalResidues)*100)
		fmt.Printf("  Least common: %c (%.2f%%)\n", report.LeastCommonAA,
			float64(report.AminoAcidCounts[report.LeastCommonAA])/float64(report.TotalResidues)*100)
	}	

	if len(report.AmbiguousResidues) > 0 {
		fmt.Println("\nAmbiguous amino acid codes detected:")
		for aa, count := range report.AmbiguousResidues {
			fmt.Printf("  %c: %d time(s)\n", aa, count)
		}
	} else {
		fmt.Println("\nNo ambiguous amino acid codes detected")
	}	

	if len(report.Warnings) > 0 {
		fmt.Println("\nWarnings:")
		for _, w := range report.Warnings {
			fmt.Println("  -", w)
		}
	} else {
		fmt.Println("\nNo additional warnings")
	}
}
