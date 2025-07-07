package orf_to_faa

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"io"

	"lab_buddy_go/tools/fasta_indexer"
	"lab_buddy_go/utils"
)

type FastaIndex struct {
    SeqID        string
    SeqLen       int
    Offset       int64
    BasesPerLine int
    BytesPerLine int
}

type ORF struct {
    SeqID     string
    Start     int
    End       int
    Strand    string
    UniqueID  string
}

var codonMap = map[string]rune{
	// Phenylalanine
	"TTT": 'F', "TTC": 'F',
	// Leucine
	"TTA": 'L', "TTG": 'L', "CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	// Isoleucine
	"ATT": 'I', "ATC": 'I', "ATA": 'I',
	// Methionine (Start)
	"ATG": 'M',
	// Valine
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	// Serine
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', "AGT": 'S', "AGC": 'S',
	// Proline
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	// Threonine
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	// Alanine
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	// Tyrosine
	"TAT": 'Y', "TAC": 'Y',
	// Histidine
	"CAT": 'H', "CAC": 'H',
	// Glutamine
	"CAA": 'Q', "CAG": 'Q',
	// Asparagine
	"AAT": 'N', "AAC": 'N',
	// Lysine
	"AAA": 'K', "AAG": 'K',
	// Aspartic Acid
	"GAT": 'D', "GAC": 'D',
	// Glutamic Acid
	"GAA": 'E', "GAG": 'E',
	// Cysteine
	"TGT": 'C', "TGC": 'C',
	// Tryptophan
	"TGG": 'W',
	// Arginine
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R', "AGA": 'R', "AGG": 'R',
	// Glycine
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	// Stop codons
	"TAA": '*', "TAG": '*', "TGA": '*',
}

type ProteinResult struct {
	UniqueID string
	SeqID   string
	Start   int
	End     int
	Strand  string
	Protein string
}

func extractAndTranslateORFs(fasta string, index map[string]FastaIndex, orfList []ORF) ([]ProteinResult, error) {
	f, err := os.Open(fasta)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer f.Close()

	if strings.HasSuffix(fasta, "gz") {
		return nil, fmt.Errorf("gzipped files are not supported by this tool. gunzip and try again")
	}

	var results []ProteinResult

	for _, orf := range orfList {
		entry, ok := index[orf.SeqID]
		if !ok {
			return nil, fmt.Errorf("sequence %s not found in index", orf.SeqID)
		}

		lineNum := (orf.Start - 1) / entry.BasesPerLine
		offsetInLine := (orf.Start - 1) % entry.BasesPerLine
		byteOffset := entry.Offset + int64(lineNum*entry.BytesPerLine+offsetInLine)		
		baseCount := orf.End - orf.Start + 1
		linesToRead := (baseCount + entry.BasesPerLine - 1) / entry.BasesPerLine
		bytesToRead := linesToRead * entry.BytesPerLine

		_, err = f.Seek(byteOffset, io.SeekStart)
		if err != nil {
			return nil, fmt.Errorf("failed to seek: %w", err)
		}

		readBuf := make([]byte, bytesToRead)
		_, err = io.ReadFull(f, readBuf)
		if err != nil {
			return nil, fmt.Errorf("failed to read from FASTA: %w", err)
		}

		cleaned := strings.Map(func(r rune) rune {
			if r == '\n' || r == '\r' {
				return -1
			}
			return r
		}, string(readBuf)) 	
		if len(cleaned) > baseCount {
			cleaned = cleaned[:baseCount]
		}			

		if orf.Strand == "-" {
			cleaned = common.ReverseComplement(cleaned)
		}

		var protein []rune
		for i := 0; i+3 <= len(cleaned); i += 3 {
			codon := cleaned[i : i+3]
			aa, ok := codonMap[codon]
			if !ok {
				aa = 'X'
			}
			protein = append(protein, aa)
		}

		results = append(results, ProteinResult{
			UniqueID: orf.UniqueID,
			SeqID:   orf.SeqID,
			Start:   orf.Start,
			End:     orf.End,
			Strand:  orf.Strand,
			Protein: string(protein),
		})
	}

	return results, nil
}

func parseFai(file string) (map[string]FastaIndex, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open index file: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	index := make(map[string]FastaIndex)

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			return nil, fmt.Errorf("invalid .fai line: %q", line)
		}

		seqID := fields[0]
		seqLen, err1 := strconv.Atoi(fields[1])
		offset, err2 := strconv.ParseInt(fields[2], 10, 64)
		basesPerLine, err3 := strconv.Atoi(fields[3])
		bytesPerLine, err4 := strconv.Atoi(fields[4])

		if err := firstError(err1, err2, err3, err4); err != nil {
			return nil, fmt.Errorf("failed parsing line %q: %w", line, err)
		}

		index[seqID] = FastaIndex{
			SeqID:        seqID,
			SeqLen:       seqLen,
			Offset:       offset,
			BasesPerLine: basesPerLine,
			BytesPerLine: bytesPerLine,
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %w", err)
	}

	return index, nil
}

func firstError(errs ...error) error {
	for _, err := range errs {
		if err != nil {
			return err
		}
	}
	return nil
}

func parseGFF3(file string) ([]ORF, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open gff3 file: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var orfs []ORF

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 9 {
			return nil, fmt.Errorf("invalid gff3 line: %q", line)
		}

		if fields[2] != "ORF" {
			continue
		}

		seqID := fields[0]
		start, err := strconv.Atoi(fields[3])
		if err != nil {
			return nil, fmt.Errorf("invalid start position: %w", err)
		}
		end, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, fmt.Errorf("invalid end position: %w", err)
		}
		if start < -1 || end < -1 {
			continue
		}
		directionality := fields[6]
		var uniqueID string
		for _, field := range strings.Split(fields[8], ";") {
			field = strings.TrimSpace(field)
			if strings.HasPrefix(field, "ID=") {
				parts := strings.SplitN(field, "=", 2)
				if len(parts) == 2 {
					uniqueID = parts[1]
					break
				}
			}
		}
		if uniqueID == "" {
			uniqueID = "unknown"
		}

		orfs = append(orfs, ORF{
			SeqID:    seqID,
			Start:    start,
			End:      end,
			Strand:   directionality,
			UniqueID: uniqueID,
		})
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %w", err)
	}

	return orfs, nil
}

func writeFaa(results []ProteinResult, outPath string) error {
	var writer *bufio.Writer
	var file *os.File
	var err error

	if outPath != "" {
		file, err = os.Create(outPath)
		if err != nil {
			return fmt.Errorf("failed to create .faa file: %w", err)
		}
		defer file.Close()
		writer = bufio.NewWriter(file)
		defer writer.Flush()
	} else {
		writer = bufio.NewWriter(os.Stdout)
		defer writer.Flush()
	}

	for _, res := range results {
		fmt.Fprintf(writer, ">%s|%s:%d-%d [%s]\n", res.UniqueID, res.SeqID, res.Start, res.End, res.Strand)

		prot := res.Protein
		lineWidth := 60
		for i := 0; i < len(prot); i += lineWidth {
			end := i + lineWidth
			if end > len(prot) {
				end = len(prot)
			}
			fmt.Fprintln(writer, prot[i:end])
		}
	}

	return nil
}


func Orf_to_faa_Run(args []string) {
	fs := flag.NewFlagSet("orf_to_faa", flag.ExitOnError)
	inputFile := fs.String("in_file", "", "Input FASTA file")
	gffFile := fs.String("orf_file", "", "GFF3 file with ORFs")
	outFile := fs.String("out_file", "", "Output .faa file (default: stdout)")
	fs.Parse(args)

	if *inputFile == "" || *gffFile == "" {
		log.Fatal("Error: -in_file and -orf_file are required")
	}
	if len(fs.Args()) > 0 {
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}

	// Always regenerate the index before proceeding
	fasta_indexer.FastaIndex_Run([]string{"-in_file", *inputFile})
	indexPath := *inputFile + ".fai"

	// Check if the index is fresh
	if err := common.CheckIndexFreshness(*inputFile, indexPath); err != nil {
		log.Fatalf("Index freshness check failed: %v", err)
	}

	// Parse the index
	index, err := parseFai(indexPath)
	if err != nil {
		log.Fatalf("Failed to parse FASTA index: %v", err)
	}

	// Parse the ORF list
	orfs, err := parseGFF3(*gffFile)
	if err != nil {
		log.Fatalf("Failed to parse GFF3: %v", err)
	}

	var results []ProteinResult
	
	// Extract and translate (to be implemented)
	results, err = extractAndTranslateORFs(*inputFile, index, orfs)
	if err != nil {
		log.Fatalf("Translation failed: %v", err)
	}

	err = writeFaa(results, *outFile)
	if err != nil {
		log.Fatalf("Failed to write output: %v", err)
	}

	if *outFile != "" {
		fmt.Printf("Wrote %d proteins to %s\n", len(results), *outFile)
	}
	
}
