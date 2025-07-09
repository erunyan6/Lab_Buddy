package seq_sim

import (
	"bufio"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

type IndexRecord struct {
	SeqID        string
	SeqLen       int
	Offset       int64
	BasesPerLine int
	BytesPerLine int
}

// Read FASTA Index into map			map[string]IndexRecord
func parse_fai(index_file string) (map[string]IndexRecord, error) {
	index := make(map[string]IndexRecord)

	f, err := os.Open(index_file)
	if err != nil {
		return nil, fmt.Errorf("failed to open index file: %w", err)
	}
	defer f.Close() // <-- add this

	scanner := bufio.NewScanner(f)
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

		index[seqID] = IndexRecord{
			SeqID:        seqID,
			SeqLen:       seqLen,
			Offset:       offset,
			BasesPerLine: basesPerLine,
			BytesPerLine: bytesPerLine,
		}
	}

	// Scanner error check
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error reading index file: %w", err)
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

////////////////////////////////////////////


func calcByteOffset(basePos int, rec IndexRecord) int64 {
	lineCount := basePos / rec.BasesPerLine
	extraBytes := lineCount * (rec.BytesPerLine - rec.BasesPerLine)
	return rec.Offset + int64(basePos) + int64(extraBytes)
}


func simulateRegion(fasta_file string, index_map map[string]IndexRecord, fasta_header string, start int, end int,
	readLen int, coverageDepth int, writer io.Writer) error {

	// Open FASTA file
	f, err := os.Open(fasta_file)
	if err != nil {
		return fmt.Errorf("failed to open fasta file: %w", err)
	}
	defer f.Close()

	// Allocate reusable buffers
	buf := make([]byte, readLen*2)              // Buffer for FASTA slice
	qualityBuf := make([]byte, readLen)         // Reusable quality string
	for i := range qualityBuf {
		qualityBuf[i] = 'I'
	}

	// Get index record
	rec, ok := index_map[fasta_header]
	if !ok {
		return fmt.Errorf("fasta header %q not found in index", fasta_header)
	}

	// Validate region length
	if end-start < readLen {
		return fmt.Errorf("region %s:%d-%d too short for readLen %d", fasta_header, start, end, readLen)
	}

	// Calculate total reads
	regionLen := end - start
	numReads := (regionLen * coverageDepth) / readLen

	for i := 0; i < numReads; i++ {
		baseStart := rand.Intn(regionLen - readLen + 1) + start
		baseEnd := baseStart + readLen

		byteStart := calcByteOffset(baseStart, rec)
		byteEnd := calcByteOffset(baseEnd, rec)

		rawSeq, err := extractSequence(f, byteStart, byteEnd, buf)
		if err != nil {
			return fmt.Errorf("read #%d failed: %w", i, err)
		}

		// Strand flip
		strand := "+"
		if rand.Float64() < 0.5 {
			rawSeq = reverseComplementBytes(rawSeq)
			strand = "-"
		}

		readID := fmt.Sprintf("@%s_%d_%d_(%s)", fasta_header, baseStart, baseEnd, strand)

		// Write FASTQ entry
		fmt.Fprintf(writer, "%s\n%s\n+\n%s\n", readID, rawSeq, qualityBuf[:len(rawSeq)])
	}

	return nil
}


func extractSequence(f *os.File, byteStart, byteEnd int64, buf []byte) ([]byte, error) {
	readLen := byteEnd - byteStart
	if int64(cap(buf)) < readLen {
		return nil, fmt.Errorf("buffer too small for read")
	}

	_, err := f.Seek(byteStart, io.SeekStart)
	if err != nil {
		return nil, fmt.Errorf("seek failed: %w", err)
	}

	n, err := f.Read(buf[:readLen])
	if err != nil && err != io.EOF {
		return nil, fmt.Errorf("read failed: %w", err)
	}

	// In-place filtering: remove \n and \r
	clean := buf[:0] // reuse buf but reset length
	for _, b := range buf[:n] {
		if b != '\n' && b != '\r' {
			clean = append(clean, b)
		}
	}

	return clean, nil
}

func reverseComplementBytes(seq []byte) []byte {
	// Allocate a new slice to hold the result
	rc := make([]byte, len(seq))
	last := len(seq) - 1

	for i, b := range seq {
		rc[last-i] = complement(b)
	}

	return rc
}

func complement(b byte) byte {
	switch b {
	case 'A', 'a':
		return 'T'
	case 'T', 't':
		return 'A'
	case 'C', 'c':
		return 'G'
	case 'G', 'g':
		return 'C'
	case 'N', 'n':
		return 'N'
	default:
		return 'N' // fallback for ambiguous/invalid bases
	}
}