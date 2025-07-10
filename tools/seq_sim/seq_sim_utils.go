package seq_sim

import (
	"bufio"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"math"
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


func injectSequencingErrors(
	seq []byte,
	subRate, indelRate, ambigRate float64,
	clusterBias, gcBoost float64,
	maxIndelLen int,
	homoMult float64,
) ([]byte, []bool, []string) {

	var result []byte
	var errorMask []bool
	var mutationLog []string

	window := 7
	lastError := false
	prev := byte(0)
	homoLen := 1

	for i := 0; i < len(seq); i++ {
		b := seq[i]

		// Homopolymer tracking
		if b == prev {
			homoLen++
		} else {
			homoLen = 1
		}
		prev = b

		// Local GC window boost
		start := max(0, i-window)
		end := min(len(seq), i+window+1)
		gcCount := 0
		for j := start; j < end; j++ {
			if seq[j] == 'G' || seq[j] == 'C' || seq[j] == 'g' || seq[j] == 'c' {
				gcCount++
			}
		}
		gcFrac := float64(gcCount) / float64(end-start)

		localSubRate := subRate
		localIndelRate := indelRate

		// GC boost
		if gcFrac > 0.6 {
			localSubRate *= gcBoost
		}

		// Homopolymer indel boost
		if homoLen >= 3 {
			localIndelRate *= homoMult
		}

		// Error momentum boost
		if lastError {
			localSubRate *= clusterBias
			localIndelRate *= clusterBias
		}

		// Ambiguous base
		if ambigRate > 0 && rand.Float64() < ambigRate {
			result = append(result, 'N')
			errorMask = append(errorMask, true)
			mutationLog = append(mutationLog, fmt.Sprintf("%c → N @%d", b, i))
			lastError = true
			continue
		}

		// Substitution
		if localSubRate > 0 && rand.Float64() < localSubRate {
			mut := randBase(b)
			result = append(result, mut)
			errorMask = append(errorMask, true)
			mutationLog = append(mutationLog, fmt.Sprintf("%c → %c @%d", b, mut, i))
			lastError = true
			continue
		}

		// Indels
		if localIndelRate > 0 {
			r := rand.Float64()
			if r < localIndelRate/2 {
				// Deletion
				delLen := min(maxIndelLen, len(seq)-i)
				mutationLog = append(mutationLog, fmt.Sprintf("del @%d: %s", i, seq[i:i+delLen]))
				lastError = true
				i += delLen - 1 // skip ahead
				continue
			} else if r < localIndelRate {
				// Insertion
				insLen := 1 + rand.Intn(maxIndelLen)
				inserted := make([]byte, insLen)
				for j := range inserted {
					inserted[j] = randBase(0)
				}
				result = append(result, inserted...)
				for j := 0; j < insLen; j++ {
					errorMask = append(errorMask, true)
				}
				mutationLog = append(mutationLog, fmt.Sprintf("ins @%d: %s", i, inserted))
				lastError = true
			}
		}

		// Normal base
		result = append(result, b)
		errorMask = append(errorMask, false)
		lastError = false
	}

	return result, errorMask, mutationLog
}


func randBase(exclude byte) byte {
	bases := []byte{'A', 'C', 'G', 'T'}
	for {
		b := bases[rand.Intn(4)]
		if b != exclude {
			return b
		}
	}
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func randReadLen(mean, stddev, min, max int) int {
	if stddev == 0 {
		return mean
	}
	for {
		// Draw from normal distribution using Box-Muller transform
		u1 := rand.Float64()
		u2 := rand.Float64()
		n := math.Sqrt(-2*math.Log(u1)) * math.Cos(2*math.Pi*u2)
		length := int(n*float64(stddev)) + mean
		if length >= min && length <= max {
			return length
		}
	}
}


func simulateRegion(
	fasta_file string,
	index_map map[string]IndexRecord,
	fasta_header string,
	start int,
	end int,
	readLenMean, readLenStdDev, readLenMin, readLenMax int,
	coverageDepth int,
	writer io.Writer,
	errorRate, indelRate, ambigRate float64,
	qualityProfile string, logErrors bool,
	clusterBias, gcBoost float64,
	maxIndelLen int,
	homopolymerMultiplier float64,
) error {

	// Open FASTA file
	f, err := os.Open(fasta_file)
	if err != nil {
		return fmt.Errorf("failed to open fasta file: %w", err)
	}
	defer f.Close()

	// Allocate reusable buffers
	maxReadLen := readLenMax
	buf := make([]byte, maxReadLen*2)
	qualityBuf := make([]byte, maxReadLen)
	for i := range qualityBuf {
		qualityBuf[i] = 'I'
	}

	// Get index record
	rec, ok := index_map[fasta_header]
	if !ok {
		return fmt.Errorf("fasta header %q not found in index", fasta_header)
	}

	regionLen := end - start
	if regionLen < readLenMin {
		return fmt.Errorf("region %s:%d-%d too short for minimum read length %d", fasta_header, start, end, readLenMin)
	}

	// Simulate reads until target coverage is reached
	targetBases := regionLen * coverageDepth
	basesSimulated := 0

	for basesSimulated < targetBases {
		readLen := randReadLen(readLenMean, readLenStdDev, readLenMin, readLenMax)

		if regionLen < readLen {
			continue // skip if region is too short for this read
		}

		baseStart := rand.Intn(regionLen - readLen + 1) + start
		baseEnd := baseStart + readLen

		byteStart := calcByteOffset(baseStart, rec)
		byteEnd := calcByteOffset(baseEnd, rec)

		rawSeq, err := extractSequence(f, byteStart, byteEnd, buf)
		if err != nil {
			return fmt.Errorf("failed extracting read at %d-%d: %w", baseStart, baseEnd, err)
		}

		// Strand flip
		strand := "+"
		if rand.Float64() < 0.5 {
			rawSeq = reverseComplementBytes(rawSeq)
			strand = "-"
		}

		// Inject sequencing errors
		originalSeq := make([]byte, len(rawSeq))
		copy(originalSeq, rawSeq)
		
		readID := fmt.Sprintf("@%s_%d_%d_(%s)", fasta_header, baseStart, baseEnd, strand)
		
		// Now inject errors and collect errorMask + mutation log
		mutatedSeq, errorMask, mutationLog := injectSequencingErrors(
			rawSeq,
			errorRate,
			indelRate,
			ambigRate,
			clusterBias,
			gcBoost,
			maxIndelLen,
			homopolymerMultiplier,
		)
		
		if logErrors {
			for _, entry := range mutationLog {
				fmt.Fprintf(os.Stderr, "%s MUT %s\n", readID, entry)
			}
		}
		
		
		var qual []byte
		switch strings.ToLower(qualityProfile) {
		case "short":
			qual = generateShortReadQual(mutatedSeq, errorMask)
		case "long":
			qual = generateLongReadQual(mutatedSeq, errorMask)
		default:
			return fmt.Errorf("invalid quality_profile: %s (choose 'short' or 'long')", qualityProfile)
		}

		// Write FASTQ
		fmt.Fprintf(writer, "%s\n%s\n+\n%s\n", readID, mutatedSeq, qual)
		basesSimulated += readLen
	}

	return nil
}

func simulateRegionPaired(
	fasta_file string,
	index_map map[string]IndexRecord,
	fasta_header string,
	start int,
	end int,
	fragLenMean, fragLenStdDev int,
	readLenMin, readLenMax int,
	coverageDepth int,
	writer1, writer2 io.Writer,
	errorRate, indelRate, ambigRate float64,
	qualityProfile string, logErrors bool,
	clusterBias, gcBoost float64,
	maxIndelLen int,
	homopolymerMultiplier float64,
) error {
	// Open FASTA file
	f, err := os.Open(fasta_file)
	if err != nil {
		return fmt.Errorf("failed to open fasta file: %w", err)
	}
	defer f.Close()

	// Index record
	rec, ok := index_map[fasta_header]
	if !ok {
		return fmt.Errorf("fasta header %q not found in index", fasta_header)
	}

	regionLen := end - start
	if regionLen < readLenMin*2 {
		return fmt.Errorf("region %s:%d-%d too short for paired-end reads", fasta_header, start, end)
	}

	// Simulate to meet target coverage
	targetBases := regionLen * coverageDepth
	basesSimulated := 0

	fragLenMax := fragLenMean + 3*fragLenStdDev
	bufferSize := fragLenMax + readLenMax
	
	buf := make([]byte, bufferSize)
	
	for basesSimulated < targetBases {
		fragLen := randReadLen(fragLenMean, fragLenStdDev, readLenMin*2, readLenMax*2)
		if regionLen < fragLen {
			continue
		}
		fragStart := rand.Intn(regionLen-fragLen+1) + start
		fragEnd := fragStart + fragLen

		byteStart := calcByteOffset(fragStart, rec)
		byteEnd := calcByteOffset(fragEnd, rec)

		fragSeq, err := extractSequence(f, byteStart, byteEnd, buf)
		if err != nil {
			return fmt.Errorf("failed extracting fragment %d-%d: %w", fragStart, fragEnd, err)
		}

		// First read: forward from fragStart
		read1Seq := fragSeq[:readLenMin]
		read2Seq := reverseComplementBytes(fragSeq[len(fragSeq)-readLenMin:])

		readIDBase := fmt.Sprintf("@%s_%d_%d", fasta_header, fragStart, fragEnd)

		// Apply sequencing errors
		r1Mut, r1Mask, r1Log := injectSequencingErrors(
			read1Seq, errorRate, indelRate, ambigRate,
			clusterBias, gcBoost, maxIndelLen, homopolymerMultiplier,
		)
		r2Mut, r2Mask, r2Log := injectSequencingErrors(
			read2Seq, errorRate, indelRate, ambigRate,
			clusterBias, gcBoost, maxIndelLen, homopolymerMultiplier,
		)

		if logErrors {
			for _, entry := range r1Log {
				fmt.Fprintf(os.Stderr, "%s/1 MUT %s\n", readIDBase, entry)
			}
			for _, entry := range r2Log {
				fmt.Fprintf(os.Stderr, "%s/2 MUT %s\n", readIDBase, entry)
			}
		}

		var qual1, qual2 []byte
		switch strings.ToLower(qualityProfile) {
		case "short":
			qual1 = generateShortReadQual(r1Mut, r1Mask)
			qual2 = generateShortReadQual(r2Mut, r2Mask)
		case "long":
			qual1 = generateLongReadQual(r1Mut, r1Mask)
			qual2 = generateLongReadQual(r2Mut, r2Mask)
		default:
			return fmt.Errorf("invalid quality_profile: %s", qualityProfile)
		}

		// Write output
		r1ID := readIDBase + "/1"
		r2ID := readIDBase + "/2"

		if writer1 == writer2 {
			fmt.Fprintf(writer1, "%s\n%s\n+\n%s\n", r1ID, r1Mut, qual1)
			fmt.Fprintf(writer2, "%s\n%s\n+\n%s\n", r2ID, r2Mut, qual2)
		} else {
			fmt.Fprintf(writer1, "%s\n%s\n+\n%s\n", r1ID, r1Mut, qual1)
			fmt.Fprintf(writer2, "%s\n%s\n+\n%s\n", r2ID, r2Mut, qual2)
		}

		basesSimulated += fragLen
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


func generateShortReadQual(seq []byte, errorMask []bool) []byte {
	q := make([]byte, len(seq))
	readLen := len(seq)

	for i := 0; i < readLen; i++ {
		if errorMask[i] {
			q[i] = byte(33 + 10 + rand.Intn(6)) // Q10–Q15 for errors
			continue
		}

		pos := float64(i)
		length := float64(readLen)

		// Mild rise from Q30 to Q40 over the first 20 bases
		if pos < 20 {
			score := 30.0 + (10.0 * pos / 20.0) // Q30 → Q40
			q[i] = byte(33 + int(score))
		} else if pos < 50 {
			q[i] = byte(33 + 40) // plateau at Q40
		} else {
			// Gentle decay from Q40 → Q30
			score := 40.0 - ((pos - 50.0) / (length - 50.0) * 10.0)
			if score < 30.0 {
				score = 30.0
			}
			q[i] = byte(33 + int(score))
		}
	}
	return q
}


func generateLongReadQual(seq []byte, errorMask []bool) []byte {
	q := make([]byte, len(seq))

	for i := 0; i < len(seq); i++ {
		if errorMask[i] {
			q[i] = byte(33 + 7 + rand.Intn(4)) // Q7–Q10
			continue
		}

		// Simulate ONT bumpiness
		baseQ := 10 + rand.Intn(10) // Q10–Q20
		if rand.Float64() < 0.02 {
			baseQ -= rand.Intn(6) // occasional dip
		}
		if baseQ < 5 {
			baseQ = 5
		}
		q[i] = byte(33 + baseQ)
	}

	return q
}


func homopolymerLength(seq []byte, pos int) int {
	base := seq[pos]
	length := 1

	// Look backward
	for i := pos - 1; i >= 0 && seq[i] == base; i-- {
		length++
	}

	// Look forward
	for i := pos + 1; i < len(seq) && seq[i] == base; i++ {
		length++
	}

	return length
}


