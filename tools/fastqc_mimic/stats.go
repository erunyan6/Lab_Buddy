package fastqc_mimic

import (
	"math"
)

func ExtendedStats(records []FastqRecord) FastqStats {
	totalReads := len(records)
	if totalReads == 0 {
		return FastqStats{}
	}

	var (
		totalLen, totalGC, totalN int
		lengths                   []int
		gcPerRead                 []float64
		qualMeans                 []float64
		minLen, maxLen            = math.MaxInt32, 0
		maxHomopolymer            = 0
		homopolymerTotals         = 0
		readsWithN, lowQualReads  = 0, 0
		q20Bases, q30Bases        = 0, 0
		baseCounts                = map[rune]int{}
		entropySum                = 0.0
		sequenceHashes            = map[string]int{}
	)

	for _, rec := range records {
		seq := rec.Sequence
		qual := rec.Quality

		length := len(seq)
		lengths = append(lengths, length)
		totalLen += length

		if length < minLen {
			minLen = length
		}
		if length > maxLen {
			maxLen = length
		}

		// Base counts and GC/N
		gc, n := 0, 0
		counts := map[rune]int{}
		for _, b := range seq {
			switch b {
			case 'A', 'a':
				counts['A']++
				baseCounts['A']++
			case 'T', 't':
				counts['T']++
				baseCounts['T']++
			case 'C', 'c':
				counts['C']++
				baseCounts['C']++
				gc++
			case 'G', 'g':
				counts['G']++
				baseCounts['G']++
				gc++
			case 'N', 'n':
				n++
				baseCounts['N']++
			}
		}
		if n > 0 {
			readsWithN++
		}
		totalGC += gc
		totalN += n
		gcPerRead = append(gcPerRead, float64(gc)/float64(length)*100)

		// Entropy
		entropySum += shannonEntropy(counts, length)

		// Quality scores
		sumQual := 0
		for _, q := range qual {
			score := int(q) - 33
			sumQual += score
			if score >= 20 {
				q20Bases++
			}
			if score >= 30 {
				q30Bases++
			}
		}
		meanQ := float64(sumQual) / float64(length)
		qualMeans = append(qualMeans, meanQ)
		if meanQ < 20.0 {
			lowQualReads++
		}

		// Homopolymer detection
		currentChar := 'X'
		currentRun := 0
		maxRun := 0
		for _, b := range seq {
			if b == currentChar {
				currentRun++
			} else {
				currentChar = b
				currentRun = 1
			}
			if currentRun > maxRun {
				maxRun = currentRun
			}
			if currentRun > maxHomopolymer {
				maxHomopolymer = currentRun
			}
		}
		homopolymerTotals += maxRun

		// Approx duplicate detection (hash sequences)
		sequenceHashes[seq]++
	}

	// Duplicate %
	duplicateReads := 0
	for _, count := range sequenceHashes {
		if count > 1 {
			duplicateReads += count
		}
	}

	avgLen := float64(totalLen) / float64(totalReads)
	gcContent := float64(totalGC) / float64(totalLen) * 100
	nContent := float64(totalN) / float64(totalLen) * 100
	meanQual := mean(qualMeans)
	lengthStd := stddevFloat(intsToFloats(lengths))
	gcStd := stddevFloat(gcPerRead)
	qualStd := stddevFloat(qualMeans)

	// Base percentages
	totalATCG := baseCounts['A'] + baseCounts['T'] + baseCounts['C'] + baseCounts['G']
	avgA := percent(baseCounts['A'], totalATCG)
	avgT := percent(baseCounts['T'], totalATCG)
	avgC := percent(baseCounts['C'], totalATCG)
	avgG := percent(baseCounts['G'], totalATCG)

	return FastqStats{
		TotalReads:             totalReads,
		AvgLength:              avgLen,
		MinLength:              minLen,
		MaxLength:              maxLen,
		LengthStdDev:           lengthStd,
		GCContent:              gcContent,
		GCStdDev:               gcStd,
		NContent:               nContent,
		MeanQual:               meanQual,
		StdQual:                qualStd,
		MaxHomopolymer:         maxHomopolymer,
		ReadsWithNPercent:      percent(readsWithN, totalReads),
		LowQualityReadPercent:  percent(lowQualReads, totalReads),
		Q20BasePercent:         percent(q20Bases, totalLen),
		Q30BasePercent:         percent(q30Bases, totalLen),
		AvgAContent:            avgA,
		AvgTContent:            avgT,
		AvgCContent:            avgC,
		AvgGContent:            avgG,
		MeanHomopolymer:        float64(homopolymerTotals) / float64(totalReads),
		ApproxDuplicatePercent: percent(duplicateReads, totalReads),
		MeanEntropy:            entropySum / float64(totalReads),
	}
}

func shannonEntropy(counts map[rune]int, length int) float64 {
	if length == 0 {
		return 0
	}
	entropy := 0.0
	for _, c := range counts {
		p := float64(c) / float64(length)
		if p > 0 {
			entropy -= p * math.Log2(p)
		}
	}
	return entropy
}

func percent(part, total int) float64 {
	if total == 0 {
		return 0
	}
	return float64(part) / float64(total) * 100
}


func mean(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}
	sum := 0.0
	for _, v := range values {
		sum += v
	}
	return sum / float64(len(values))
}

func stddevFloat(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}
	m := mean(values)
	var sumSq float64
	for _, v := range values {
		diff := v - m
		sumSq += diff * diff
	}
	return math.Sqrt(sumSq / float64(len(values)))
}

func intsToFloats(ints []int) []float64 {
	out := make([]float64, len(ints))
	for i, v := range ints {
		out[i] = float64(v)
	}
	return out
}
