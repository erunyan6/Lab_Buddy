package fastqc_mimic

import (
	"math"
	"sync"
	"runtime"
)

type PerReadStat struct {
	Length           int
	GC, N            int
	HomopolymerMax   int
	Entropy          float64
	MeanQual         float64
	BaseCounts       map[rune]int
	Sequence         string
	ReadsWithN       bool
	LowQuality       bool
	Q20Bases, Q30Bases int
}

func ExtendedStats(records []FastqRecord) FastqStats {
	totalReads := len(records)
	if totalReads == 0 {
		return FastqStats{}
	}

	numWorkers := runtime.NumCPU()
	recordChan := make(chan FastqRecord, numWorkers*2)
	statChan := make(chan PerReadStat, numWorkers*2)

	var wg sync.WaitGroup

	// Worker pool
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for rec := range recordChan {
				statChan <- analyzeRecord(rec)
			}
		}()
	}

	// Feed records
	go func() {
		for _, rec := range records {
			recordChan <- rec
		}
		close(recordChan)
		wg.Wait()
		close(statChan)
	}()

	// Aggregate results
	return aggregateStats(statChan, totalReads)
}

func analyzeRecord(rec FastqRecord) PerReadStat {
	seq := rec.Sequence
	qual := rec.Quality
	length := len(seq)

	gc, n := 0, 0
	counts := map[rune]int{}

	for _, b := range seq {
		switch b {
		case 'A', 'a':
			counts['A']++
		case 'T', 't':
			counts['T']++
		case 'C', 'c':
			counts['C']++
			gc++
		case 'G', 'g':
			counts['G']++
			gc++
		case 'N', 'n':
			counts['N']++
			n++
		}
	}

	sumQual := 0
	q20, q30 := 0, 0
	for _, q := range qual {
		score := int(q) - 33
		sumQual += score
		if score >= 20 {
			q20++
		}
		if score >= 30 {
			q30++
		}
	}
	meanQual := float64(sumQual) / float64(length)

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
	}

	return PerReadStat{
		Length:           length,
		GC:               gc,
		N:                n,
		HomopolymerMax:   maxRun,
		Entropy:          shannonEntropy(counts, length),
		MeanQual:         meanQual,
		BaseCounts:       counts,
		Sequence:         seq,
		ReadsWithN:       n > 0,
		LowQuality:       meanQual < 20.0,
		Q20Bases:         q20,
		Q30Bases:         q30,
	}
}

func aggregateStats(statsChan <-chan PerReadStat, totalReads int) FastqStats {
	var (
		totalLen, totalGC, totalN, totalQ20, totalQ30 int
		minLen, maxLen = math.MaxInt32, 0
		homopolymerTotals, maxHomopolymer            = 0, 0
		lowQualReads, readsWithN, duplicateReads     = 0, 0, 0
		entropySum                                    float64
		qualMeans, gcPerRead, lengths                 []float64
		baseCounts                                    = map[rune]int{}
		sequenceHashes                                = map[string]int{}
	)

	for stat := range statsChan {
		lengths = append(lengths, float64(stat.Length))
		totalLen += stat.Length
		totalGC += stat.GC
		totalN += stat.N
		totalQ20 += stat.Q20Bases
		totalQ30 += stat.Q30Bases
		entropySum += stat.Entropy
		qualMeans = append(qualMeans, stat.MeanQual)
		gcPerRead = append(gcPerRead, float64(stat.GC)/float64(stat.Length)*100)

		if stat.Length < minLen {
			minLen = stat.Length
		}
		if stat.Length > maxLen {
			maxLen = stat.Length
		}
		if stat.HomopolymerMax > maxHomopolymer {
			maxHomopolymer = stat.HomopolymerMax
		}
		homopolymerTotals += stat.HomopolymerMax

		if stat.LowQuality {
			lowQualReads++
		}
		if stat.ReadsWithN {
			readsWithN++
		}

		for base, count := range stat.BaseCounts {
			baseCounts[base] += count
		}
		sequenceHashes[stat.Sequence]++
	}

	for _, count := range sequenceHashes {
		if count > 1 {
			duplicateReads += count
		}
	}

	totalATCG := baseCounts['A'] + baseCounts['T'] + baseCounts['C'] + baseCounts['G']

	return FastqStats{
		TotalReads:             totalReads,
		AvgLength:              float64(totalLen) / float64(totalReads),
		MinLength:              minLen,
		MaxLength:              maxLen,
		LengthStdDev:           stddevFloat(lengths),
		GCContent:              percent(totalGC, totalLen),
		GCStdDev:               stddevFloat(gcPerRead),
		NContent:               percent(totalN, totalLen),
		MeanQual:               mean(qualMeans),
		StdQual:                stddevFloat(qualMeans),
		MaxHomopolymer:         maxHomopolymer,
		ReadsWithNPercent:      percent(readsWithN, totalReads),
		LowQualityReadPercent:  percent(lowQualReads, totalReads),
		Q20BasePercent:         percent(totalQ20, totalLen),
		Q30BasePercent:         percent(totalQ30, totalLen),
		AvgAContent:            percent(baseCounts['A'], totalATCG),
		AvgTContent:            percent(baseCounts['T'], totalATCG),
		AvgCContent:            percent(baseCounts['C'], totalATCG),
		AvgGContent:            percent(baseCounts['G'], totalATCG),
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
