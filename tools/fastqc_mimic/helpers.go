package fastqc_mimic

import (
	"strings"
	"sort"
	"math/rand"
	"sync"
)

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

func computeMeanQuals(records []FastqRecord) []float64 {
	means := make([]float64, 0, len(records))
	for _, rec := range records {
		sum := 0
		for _, q := range rec.Quality { 
			sum += int(q) - 33 
		}
		if len(rec.Quality) == 0 {
			continue
		}		
		if len(rec.Quality) > 0 {
			mean := float64(sum) / float64(len(rec.Quality))
			means = append(means, mean)
		}
	}
	return means
}

func ComputePerBaseSequenceContent(records []FastqRecord, maxLen int) map[rune][]float64 {
	const numWorkers = 8 // You can make this dynamic based on runtime.NumCPU()

	// Final tallies
	counts := map[rune][]int{
		'A': make([]int, maxLen),
		'C': make([]int, maxLen),
		'G': make([]int, maxLen),
		'T': make([]int, maxLen),
		'N': make([]int, maxLen),
	}
	total := make([]int, maxLen)

	// Lock for shared write
	var mu sync.Mutex
	var wg sync.WaitGroup
	chunkSize := (len(records) + numWorkers - 1) / numWorkers

	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if end > len(records) {
			end = len(records)
		}
		wg.Add(1)
		go func(subset []FastqRecord) {
			defer wg.Done()

			// Local tallies
			localCounts := map[rune][]int{
				'A': make([]int, maxLen),
				'C': make([]int, maxLen),
				'G': make([]int, maxLen),
				'T': make([]int, maxLen),
				'N': make([]int, maxLen),
			}
			localTotal := make([]int, maxLen)

			for _, rec := range subset {
				seq := strings.ToUpper(rec.Sequence)
				length := len(seq)
				if length > maxLen {
					length = maxLen
				}
				for i := 0; i < length; i++ {
					base := rune(seq[i])
					switch base {
					case 'A', 'C', 'G', 'T':
						localCounts[base][i]++
					default:
						localCounts['N'][i]++
					}
					localTotal[i]++
				}
			}

			// Merge local → global
			mu.Lock()
			defer mu.Unlock()
			for base, vec := range localCounts {
				for i := range vec {
					counts[base][i] += vec[i]
				}
			}
			for i := range localTotal {
				total[i] += localTotal[i]
			}
		}(records[start:end])
	}

	wg.Wait()

	// Convert to percentages
	result := make(map[rune][]float64)
	for base, vals := range counts {
		result[base] = make([]float64, maxLen)
		for i := 0; i < maxLen; i++ {
			if total[i] > 0 {
				result[base][i] = float64(vals[i]) / float64(total[i]) * 100.0
			}
		}
	}
	return result
}

func ComputeDuplicationLevels(records []FastqRecord, maxReads int) map[int]int {
	counts := make(map[string]int)
	for _, rec := range records {
		counts[rec.Sequence]++
	}
	limit := maxReads
	if len(records) < maxReads {
		limit = len(records)
	}

	for i := 0; i < limit; i++ {
		seq := records[i].Sequence
		counts[seq]++
	}

	// Bucket by duplication level
	dupBuckets := make(map[int]int)
	for _, count := range counts {
		dupBuckets[count]++
	}
	return dupBuckets
}

func CountKmerPositions(records []FastqRecord, k int, maxReads int, trueMaxLen int) (map[string][]int, int) {
	kmerPosCounts := make(map[string][]int)
	limit := maxReads
	if len(records) < maxReads {
		limit = len(records)
	}
	maxPos := trueMaxLen - k + 1

	for i := 0; i < limit; i++ {
		seq := strings.ToUpper(records[i].Sequence)
		for j := 0; j <= len(seq)-k; j++ {
			kmer := seq[j : j+k]

			if _, ok := kmerPosCounts[kmer]; !ok {
				kmerPosCounts[kmer] = make([]int, maxPos)
			}
			if j >= len(kmerPosCounts[kmer]) {
				continue // safety check, shouldn't trigger with correct maxLen
			}
			kmerPosCounts[kmer][j]++
		}
	}

	// After position counts are done
	kmerTotals := make(map[string]int)
	for k, arr := range kmerPosCounts {
		sum := 0
		for _, v := range arr {
			sum += v
		}
		kmerTotals[k] = sum
	}

	return kmerPosCounts, maxPos
}


func GetTopPositionalKmers(kmerCounts map[string][]int, topN int) []string {
	type kv struct {
		Kmer  string
		Count int
	}
	var all []kv
	for k, v := range kmerCounts {
		sum := 0
		for _, c := range v {
			sum += c
		}
		all = append(all, kv{k, sum})
	}
	sort.Slice(all, func(i, j int) bool {
		return all[i].Count > all[j].Count
	})
	top := make([]string, 0, topN)
	for i := 0; i < topN && i < len(all); i++ {
		top = append(top, all[i].Kmer)
	}
	return top
}

func GetMaxReadLength(records []FastqRecord, maxReads int) int {
	maxLen := 0
	limit := maxReads
	if len(records) < maxReads {
		limit = len(records)
	}
	for i := 0; i < limit; i++ {
		if len(records[i].Sequence) > maxLen {
			maxLen = len(records[i].Sequence)
		}
	}
	return maxLen
}

func ComputeKmerEnrichment(
	kmerCounts map[string][]int,
	kmerTotals map[string]int,
	_ []int, // no longer used
	topKmers []string,
	maxLen int,
) map[string][]float64 {

	enrichment := make(map[string][]float64)

	for _, kmer := range topKmers {
		counts := kmerCounts[kmer]
		enrich := make([]float64, maxLen)

		total := float64(kmerTotals[kmer])
		if total == 0 {
			total = 1 // prevent div-by-zero
		}

		for i := 0; i < maxLen && i < len(counts); i++ {
			enrich[i] = (float64(counts[i]) / total) * 100.0
		}
		enrichment[kmer] = enrich
	}

	return enrichment
}




func CountReadsPerPosition(records []FastqRecord, maxLen int) []int {
	counts := make([]int, maxLen)
	for _, rec := range records {
		l := len(rec.Sequence)
		if l > maxLen {
			l = maxLen
		}
		for i := 0; i < l; i++ {
			counts[i]++
		}
	}
	return counts
}

func ComputePerBaseGCContent(records []FastqRecord, maxLen int) []float64 {
	const numWorkers = 8

	gcCounts := make([]int, maxLen)
	totalCounts := make([]int, maxLen)

	var mu sync.Mutex
	var wg sync.WaitGroup
	chunkSize := (len(records) + numWorkers - 1) / numWorkers

	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if end > len(records) {
			end = len(records)
		}
		wg.Add(1)
		go func(subset []FastqRecord) {
			defer wg.Done()

			localGC := make([]int, maxLen)
			localTotal := make([]int, maxLen)

			for _, rec := range subset {
				seq := strings.ToUpper(rec.Sequence)
				readLen := len(seq)
				if readLen > maxLen {
					readLen = maxLen
				}
				for i := 0; i < readLen; i++ {
					base := seq[i]
					if base == 'G' || base == 'C' {
						localGC[i]++
					}
					localTotal[i]++
				}
			}

			mu.Lock()
			defer mu.Unlock()
			for i := 0; i < maxLen; i++ {
				gcCounts[i] += localGC[i]
				totalCounts[i] += localTotal[i]
			}
		}(records[start:end])
	}

	wg.Wait()

	gcPercent := make([]float64, maxLen)
	for i := 0; i < maxLen; i++ {
		if totalCounts[i] > 0 {
			gcPercent[i] = float64(gcCounts[i]) / float64(totalCounts[i]) * 100.0
		}
	}
	return gcPercent
}

// SampleReads randomly selects up to n reads for plotting
func SampleReads(records []FastqRecord, n int) []FastqRecord {
	if len(records) <= n {
		return records
	}
	sampled := make([]FastqRecord, 0, n)
	perm := rand.Perm(len(records))
	for i := 0; i < n; i++ {
		sampled = append(sampled, records[perm[i]])
	}
	return sampled
}
