package fastqc_mimic

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"crypto/md5"
	"encoding/hex"
	"strings"
)

type FastqStats struct {
	TotalReads             int
	AvgLength              float64
	MinLength              int
	MaxLength              int
	LengthStdDev           float64
	GCContent              float64
	GCStdDev               float64
	NContent               float64
	MeanQual               float64
	StdQual                float64
	MaxHomopolymer         int
	ReadsWithNPercent      float64
	LowQualityReadPercent  float64
	Q20BasePercent         float64
	Q30BasePercent         float64
	AvgAContent            float64
	AvgTContent            float64
	AvgCContent            float64
	AvgGContent            float64
	MeanHomopolymer        float64
	ApproxDuplicatePercent float64
	MeanEntropy            float64
}

func WriteCSVReport(filename string, stats FastqStats) error {
	f, err := os.Create(filename + ".csv")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := csv.NewWriter(f)
	defer writer.Flush()

	headers := []string{
		"TotalReads", "AvgLength", "MinLength", "MaxLength", "LengthStdDev",
		"GCContent", "GCStdDev", "NContent", "ReadsWithNPercent", "LowQualityReadPercent",
		"Q20BasePercent", "Q30BasePercent", "MeanQual", "StdQual", "MaxHomopolymer",
		"MeanHomopolymer", "ApproxDuplicatePercent", "MeanEntropy",
		"AvgAContent", "AvgTContent", "AvgCContent", "AvgGContent",
	}

	values := []string{
		strconv.Itoa(stats.TotalReads),
		fmt.Sprintf("%.2f", stats.AvgLength),
		strconv.Itoa(stats.MinLength),
		strconv.Itoa(stats.MaxLength),
		fmt.Sprintf("%.2f", stats.LengthStdDev),
		fmt.Sprintf("%.2f", stats.GCContent),
		fmt.Sprintf("%.2f", stats.GCStdDev),
		fmt.Sprintf("%.2f", stats.NContent),
		fmt.Sprintf("%.2f", stats.ReadsWithNPercent),
		fmt.Sprintf("%.2f", stats.LowQualityReadPercent),
		fmt.Sprintf("%.2f", stats.Q20BasePercent),
		fmt.Sprintf("%.2f", stats.Q30BasePercent),
		fmt.Sprintf("%.2f", stats.MeanQual),
		fmt.Sprintf("%.2f", stats.StdQual),
		strconv.Itoa(stats.MaxHomopolymer),
		fmt.Sprintf("%.2f", stats.MeanHomopolymer),
		fmt.Sprintf("%.2f", stats.ApproxDuplicatePercent),
		fmt.Sprintf("%.2f", stats.MeanEntropy),
		fmt.Sprintf("%.2f", stats.AvgAContent),
		fmt.Sprintf("%.2f", stats.AvgTContent),
		fmt.Sprintf("%.2f", stats.AvgCContent),
		fmt.Sprintf("%.2f", stats.AvgGContent),
	}

	writer.Write(headers)
	writer.Write(values)
	return nil
}




func WritePerReadCSV(filename string, records []FastqRecord) error {
	f, err := os.Create(filename + "_per_read.csv")
	if err != nil {
		return err
	}
	defer f.Close()

	writer := csv.NewWriter(f)
	defer writer.Flush()

	headers := []string{
		"ReadID", "Length", "GCContent", "NCount", "HomopolymerMax",
		"Entropy", "MeanQual", "StdQual", "MinQual", "MaxQual",
		"Q20Bases", "Q30Bases", "GCStart", "GCEnd", "GCDelta",
		"GCSkewStart", "GCSkewEnd", "QualDrop3Prime", "ATSkew", "CGSkew",
		"AmbiguousRatio", "ReadHash", "HasLowComplexity",
	}

	writer.Write(headers)

	for _, rec := range records {
		seq := rec.Sequence
		qual := rec.Quality
		length := len(seq)

		gc, n := 0, 0
		counts := map[rune]int{}
		maxRun, curRun := 0, 0
		prev := rune(-1)

		for _, base := range seq {
			counts[base]++
			switch base {
			case 'G', 'g', 'C', 'c':
				gc++
			case 'N', 'n':
				n++
			}
			if base == prev {
				curRun++
			} else {
				curRun = 1
				prev = base
			}
			if curRun > maxRun {
				maxRun = curRun
			}
		}

		gcContent := percent(gc, length)
		entropy := shannonEntropy(counts, length)
		lowComplexity := entropy < 1.5

		head := seq[:20]
		tail := seq[len(seq)-20:]
		gcStart := calcGC(head)
		gcEnd := calcGC(tail)
		gcDelta := gcEnd - gcStart
		gcsSkew := calcGCSkew(head)
		gceSkew := calcGCSkew(tail)

		atSkew := calcSkew(counts['A'], counts['T'])
		cgSkew := calcSkew(counts['C'], counts['G'])

		var qsum, q20, q30, qmin, qmax int
		qscores := make([]float64, length)
		for i, q := range qual {
			qi := int(q) - 33
			qscores[i] = float64(qi)
			qsum += qi
			if qi >= 20 {
				q20++
			}
			if qi >= 30 {
				q30++
			}
			if i == 0 || qi < qmin {
				qmin = qi
			}
			if i == 0 || qi > qmax {
				qmax = qi
			}
		}
		qmean := float64(qsum) / float64(length)
		qstd := stddevFloat(qscores)
		qualDrop := mean(qscores[:20]) - mean(qscores[length-20:])

		hash := md5.Sum([]byte(seq))
		readHash := hex.EncodeToString(hash[:])

		values := []string{
			rec.Header,
			strconv.Itoa(length),
			fmt.Sprintf("%.4f", gcContent),
			strconv.Itoa(n),
			strconv.Itoa(maxRun),
			fmt.Sprintf("%.3f", entropy),
			fmt.Sprintf("%.2f", qmean),
			fmt.Sprintf("%.2f", qstd),
			strconv.Itoa(qmin),
			strconv.Itoa(qmax),
			fmt.Sprintf("%.2f", percent(q20, length)),
			fmt.Sprintf("%.2f", percent(q30, length)),
			fmt.Sprintf("%.2f", gcStart),
			fmt.Sprintf("%.2f", gcEnd),
			fmt.Sprintf("%.2f", gcDelta),
			fmt.Sprintf("%.2f", gcsSkew),
			fmt.Sprintf("%.2f", gceSkew),
			fmt.Sprintf("%.2f", qualDrop),
			fmt.Sprintf("%.2f", atSkew),
			fmt.Sprintf("%.2f", cgSkew),
			fmt.Sprintf("%.2f", percent(n, length)),
			readHash,
			strings.ToUpper(strconv.FormatBool(lowComplexity)),
		}

		writer.Write(values)
	}

	return nil
}


func calcGC(seq string) float64 {
	gc := 0
	for _, base := range seq {
		if base == 'G' || base == 'g' || base == 'C' || base == 'c' {
			gc++
		}
	}
	return percent(gc, len(seq))
}

func calcGCSkew(seq string) float64 {
	g, c := 0, 0
	for _, base := range seq {
		switch base {
		case 'G', 'g':
			g++
		case 'C', 'c':
			c++
		}
	}
	return calcSkew(g, c)
}

func calcSkew(a, b int) float64 {
	if a+b == 0 {
		return 0
	}
	return float64(a-b) / float64(a+b)
}