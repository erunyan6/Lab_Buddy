package fastqc_mimic

import (
	"bytes"
	"image/color"
	"fmt"
	"math"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

type IntegerTicks struct{}

func (IntegerTicks) Ticks(min, max float64) []plot.Tick {
	var ticks []plot.Tick
	for i := int(math.Ceil(min)); i <= int(math.Floor(max)); i++ {
		ticks = append(ticks, plot.Tick{
			Value: float64(i),
			Label: fmt.Sprintf("%d", i),
		})
	}
	return ticks
}

func GenerateLengthLinePlotSVG(lengths []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "Read Length Distribution"
	p.X.Label.Text = "Read Length"
	p.Y.Label.Text = "Read Count"

	p.X.Tick.Marker = IntegerTicks{}
	
	// Bin size setup
	binCount := 100
	minLen := int(lengths[0])
	maxLen := int(lengths[0])
	for _, l := range lengths {
		if int(l) < minLen {
			minLen = int(l)
		}
		if int(l) > maxLen {
			maxLen = int(l)
		}
	}

	binWidth := float64(maxLen-minLen+1) / float64(binCount)
	counts := make([]float64, binCount)

	for _, val := range lengths {
		bin := int((val - float64(minLen)) / binWidth)
		if bin >= binCount {
			bin = binCount - 1
		}
		counts[bin]++
	}

	// Build line plot points
	points := make(plotter.XYs, binCount)
	for i := 0; i < binCount; i++ {
		x := float64(minLen) + binWidth*float64(i) + binWidth/2
		y := counts[i]
		points[i].X = x
		points[i].Y = y
	}

	line, err := plotter.NewLine(points)
	if err != nil {
		return "", err
	}
	line.LineStyle.Color = color.RGBA{R: 50, G: 100, B: 200, A: 255}
	line.LineStyle.Width = vg.Points(2)
	p.Add(line)
	p.Legend.Add("Read Count", line)
	p.Legend.Top = true

	// Write to SVG
	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}


func GenerateGCContentLinePlot(gcValues []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "Per Sequence GC Content"
	p.X.Label.Text = "GC Content (%)"
	p.Y.Label.Text = "Read Count"

	// A. Build observed GC histogram
	binCount := 100
	binWidth := 100.0 / float64(binCount)
	observed := make([]float64, binCount)

	for _, val := range gcValues {
		bin := int(val / binWidth)
		if bin >= binCount {
			bin = binCount - 1
		}
		observed[bin]++
	}

	// B. Compute mean and stddev of GC
	mean := stat.Mean(gcValues, nil)
	stddev := stat.StdDev(gcValues, nil)

	// C. Build expected normal curve (normalized to same total)
	totalReads := float64(len(gcValues))
	expected := make([]float64, binCount)
	normDist := distuv.Normal{Mu: mean, Sigma: stddev}
	scaleFactor := totalReads * binWidth // for normalization to observed scale

	for i := 0; i < binCount; i++ {
		x := binWidth*float64(i) + binWidth/2
		expected[i] = normDist.Prob(x) * scaleFactor
	}

	// D. Convert to line plots
	observedXY := make(plotter.XYs, binCount)
	expectedXY := make(plotter.XYs, binCount)
	for i := 0; i < binCount; i++ {
		x := binWidth*float64(i) + binWidth/2
		observedXY[i].X = x
		observedXY[i].Y = observed[i]
		expectedXY[i].X = x
		expectedXY[i].Y = expected[i]
	}

	// E. Add lines
	obsLine, err := plotter.NewLine(observedXY)
	if err != nil {
		return "", err
	}
	obsLine.Color = color.RGBA{B: 255, A: 255}
	obsLine.Width = vg.Points(2)
	obsLine.Dashes = []vg.Length{}

	expLine, err := plotter.NewLine(expectedXY)
	if err != nil {
		return "", err
	}
	expLine.Color = color.RGBA{R: 255, G: 100, B: 100, A: 255}
	expLine.Width = vg.Points(2)
	expLine.Dashes = []vg.Length{vg.Points(3), vg.Points(3)}

	p.Add(obsLine, expLine)
	p.Legend.Add("Observed", obsLine)
	p.Legend.Add("Modelled Normal", expLine)
	p.Legend.Top = true

	// F. Export SVG
	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}


func GeneratePerBaseGCPlot(gcPercent []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "Per Base GC Content"
	p.X.Label.Text = "Position in Read (bp)"
	p.Y.Label.Text = "GC Content (%)"
	p.Y.Min = 0
	p.Y.Max = 100

	pts := make(plotter.XYs, len(gcPercent))
	for i, val := range gcPercent {
		pts[i].X = float64(i + 1)
		pts[i].Y = val
	}

	line, err := plotter.NewLine(pts)
	if err != nil {
		return "", err
	}
	line.LineStyle.Color = color.RGBA{B: 200, A: 255}
	line.LineStyle.Width = vg.Points(2)
	p.Add(line)

	p.Legend.Add("GC %", line)
	p.Legend.Top = true

	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}



func GeneratePerBaseQualityBoxPlot(reads []FastqRecord) (string, error) {
	p := plot.New()
	p.Title.Text = "Per-Base Quality (Boxplot)"
	p.X.Label.Text = "Base Position"
	p.Y.Label.Text = "Quality Score"
	p.Y.Min = 0
	p.Y.Max = 45

	// Determine max read length
	maxLength := 0
	for _, r := range reads {
		if len(r.Quality) > maxLength {
			maxLength = len(r.Quality)
		}
	}

	p.Add(plotter.NewGrid())

	// Collect quality scores per base position
	perBaseQuals := make([]plotter.Values, maxLength)
	for _, r := range reads {
		for j, qChar := range r.Quality {
			qScore := float64(qChar - 33)
			perBaseQuals[j] = append(perBaseQuals[j], qScore)
		}
	}

	// Add boxplots per position
	for i, values := range perBaseQuals {
		if len(values) == 0 {
			continue
		}
		box, err := plotter.NewBoxPlot(vg.Points(4), float64(i+1), values)
		if err != nil {
			return "", err
		}
		box.FillColor = color.RGBA{R: 255, G: 215, B: 0, A: 255} // gold-yellow like FastQC
		p.Add(box)
	}

	// Optionally overlay a mean line
	meanLine := make(plotter.XYs, maxLength)
	for i, values := range perBaseQuals {
		if len(values) == 0 {
			continue
		}
		sum := 0.0
		for _, v := range values {
			sum += v
		}
		mean := sum / float64(len(values))
		meanLine[i].X = float64(i + 1)
		meanLine[i].Y = mean
	}
	line, err := plotter.NewLine(meanLine)
	if err != nil {
		return "", err
	}
	line.Color = color.RGBA{B: 255, A: 255} // blue line
	line.Width = vg.Points(1)
	p.Add(line)

	// Output SVG
	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg") // consistent sizing
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
	
}



func GeneratePerReadQualityLinePlot(means []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "Per Sequence Quality Scores"
	p.X.Label.Text = "Mean Quality Score"
	p.Y.Label.Text = "Number of Reads"

	// A. Bin the scores
	binCount := 50
	minScore := math.Floor(minFloat64(means))
	maxScore := math.Ceil(maxFloat64(means))
	binWidth := (maxScore - minScore) / float64(binCount)

	counts := make([]float64, binCount)

	for _, val := range means {
		bin := int((val - minScore) / binWidth)
		if bin >= binCount {
			bin = binCount - 1
		}
		counts[bin]++
	}

	// B. Convert to line points
	pts := make(plotter.XYs, binCount)
	for i := 0; i < binCount; i++ {
		x := minScore + binWidth*float64(i) + binWidth/2
		pts[i].X = x
		pts[i].Y = counts[i]
	}

	line, err := plotter.NewLine(pts)
	if err != nil {
		return "", err
	}
	line.LineStyle.Width = vg.Points(2)
	line.LineStyle.Color = color.RGBA{R: 200, G: 100, B: 100, A: 255}
	p.Add(line)
	p.Legend.Add("Mean Read Quality", line)
	p.Legend.Top = true

	// Save as SVG
	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}

func minFloat64(vals []float64) float64 {
	min := vals[0]
	for _, v := range vals {
		if v < min {
			min = v
		}
	}
	return min
}

func maxFloat64(vals []float64) float64 {
	max := vals[0]
	for _, v := range vals {
		if v > max {
			max = v
		}
	}
	return max
}


func GeneratePerBaseSeqContentPlot(data map[rune][]float64, maxLen int) (string, error) {
	p := plot.New()
	p.Title.Text = "Per Base Sequence Content"
	p.X.Label.Text = "Position in Read"
	p.Y.Label.Text = "Base Composition (%)"
	p.Y.Min = 0
	p.Y.Max = 100
	p.Legend.Top = true
	p.Legend.XOffs = -10 // optional offset to pull it in from edge
	p.Legend.YOffs = 0   // no vertical offset needed


	colors := map[rune]color.RGBA{
		'A': {R: 255, A: 255},
		'C': {G: 200, A: 255},
		'G': {B: 255, A: 255},
		'T': {R: 255, G: 165, A: 255},
		'N': {R: 150, G: 150, A: 255},
	}

	for base, percentages := range data {
		pts := make(plotter.XYs, maxLen)
		for i := 0; i < maxLen; i++ {
			pts[i].X = float64(i + 1)
			pts[i].Y = percentages[i]
		}
		line, err := plotter.NewLine(pts)
		if err != nil {
			return "", err
		}
		line.LineStyle.Width = vg.Points(1.3)
		line.LineStyle.Color = colors[base]
		p.Add(line)
		p.Legend.Add(string(base), line)
	}

	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}

func DuplicationBucketsToPlotData(dupBuckets map[int]int, total int) plotter.Values {
	maxBucket := 0
	for level := range dupBuckets {
		if level > maxBucket {
			maxBucket = level
		}
	}

	values := make(plotter.Values, maxBucket)
	for i := 1; i <= maxBucket; i++ {
		count := dupBuckets[i]
		percent := float64(count) / float64(total) * 100.0
		values[i-1] = percent // bucket i at index i-1
	}
	return values
}

func GenerateDuplicationLinePlot(values plotter.Values) (string, error) {
	p := plot.New()
	p.Title.Text = "Sequence Duplication Levels"
	p.X.Label.Text = "Duplication Count"
	p.Y.Label.Text = "Percent of Reads"
	p.Y.Max = 100
	p.X.Tick.Marker = IntegerTicks{}

	pts := make(plotter.XYs, len(values))
	for i := range values {
		pts[i].X = float64(i + 1) // duplication level starts at 1
		pts[i].Y = values[i]
	}

	line, err := plotter.NewLine(pts)
	if err != nil {
		return "", err
	}
	line.LineStyle.Width = vg.Points(2)
	line.LineStyle.Color = color.RGBA{R: 100, G: 180, B: 255, A: 255}
	p.Add(line)
	p.Legend.Add("Duplication %", line)
	p.Legend.Top = true

	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}

func GenerateKmerEnrichmentPlot(enrichment map[string][]float64, topKmers []string) (string, error) {
	p := plot.New()
	p.Title.Text = "Relative enrichment over read length"
	p.X.Label.Text = "Position in Read (bp)"
	p.Y.Label.Text = "Relative Enrichment (%)"
	p.Y.Min = 0
	p.Y.Max = 100
	p.Legend.Top = true
	p.Legend.XOffs = -10

	colors := []color.RGBA{
		{R: 255, G: 0, B: 0, A: 255},    // red
		{G: 200, A: 255},                // green
		{B: 255, A: 255},                // blue
		{R: 255, G: 165, A: 255},        // orange
		{R: 255, B: 255, A: 255},        // pink
		{R: 200, G: 200, A: 255},        // gray
	}

	for i, kmer := range topKmers {
		values := enrichment[kmer]
		pts := make(plotter.XYs, len(values))
		for j, val := range values {
			pts[j].X = float64(j + 1)
			pts[j].Y = val
		}
		line, err := plotter.NewLine(pts)
		if err != nil {
			return "", err
		}
		line.LineStyle.Width = vg.Points(2)
		line.LineStyle.Color = colors[i%len(colors)]
		p.Add(line)
		p.Legend.Add(kmer, line)
	}

	var buf bytes.Buffer
	writer, err := p.WriterTo(10*vg.Inch, 4*vg.Inch, "svg")
	if err != nil {
		return "", err
	}
	_, err = writer.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}


func SmoothCounts(raw []int, window int) []float64 {
	smoothed := make([]float64, len(raw))
	for i := range raw {
		start := max(0, i-window)
		end := min(len(raw), i+window+1)
		sum := 0
		for j := start; j < end; j++ {
			sum += raw[j]
		}
		smoothed[i] = float64(sum) / float64(end-start)
	}
	return smoothed
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
