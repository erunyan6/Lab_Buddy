package fastqc_mimic

import (
	"bytes"
	"image/color"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/vgsvg"
	"gonum.org/v1/plot/vg/draw"
)

func GenerateLengthHistogramSVG(lengths []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "Read Length Distribution"
	p.X.Label.Text = "Length"
	p.Y.Label.Text = "Count"

	h, err := plotter.NewHist(plotter.Values(lengths), 20) // 20 bins
	if err != nil {
		return "", err
	}
	p.Add(h)

	var buf bytes.Buffer
	canvas := vgsvg.New(vg.Points(500), vg.Points(300))
	dc := draw.New(canvas)
	p.Draw(dc)

	_, err = canvas.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}


func GenerateGCHistogramSVG(gcValues []float64) (string, error) {
	p := plot.New()
	p.Title.Text = "GC Content Distribution"
	p.X.Label.Text = "GC Content (%)"
	p.Y.Label.Text = "Count"

	h, err := plotter.NewHist(plotter.Values(gcValues), 20)
	if err != nil {
		return "", err
	}
	p.Add(h)

	var buf bytes.Buffer
	canvas := vgsvg.New(vg.Points(500), vg.Points(300))
	dc := draw.New(canvas)
	p.Draw(dc)
	_, err = canvas.WriteTo(&buf)
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
	canvas := vgsvg.New(vg.Points(500), vg.Points(300))
	p.Draw(draw.New(canvas))
	_, err = canvas.WriteTo(&buf)
	if err != nil {
		return "", err
	}
	return buf.String(), nil
}
