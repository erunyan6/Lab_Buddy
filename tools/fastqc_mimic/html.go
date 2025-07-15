package fastqc_mimic

import (
	"fmt"
	"os"
)

// WriteHTMLReport writes the final report.
func WriteHTMLReport(
	filename string,
	stats FastqStats,
	svgLength string,
	svgGC string,
	svgPQual string,
	svgRQuality string,
	svgBaseContent string,
	svgDuplication string,
	svgKmerEnrichment string,
	svgGCBase string,
) error {
	f, err := os.Create(filename + ".html")
	if err != nil {
		return err
	}
	defer f.Close()

	html := fmt.Sprintf(`
<!DOCTYPE html>
<html>
<head>
	<meta charset="UTF-8">
	<title>FASTQC Mimic Report</title>
	<style>
		body { font-family: Arial, sans-serif; padding: 20px; background: #f9f9f9; }
		h1, h2 { color: #333; }
		table { border-collapse: collapse; margin: 20px 0; }
		th, td { border: 1px solid #ccc; padding: 8px 12px; }
		th { background: #eee; }
		svg { background: #fff; border: 1px solid #ccc; margin: 10px 0; }
	</style>
</head>
<body>
	<h1>FASTQC Mimic Report</h1>

	<h2>Summary Statistics</h2>
	<table>
		<tr><th>Metric</th><th>Value</th></tr>
		<tr><td>Total Reads</td><td>%d</td></tr>
		<tr><td>Average Read Length</td><td>%.2f</td></tr>
		<tr><td>Min Read Length</td><td>%d</td></tr>
		<tr><td>Max Read Length</td><td>%d</td></tr>
		<tr><td>Length StdDev</td><td>%.2f</td></tr>
		<tr><td>GC Content</td><td>%.2f%%</td></tr>
		<tr><td>GC Content StdDev</td><td>%.2f</td></tr>
		<tr><td>N Content</td><td>%.2f%%</td></tr>
		<tr><td>Reads with N</td><td>%.2f%%</td></tr>
		<tr><td>Reads with Mean Q&lt;20</td><td>%.2f%%</td></tr>
		<tr><td>Bases with Q≥20</td><td>%.2f%%</td></tr>
		<tr><td>Bases with Q≥30</td><td>%.2f%%</td></tr>
		<tr><td>Mean Quality Score</td><td>%.2f</td></tr>
		<tr><td>Quality Score StdDev</td><td>%.2f</td></tr>
		<tr><td>Max Homopolymer Run</td><td>%d</td></tr>
		<tr><td>Mean Homopolymer Run</td><td>%.2f</td></tr>
		<tr><td>Approx Duplicate Reads</td><td>%.2f%%</td></tr>
		<tr><td>Mean Shannon Entropy</td><td>%.2f</td></tr>
		<tr><td>Average A Content</td><td>%.2f%%</td></tr>
		<tr><td>Average T Content</td><td>%.2f%%</td></tr>
		<tr><td>Average C Content</td><td>%.2f%%</td></tr>
		<tr><td>Average G Content</td><td>%.2f%%</td></tr>
	</table>

	<h2>Read Length Distribution</h2>
	<div>%s</div>

	<h2>Per-Base GC Content</h2>
	<p>This plot shows the GC percentage at each base position across all reads.</p>
	<div>%s</div>

	<h2>Per Sequence GC Content</h2>
	<p>This plot compares observed per-read GC content to a modeled normal distribution.</p>
	<div>%s</div>

	<h2>Per Base Quality Scores</h2>
	<p>Boxplots of base qualities across all reads.</p>
	<div>%s</div>

	<h2>Per Read Mean Quality</h2>
	<p>Distribution of average quality scores per read.</p>
	<div>%s</div>

	<h2>Per Base Sequence Content</h2>
	<p>Proportion of A, C, G, T, and N bases at each position.</p>
	<div>%s</div>

	<h2>Sequence Duplication Levels</h2>
	<p>Proportion of reads with different duplication counts.</p>
	<div>%s</div>

	<h2>K-mer Enrichment</h2>
	<p>Relative enrichment of the most common k-mers across read positions.</p>
	<div>%s</div>
</body>
</html>`,
		stats.TotalReads,
		stats.AvgLength,
		stats.MinLength,
		stats.MaxLength,
		stats.LengthStdDev,
		stats.GCContent,
		stats.GCStdDev,
		stats.NContent,
		stats.ReadsWithNPercent,
		stats.LowQualityReadPercent,
		stats.Q20BasePercent,
		stats.Q30BasePercent,
		stats.MeanQual,
		stats.StdQual,
		stats.MaxHomopolymer,
		stats.MeanHomopolymer,
		stats.ApproxDuplicatePercent,
		stats.MeanEntropy,
		stats.AvgAContent,
		stats.AvgTContent,
		stats.AvgCContent,
		stats.AvgGContent,
		svgLength,
		svgGCBase,
		svgGC,
		svgPQual,
		svgRQuality,
		svgBaseContent,
		svgDuplication,
		svgKmerEnrichment,
	)

	_, err = f.WriteString(html)
	return err
}
