package fastqc_mimic

import (
	"fmt"
	"os"
)

func WriteHTMLReport(filename string, stats FastqStats, svgLength string, svgGC string, svgPQual string,) error {
	f, err := os.Create(filename + ".html")
	if err != nil {
		return err
	}
	defer f.Close()

	html := fmt.Sprintf(`
<!DOCTYPE html>
<html>
<head>
	<title>FASTQC Mimic Report</title>
	<meta charset="UTF-8">
	<style>
		body { font-family: Arial, sans-serif; padding: 20px; background-color: #f9f9f9; }
		h1 { color: #333; }
		.metric { margin: 10px 0; }
		table { border-collapse: collapse; margin-top: 20px; }
		th, td { padding: 8px 12px; border: 1px solid #ccc; text-align: left; }
		th { background-color: #eee; }
	</style>
</head>
<body>
	<h1>FASTQC Mimic Report</h1>
	<table>
		<tr><th>Metric</th><th>Value</th></tr>
		<tr><td>Total Reads</td><td>%d</td></tr>
		<tr><td>Average Read Length</td><td>%.2f</td></tr>
		<tr><td>Min Read Length</td><td>%d</td></tr>
		<tr><td>Max Read Length</td><td>%d</td></tr>
		<tr><td>Read Length StdDev</td><td>%.2f</td></tr>
		<tr><td>GC Content</td><td>%.2f%%</td></tr>
		<tr><td>GC Content StdDev</td><td>%.2f</td></tr>
		<tr><td>N Content</td><td>%.2f%%</td></tr>
		<tr><td>%% Reads with N</td><td>%.2f%%</td></tr>
		<tr><td>%% Reads with Mean Q &lt; 20</td><td>%.2f%%</td></tr>
		<tr><td>%% Bases with Q ≥ 20</td><td>%.2f%%</td></tr>
		<tr><td>%% Bases with Q ≥ 30</td><td>%.2f%%</td></tr>
		<tr><td>Mean Quality Score</td><td>%.2f</td></tr>
		<tr><td>Quality Score StdDev</td><td>%.2f</td></tr>
		<tr><td>Max Homopolymer Run</td><td>%d</td></tr>
		<tr><td>Mean Homopolymer Run</td><td>%.2f</td></tr>
		<tr><td>%% Duplicate Reads (approx.)</td><td>%.2f%%</td></tr>
		<tr><td>Mean Shannon Entropy</td><td>%.2f</td></tr>
		<tr><td>Average A Content</td><td>%.2f%%</td></tr>
		<tr><td>Average T Content</td><td>%.2f%%</td></tr>
		<tr><td>Average C Content</td><td>%.2f%%</td></tr>
		<tr><td>Average G Content</td><td>%.2f%%</td></tr>
	</table>
	<h2>Read Length Distribution</h2>
	<div>%s</div>
	<h2>GC Content Distribution</h2>
	<div>%s</div>
	<h2>Per-Base Read Quality</h2>
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
		svgGC,
		svgPQual,
	)

	_, err = f.WriteString(html)
	return err
}

