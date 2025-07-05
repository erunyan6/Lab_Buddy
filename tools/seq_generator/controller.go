package seq_generator

import (
	"compress/gzip"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
	"bufio"
	"io"
)

// For repeated -seq arguments
type SequenceRequest struct {
	ID     string
	Length int
	GCBias float64
}

type MultiSeqFlag []SequenceRequest

func (m *MultiSeqFlag) String() string { return fmt.Sprint(*m) }
func (m *MultiSeqFlag) Set(value string) error {
	parts := strings.Split(value, ",")
	if len(parts) < 2 {
		return fmt.Errorf("expected format: name,length[,gc_bias]")
	}
	length, err := strconv.Atoi(parts[1])
	if err != nil {
		return fmt.Errorf("invalid length")
	}
	gc := 0.5
	if len(parts) == 3 {
		gc, err = strconv.ParseFloat(parts[2], 64)
		if err != nil || gc < 0.0 || gc > 1.0 {
			return fmt.Errorf("invalid gc_bias")
		}
	}
	*m = append(*m, SequenceRequest{ID: parts[0], Length: length, GCBias: gc})
	return nil
}

func WrapFastaToWriter(w io.Writer, seq string, width int) {
	for i := 0; i < len(seq); i += width {
		end := i + width
		if end > len(seq) {
			end = len(seq)
		}
		w.Write([]byte(seq[i:end]))
		w.Write([]byte("\n"))
	}
}

func Run(args []string) {
	fs := flag.NewFlagSet("seq_generator", flag.ExitOnError)

	mode := fs.String("mode", "dna", "Sequence type: dna, rna, or protein")
	name := fs.String("name", "random_seq", "Sequence name")
	length := fs.Int("length", 100, "Sequence length")
	gc := fs.Float64("gc_bias", 0.5, "GC bias for DNA/RNA")
	seed := fs.Int64("seed", 0, "Random seed")
	outFile := fs.String("out_file", "", "Output FASTA file (omit to write to stdout)")
	gzipPreset := fs.String("gzip_preset", "none", "Compression preset: fast, balanced, archival, none")

	var multiSeq MultiSeqFlag
	fs.Var(&multiSeq, "seq", "Use format name,length[,gc_bias] (repeatable)")

	err := fs.Parse(args)
	if err != nil {
		fmt.Println("Error parsing flags:", err)
		os.Exit(1)
	}
	if len(fs.Args()) > 0 {
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}
	if *mode == "protein" && *gc != 0.5 {
		fmt.Fprintln(os.Stderr, "Error: -gc_bias is not applicable in protein mode.")
		os.Exit(1)
	}

	if *mode == "protein" && *gc != 0.5 {
		fmt.Fprintln(os.Stderr, "Warning: -gc_bias has no effect in protein mode and will be ignored.")
	}	

	// Handle compression preset
	var useGzip bool
	var gzipLevel int
	switch strings.ToLower(*gzipPreset) {
	case "fast":
		useGzip = true
		gzipLevel = gzip.BestSpeed // 1
	case "balanced":
		useGzip = true
		gzipLevel = 5
	case "archival":
		useGzip = true
		gzipLevel = 8
	case "none":
		useGzip = false
	default:
		fmt.Fprintf(os.Stderr, "Unknown gzip_preset: %q. Valid options are: fast, balanced, archival, none\n", *gzipPreset)
		os.Exit(1)
	}

	// Handle random seed
	if *seed == 0 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(*seed)
	}

	// Define sequence generation function
	makeSeq := func(length int, gc float64) string {
		switch *mode {
		case "dna":
			return GenerateDNA(length, gc, false)
		case "rna":
			return GenerateDNA(length, gc, true)
		case "protein":
			return GenerateProtein(length)
		default:
			fmt.Fprintf(os.Stderr, "Unknown mode: %s\n", *mode)
			os.Exit(1)
			return ""
		}
	}

	// ===========================
	// OUTPUT TO STDOUT (NO FILE)
	// ===========================
	if *outFile == "" {
		if useGzip {
			fmt.Fprintln(os.Stderr, "Error: cannot write gzipped data to stdout. Please specify -out_file.")
			os.Exit(1)
		}

		writer := bufio.NewWriter(os.Stdout)
		defer writer.Flush()

		if len(multiSeq) > 0 {
			for _, req := range multiSeq {
				fmt.Fprintf(writer, ">%s\n", req.ID)
				WrapFastaToWriter(writer, makeSeq(req.Length, req.GCBias), 60)
			}
		} else {
			fmt.Fprintf(writer, ">%s\n", *name)
			WrapFastaToWriter(writer, makeSeq(*length, *gc), 60)
		}

		return
	}

	// ===========================
	// OUTPUT TO FILE
	// ===========================
	path := *outFile
	if useGzip {
		path += ".gz"
	}

	file, err := os.Create(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating file: %v\n", err)
		os.Exit(1)
	}
	defer file.Close()

	var baseWriter io.Writer
	if useGzip {
		gz, err := gzip.NewWriterLevel(file, gzipLevel)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating gzip writer: %v\n", err)
			os.Exit(1)
		}
		defer gz.Close()
		baseWriter = gz
	} else {
		baseWriter = file
	}

	writer := bufio.NewWriter(baseWriter)
	defer writer.Flush()

	if len(multiSeq) > 0 {
		for _, req := range multiSeq {
			fmt.Fprintf(writer, ">%s\n", req.ID)
			WrapFastaToWriter(writer, makeSeq(req.Length, req.GCBias), 60)
		}
	} else {
		fmt.Fprintf(writer, ">%s\n", *name)
		WrapFastaToWriter(writer, makeSeq(*length, *gc), 60)
	}

	// Final message
	if useGzip {
		fmt.Printf("Wrote compressed sequence to %s using preset %q\n", path, *gzipPreset)
	} else {
		fmt.Printf("Wrote uncompressed sequence to %s\n", path)
	}
}
