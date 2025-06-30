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

func Run(args []string) {
	fs := flag.NewFlagSet("seq_generator", flag.ExitOnError)

	mode := fs.String("mode", "dna", "Sequence type: dna, rna, or protein")
	name := fs.String("name", "random_seq", "Sequence name")
	length := fs.Int("length", 100, "Sequence length")
	gc := fs.Float64("gc_bias", 0.5, "GC bias for DNA/RNA")
	seed := fs.Int64("seed", 0, "Random seed")
	outFile := fs.String("out_file", "", "Output FASTA file")
	gzipOut := fs.Bool("gzip", false, "Compress output with gzip")

	var multiSeq MultiSeqFlag
	fs.Var(&multiSeq, "seq", "Use format name,length[,gc_bias] (repeatable)")

	err := fs.Parse(args)										// Parse inputs 
	if err != nil {
		fmt.Println("Error parsing flags:", err)				// Check for outright input failures
		os.Exit(1)												// E.g., expected int by recieved str
	}

	if len(fs.Args()) > 0 {										// If unparsed arguments remain:
		fmt.Printf("Unrecognized arguments: %v\n", fs.Args())	// Flag the error and report it
		fmt.Println("Use -h to view valid flags.")
		os.Exit(1)
	}

	if *mode == "protein" && *gc != 0.5 {
		fmt.Fprintln(os.Stderr, "Error: -gc_bias is not applicable in protein mode.")
		os.Exit(1)
	}
	

	if *seed == 0 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(*seed)
	}

	var fastaOut strings.Builder

	makeSeq := func(id string, length int, gc float64) string {
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

	if len(multiSeq) > 0 {
		for _, req := range multiSeq {
			seq := makeSeq(req.ID, req.Length, req.GCBias)
			fastaOut.WriteString(fmt.Sprintf(">%s\n%s", req.ID, WrapFasta(seq, 60)))
		}
	} else {
		seq := makeSeq(*name, *length, *gc)
		fastaOut.WriteString(fmt.Sprintf(">%s\n%s", *name, WrapFasta(seq, 60)))
	}

	output := fastaOut.String()
	if *outFile == "" {
		if *gzipOut {
			fmt.Fprintln(os.Stderr, "Cannot gzip to stdout. Specify -out_file.")
			os.Exit(1)
		}
		fmt.Print(output)
		return
	}

	path := *outFile
	if *gzipOut {
		path += ".gz"
		file, err := os.Create(path)
		if err != nil {
			fmt.Println("Error creating file:", err)
			os.Exit(1)
		}
		defer file.Close()
		gz := gzip.NewWriter(file)
		defer gz.Close()
		_, err = gz.Write([]byte(output))
		if err != nil {
			fmt.Println("Error writing compressed data:", err)
			os.Exit(1)
		}
	} else {
		err := os.WriteFile(path, []byte(output), 0644)
		if err != nil {
			fmt.Println("Error writing file:", err)
			os.Exit(1)
		}
	}
	fmt.Printf("Wrote sequence to %s\n", path)
}
