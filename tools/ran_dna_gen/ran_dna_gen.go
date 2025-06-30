package ran_dna_gen

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

// SequenceRequest holds parameters for one sequence
type SequenceRequest struct {
	ID     string
	Length int
	GCBias float64
}

// MultiSeqFlag parses multiple -seq inputs
type MultiSeqFlag []SequenceRequest

func (m *MultiSeqFlag) String() string {
	return fmt.Sprint(*m)
}

func (m *MultiSeqFlag) Set(value string) error {
	parts := strings.Split(value, ",")
	if len(parts) != 3 {
		return fmt.Errorf("expected format: name,length,gc_bias")
	}
	length, err1 := strconv.Atoi(parts[1])
	gc, err2 := strconv.ParseFloat(parts[2], 64)
	if err1 != nil || err2 != nil || gc < 0.0 || gc > 1.0 {
		return fmt.Errorf("invalid sequence format or values")
	}
	*m = append(*m, SequenceRequest{
		ID:     parts[0],
		Length: length,
		GCBias: gc,
	})
	return nil
}

// Generate a random DNA sequence of given length and GC bias (0.0–1.0)
func randSeq(seqLength int, gcBias float64) string {
	cWeight := gcBias / 2
	aWeight := (1 - gcBias) / 2
	tWeight := (1 - gcBias) / 2

	seq := make([]rune, seqLength)
	for i := 0; i < seqLength; i++ {
		r := rand.Float64()
		switch {
		case r < aWeight:
			seq[i] = 'A'
		case r < aWeight+tWeight:
			seq[i] = 'T'
		case r < aWeight+tWeight+cWeight:
			seq[i] = 'C'
		default:
			seq[i] = 'G'
		}
	}
	return string(seq)
}

// Wrap sequence every `width` characters for FASTA formatting
func wrapFasta(seq string, width int) string {
	var wrapped string
	for i := 0; i < len(seq); i += width {
		end := i + width
		if end > len(seq) {
			end = len(seq)
		}
		wrapped += seq[i:end] + "\n"
	}
	return wrapped
}

func Run(args []string) {
	fs := flag.NewFlagSet("ran_dna_gen", flag.ExitOnError)

	length := fs.Int("length", 100, "Length of generated DNA sequence")
	gc := fs.Float64("gc_bias", 0.5, "GC bias (0.0-1.0)")
	seed := fs.Int64("seed", 0, "Seed for RNG")
	outFile := fs.String("out_file", "", "Output FASTA file")
	name := fs.String("name", "random_seq", "Sequence name (FASTA header)")
	gzip_option := fs.Bool("gzip", false, "Compress output using gzip (.gz)")

	var multiSeq MultiSeqFlag
	fs.Var(&multiSeq, "seq", "Define sequence as 'name,length,gc_bias' (can be repeated)")

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
	if *gc < 0.0 || *gc > 0.99 {
		fmt.Println("GC bias must be between 0.0 and 0.99")
		os.Exit(1)
	}

	if *seed == 0 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(*seed)
	}

	var fastaOut strings.Builder

	if len(multiSeq) > 0 {
		for _, req := range multiSeq {
			seq := randSeq(req.Length, req.GCBias)
			fastaOut.WriteString(fmt.Sprintf(">%s\n%s", req.ID, wrapFasta(seq, 60)))
		}
	} else {
		seq := randSeq(*length, *gc)
		fastaOut.WriteString(fmt.Sprintf(">%s\n%s", *name, wrapFasta(seq, 60)))
	}

	fasta := fastaOut.String()

	if *outFile == "" {
		if *gzip_option {
			fmt.Fprintln(os.Stderr, "Cannot gzip to stdout directly. Please specify an output file.")
			os.Exit(1)
		}
		fmt.Print(fasta)
	} else {
		outputPath := *outFile
		if *gzip_option {
			outputPath += ".gz"
			file, err := os.Create(outputPath)
			if err != nil {
				fmt.Println("Error creating gzip file:", err)
				os.Exit(1)
			}
			defer file.Close()

			writer := gzip.NewWriter(file)
			defer writer.Close()

			_, err = writer.Write([]byte(fasta))
			if err != nil {
				fmt.Println("Error writing compressed data:", err)
				os.Exit(1)
			}
			fmt.Printf("Wrote compressed sequence to %s\n", outputPath)
		} else {
			err := os.WriteFile(outputPath, []byte(fasta), 0644)
			if err != nil {
				fmt.Println("Error writing to file:", err)
				os.Exit(1)
			}
			fmt.Printf("Wrote sequence to %s\n", outputPath)
		}
	}
}
