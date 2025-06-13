package ran_dna_gen

import (
	"flag"
	"fmt"
	"math/rand"
	"os"
	"time"
	"compress/gzip"
)

// Generate a random DNA sequence of given length and GC bias (0.0–1.0)
func randSeq(seqLength int, gcBias float64) string {
	if gcBias < 0.0 || gcBias > 1.0 {
		fmt.Println("GC bias must be between 0.0 and 1.0")
		return ""
	}

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

    fs.Parse(args)

	if *gc < 0.0 || *gc > 0.99 {
		fmt.Println("GC bias must be between 0.0 and 0.99")
		os.Exit(1)
	}

	// Set RNG seed
	if *seed == 0 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(*seed)
	}

	// Generate and wrap the sequence
	sequence := randSeq(*length, *gc)
	fasta := fmt.Sprintf(">%s\n%s", *name, wrapFasta(sequence, 60))

	// Output the result
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
