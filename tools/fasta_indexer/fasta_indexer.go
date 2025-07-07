package fasta_indexer

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"
)

type FastaIndex struct {
	SeqID        string
	SeqLen       int
	Offset       int64
	BasesPerLine int
	BytesPerLine int
}

func indexFasta(file string) ([]FastaIndex, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer f.Close()

	var reader io.Reader = f
	buf := make([]byte, 2)
	if _, err := f.Read(buf); err == nil && buf[0] == 0x1F && buf[1] == 0x8B {
		f.Seek(0, io.SeekStart)
		gr, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("failed to open gzip reader: %w", err)
		}
		defer gr.Close()
		reader = gr
	} else {
		f.Seek(0, io.SeekStart)
	}

	scanner := bufio.NewScanner(reader)

	var indexes []FastaIndex
	var current FastaIndex
	var byteCount int64 = 0
	var firstSeqLine = true
	var inSequence = false

	for scanner.Scan() {
		line := scanner.Text()
		lineLen := len(line)
		byteCount += int64(lineLen) + 1 	// Add 1 for '\n'

		if strings.HasPrefix(line, ">") {
			if inSequence {					// Save the previous entry
				indexes = append(indexes, current)
			}

			// Start a new record
			current = FastaIndex{
				SeqID: strings.TrimPrefix(line, ">"),
				SeqLen: 0,
				Offset: byteCount,
				BasesPerLine: 0,
				BytesPerLine: 0,
			}
			firstSeqLine = true
			inSequence = true
			continue
		}

		current.SeqLen += len(strings.TrimSpace(line))

		if firstSeqLine{
			current.BasesPerLine = len(strings.TrimSpace(line))
			current.BytesPerLine = lineLen + 1
			firstSeqLine = false
		}
	}

	if inSequence {
		indexes = append(indexes, current)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %w", err)
	}

	return indexes, nil

}

func FastaIndex_Run(args []string) {
	fs := flag.NewFlagSet("fasta_indexer", flag.ExitOnError)
	inFile := fs.String("in_file", "", "Input FASTA file")
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

	if *inFile == "" {
		fmt.Println("Error: -in_file is required to run the fasta_index tool")
		os.Exit(1)
	}

	indexes, err := indexFasta(*inFile)
	if err != nil {
		fmt.Println("Error:", err)
		os.Exit(1)
	}

	path := *inFile + ".fai"
	file, err := os.Create(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating file: %v\n", err)
		os.Exit(1)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	defer writer.Flush()

	// Write each index line
	for _, idx := range indexes {
		fmt.Fprintf(writer, "%s\t%d\t%d\t%d\t%d\n", idx.SeqID, idx.SeqLen, idx.Offset, idx.BasesPerLine, idx.BytesPerLine)
	}

	fmt.Printf("FASTA file %s successfully indexed (%s)\n", *inFile, path)
}

