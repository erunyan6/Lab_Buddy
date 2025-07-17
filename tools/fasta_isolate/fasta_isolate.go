package fasta_isolate

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strings"
	"strconv"
	"io"
	"compress/gzip"
	"path/filepath"

	"lab_buddy_go/tools/fasta_indexer" 
)

type multiString []string

func (s *multiString) String() string {
	return strings.Join(*s, ",")
}

func (s *multiString) Set(value string) error {
	*s = append(*s, value)
	return nil
}

type TargetSpec struct {
	Header     string
	Start, End *int // nil if full range
}

func parseTargetSpec(s string) (TargetSpec, error) {
	var ts TargetSpec
	if strings.Contains(s, ":") && strings.Contains(s, "-") {
		parts := strings.SplitN(s, ":", 2)
		rangeParts := strings.SplitN(parts[1], "-", 2)

		start, err1 := strconv.Atoi(rangeParts[0])
		end, err2 := strconv.Atoi(rangeParts[1])
		if err1 != nil || err2 != nil || start < 0 || end <= start {
			return ts, fmt.Errorf("invalid coordinate range in -seq %s", s)
		}
		ts = TargetSpec{
			Header: parts[0],
			Start:  &start,
			End:    &end,
		}
	} else {
		ts = TargetSpec{Header: s}
	}
	return ts, nil
}


func FastaIsolate_Run(args []string) {
	fs := flag.NewFlagSet("fasta_isolate", flag.ExitOnError)

	inFile := fs.String("in_file", "", "Input FASTA file")
	outFile := fs.String("out_file", "isolated.fasta", "Output FASTA file")
	useIndex := fs.Bool("use_index", false, "Use FASTA index (.fai) for faster extraction")
	var targets multiString
	fs.Var(&targets, "seq", "Header(s) to extract (can repeat -seq multiple times)")

	err := fs.Parse(args)
	if err != nil {
		fmt.Println("Error parsing flags:", err)
		os.Exit(1)
	}

	if *inFile == "" || len(targets) == 0 {
		fmt.Println("Usage: -in_file <file> -out_file <file> -seq <header1> [-seq <header2> ...] [-use_index]")
		os.Exit(1)
	}

	if *useIndex && filepath.Ext(*inFile) == ".gz" {
		fmt.Fprintln(os.Stderr, "Warning: Indexed mode not supported for gzipped files. Using buffered mode instead.")
		*useIndex = false
	}	

	targetSpecs := make(map[string]TargetSpec)
	for _, t := range targets {
		spec, err := parseTargetSpec(t)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Warning: %v (skipping)\n", err)
			continue
		}
		targetSpecs[spec.Header] = spec
	}
	

	if *useIndex {
		// Create index if not already present
		fasta_indexer.FastaIndex_Run([]string{"-in_file", *inFile})
		indexPath := *inFile + ".fai"
		err = extractWithIndex(*inFile, indexPath, *outFile, targetSpecs)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error during index-based extraction: %v\n", err)
			os.Exit(1)
		}
	} else {
		err = extractBuffered(*inFile, *outFile, targetSpecs)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error during buffered extraction: %v\n", err)
			os.Exit(1)
		}
	}
}

func extractBuffered(inPath, outPath string, targets map[string]TargetSpec) error {
	found := make(map[string]bool)
	in, scanner, err := openPossiblyGzipped(inPath)
	if err != nil {
		return err
	}
	defer in.Close()
	
	out, writer, err := createPossiblyGzipped(outPath)
	if err != nil {
		return err
	}
	defer out.Close()	

	var keep bool
	var currentHeader string
	var currentSpec TargetSpec
	var seqBuilder strings.Builder

	flushSequence := func() {
		seq := seqBuilder.String()
		if currentSpec.Start != nil && currentSpec.End != nil {
			start := *currentSpec.Start
			end := *currentSpec.End
			if start >= len(seq) {
				fmt.Fprintf(os.Stderr, "Warning: Start %d beyond length of '%s'\n", start, currentHeader)
				return
			}
			if end > len(seq) {
				end = len(seq)
			}
			seq = seq[start:end]
		}
		for i := 0; i < len(seq); i += 60 {
			end := i + 60
			if end > len(seq) {
				end = len(seq)
			}
			writer.WriteString(seq[i:end] + "\n")
		}
		seqBuilder.Reset()
	}

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			// flush previous record if needed
			if keep {
				flushSequence()
			}
			header := strings.Fields(line[1:])[0]
			spec, ok := targets[header]
			if ok {
				keep = true
				currentHeader = header
				currentSpec = spec
				found[header] = true
				writer.WriteString(">" + header + "\n")
			} else {
				keep = false
			}
		} else if keep {
			seqBuilder.WriteString(strings.TrimSpace(line))
		}
	}
	if keep {
		flushSequence()
	}
	writer.Flush()

	for k := range targets {
		if !found[k] {
			fmt.Fprintf(os.Stderr, "Warning: Header '%s' not found in FASTA\n", k)
		}
	}

	fmt.Printf("Extracted %d record(s) to %s\n", len(found), outPath)

	return nil
}


type FastaIndex struct {
	SeqID        string
	SeqLen       int
	Offset       int64
	BasesPerLine int
	BytesPerLine int
}

func extractWithIndex(fastaPath, indexPath, outPath string, targets map[string]TargetSpec) error {
	// Read index into a map
	indexFile, err := os.Open(indexPath)
	if err != nil {
		return fmt.Errorf("failed to open index file: %w", err)
	}
	defer indexFile.Close()

	indexMap := make(map[string]FastaIndex)
	scanner := bufio.NewScanner(indexFile)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) != 5 {
			continue
		}
		seqLen, _ := strconv.Atoi(fields[1])
		offset, _ := strconv.ParseInt(fields[2], 10, 64)
		basesPerLine, _ := strconv.Atoi(fields[3])
		bytesPerLine, _ := strconv.Atoi(fields[4])

		indexMap[fields[0]] = FastaIndex{
			SeqID:        fields[0],
			SeqLen:       seqLen,
			Offset:       offset,
			BasesPerLine: basesPerLine,
			BytesPerLine: bytesPerLine,
		}
	}
	if err := scanner.Err(); err != nil {
		return err
	}

	// Open FASTA file and prepare output
	fastaFile, err := os.Open(fastaPath)
	if err != nil {
		return fmt.Errorf("failed to open FASTA: %w", err)
	}
	defer fastaFile.Close()

	outFile, err := os.Create(outPath)
	if err != nil {
		return fmt.Errorf("failed to create output file: %w", err)
	}
	defer outFile.Close()
	writer := bufio.NewWriter(outFile)

	found := make(map[string]bool)

	for seqID, spec := range targets {
		idx, ok := indexMap[seqID]
		if !ok {
			fmt.Fprintf(os.Stderr, "Warning: Header '%s' not found in index\n", seqID)
			continue
		}
		found[seqID] = true
		writer.WriteString(">" + seqID + "\n")
	
		start := 0
		end := idx.SeqLen
		if spec.Start != nil && spec.End != nil {
			start = *spec.Start
			end = *spec.End
			if start >= idx.SeqLen {
				fmt.Fprintf(os.Stderr, "Warning: Start %d beyond length of '%s'\n", start, seqID)
				continue
			}
			if end > idx.SeqLen {
				end = idx.SeqLen
			}
		}
	
		startLine := start / idx.BasesPerLine
		endLine := (end - 1) / idx.BasesPerLine
		linesToRead := endLine - startLine + 1
	
		readOffset := idx.Offset + int64(startLine*idx.BytesPerLine)
		_, err := fastaFile.Seek(readOffset, io.SeekStart)
		if err != nil {
			return fmt.Errorf("failed to seek: %w", err)
		}
	
		// Read all required lines
		var seqBuilder strings.Builder
		buf := make([]byte, idx.BytesPerLine)
		for i := 0; i < linesToRead; i++ {
			n, err := fastaFile.Read(buf)
			if err != nil && err != io.EOF {
				return fmt.Errorf("failed to read sequence data: %w", err)
			}
			seqBuilder.WriteString(strings.TrimSpace(string(buf[:n])))
		}
	
		fullSeq := seqBuilder.String()
		subSeq := fullSeq[start:end]
		if len(subSeq) > (end - start) {
			subSeq = subSeq[:end-start]
		}
	
		for i := 0; i < len(subSeq); i += 60 {
			e := i + 60
			if e > len(subSeq) {
				e = len(subSeq)
			}
			writer.WriteString(subSeq[i:e] + "\n")
		}
	}	

	writer.Flush()

	// Extra warning pass (in case index exists but target not found)
	for seqID := range targets {
		if !found[seqID] {
			fmt.Fprintf(os.Stderr, "Warning: Header '%s' not found in index\n", seqID)
		}
	}

	fmt.Printf("Extracted %d record(s) to %s\n", len(found), outPath)

	return nil
}


// Detect gzip input and return buffered reader
func openPossiblyGzipped(path string) (io.ReadCloser, *bufio.Scanner, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}

	var reader io.ReadCloser = file
	if filepath.Ext(path) == ".gz" {
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, nil, err
		}
		reader = gz
	}

	scanner := bufio.NewScanner(reader)
	return reader, scanner, nil
}

// Detect gzip output and return buffered writer
func createPossiblyGzipped(path string) (io.WriteCloser, *bufio.Writer, error) {
	file, err := os.Create(path)
	if err != nil {
		return nil, nil, err
	}

	var writer io.WriteCloser = file
	if filepath.Ext(path) == ".gz" {
		gz := gzip.NewWriter(file)
		writer = gz
	}

	bufWriter := bufio.NewWriter(writer)
	return writer, bufWriter, nil
}
