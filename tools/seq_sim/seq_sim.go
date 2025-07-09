package seq_sim

import (
	"flag"
	"fmt"
	"log"
	"os"
	"io"
	"strconv"
	"strings"
	"compress/gzip"
	"bufio"

	"lab_buddy_go/tools/fasta_indexer"
	"lab_buddy_go/utils"
)

// For repeated -range arguments
type SequenceRequest struct {
	ID    string
	Start int
	Stop  int
}

type MultiSeqFlag []SequenceRequest
func (m *MultiSeqFlag) String() string { return fmt.Sprint(*m) }
func (m *MultiSeqFlag) Set(value string) error {
	parts := strings.Split(value, ",")

	if len(parts) < 1 || parts[0] == "" {
		return fmt.Errorf("missing sequence ID")
	}

	var start, stop int = -1, -1 // default values

	// Parse start if provided and not empty
	if len(parts) > 1 && parts[1] != "" {
		var err error
		start, err = strconv.Atoi(parts[1])
		if err != nil {
			return fmt.Errorf("invalid start: %v", err)
		}
	}

	// Parse stop if provided and not empty
	if len(parts) > 2 && parts[2] != "" {
		var err error
		stop, err = strconv.Atoi(parts[2])
		if err != nil {
			return fmt.Errorf("invalid stop: %v", err)
		}
	}

	if start >= 0 && stop >= 0 && start >= stop {
		return fmt.Errorf("start must be less than stop: got %d >= %d", start, stop)
	}

	*m = append(*m, SequenceRequest{ID: parts[0], Start: start, Stop: stop})
	return nil
}


func SeqSimRun(args []string) {

	// Gather Arguments
	fs := flag.NewFlagSet("seq_sim", flag.ExitOnError)
	inFile := fs.String("in_file", "", "Input FASTA file for sequencing simulation")
	outFile := fs.String("out_file", "", "Output FASTQ file (default: stdout)")
	readLen := fs.Int("read_len", 150, "Length of sequencing reads")
	coverageDepth := fs.Int("depth", 5, "Coverage depth of sequencing")

	var multiSeq MultiSeqFlag
	fs.Var(&multiSeq, "range", "Use format <Header>,[<start>,<end>] (repeatable)")

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
		log.Fatal("Error: -in_file is required")
	}

	if *readLen < 10 {
		log.Fatal("Error: readlen must be a whole integer higher than 10")
	}
	if *coverageDepth < 1 {
		log.Fatal("Error: depth must be a whole integer higher than 1")
	}
	
	// Index FASTA
	fasta_indexer.FastaIndex_Run([]string{"-in_file", *inFile})
	fasta_index := *inFile + ".fai"

	// Check FASTA Index freshness
	common.CheckIndexFreshness(*inFile, fasta_index)

	// Parse FASTA index into map
	index_map, err := parse_fai(fasta_index)
	if err != nil {
		log.Fatalf("failed to parse FASTA index file: %v", err)
	}

	// If no -range provided, simulate entire FASTA
	if len(multiSeq) == 0 {
		fmt.Println("No -range provided, simulating entire FASTA file...")
		for id, rec := range index_map {
			multiSeq = append(multiSeq, SequenceRequest{
				ID:    id,
				Start: 0,
				Stop:  rec.SeqLen,
			})
		}
	}

	var out io.Writer
	var outFileHandle *os.File
	
	if *outFile != "" {
		file, err := os.Create(*outFile)
		if err != nil {
			log.Fatalf("failed to create output file: %v", err)
		}
		outFileHandle = file
		defer outFileHandle.Close()
	
		if strings.HasSuffix(*outFile, ".gz") {
			gz := gzip.NewWriter(file)
			defer gz.Close()
			out = gz
		} else {
			out = file
		}
	} else {
		out = os.Stdout
	}
	
	// Wrap in bufio for speed
	bufOut := bufio.NewWriter(out)
	defer bufOut.Flush()

	// simulate region function here
	for _, region := range multiSeq {
		idx, ok := index_map[region.ID]
		if !ok {
			log.Printf("Warning: ID %s is not found in FASTA index. skipping.\n", region.ID)
			continue
		}
		start := region.Start
		stop := region.Stop

		// Default to full region if start/stop are not set
		if start == -1 {
			start = 0
		}
		if stop == -1 || stop > idx.SeqLen {
			stop = idx.SeqLen
		}
				
		err := simulateRegion(*inFile, index_map, region.ID, start, stop, *readLen, *coverageDepth, bufOut)
		if err != nil {
			log.Printf("Simulation failed for %s [%d-%d]: %v\n", region.ID, start, stop, err)
		}
	}
	fmt.Printf("Completed simulation for %d region(s).\n", len(multiSeq))
}
