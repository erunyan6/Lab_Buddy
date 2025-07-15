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
	ambigRate := fs.Float64("ambig_rate", 0.0, "Probability of substituting a base with 'N'")
	errorRate := fs.Float64("error_rate", 0.0, "Base substitution error rate (e.g., 0.01 for 1%)")
	indelRate := fs.Float64("indel_rate", 0.0, "Insertion/deletion rate (e.g., 0.001 for 0.1%)")
	readLenMean := fs.Int("read_len_mean", 150, "Mean read length")
	readLenStdDev := fs.Int("read_len_stddev", 0, "Standard deviation for read length (0 = fixed)")
	readLenMin := fs.Int("read_len_min", 50, "Minimum read length")
	readLenMax := fs.Int("read_len_max", 50000, "Maximum read length")
	qualityProfile := fs.String("quality_profile", "short", "Quality score profile: short (Illumina-style) or long (PacBio-style)")
	logErrors := fs.Bool("log", false, "Log sequencing error coordinates and mutations")
	clusterBias := fs.Float64("cluster_bias", 2.0, "Multiplier for error rate after a previous error")
	gcBoost := fs.Float64("sub_rate_gc_boost", 1.5, "Multiplier for substitution rate in high-GC windows")
	maxIndel := fs.Int("max_indel_len", 3, "Maximum indel length (insertions and deletions)")
	homoBoost := fs.Float64("homopolymer_multiplier", 2.0, "Indel rate multiplier in homopolymer regions")
	paired := fs.Bool("paired", false, "Enable paired-end sequencing simulation")
	fragLenMean := fs.Int("frag_len_mean", 600, "Mean DNA fragment length for paired-end sequencing")
	fragLenStddev := fs.Int("frag_len_stddev", 150, "Standard deviation of fragment length")
	splitReads := fs.Bool("split_reads", false, "Output paired-end reads into separate files (R1 and R2)")

	platform := fs.String("platform", "", "Preset platform type (e.g., illumina_hiseq, pacbio_hifi, ont_minion, etc.)")

	var multiSeq MultiSeqFlag
	fs.Var(&multiSeq, "range", "Use format <Header>,[<start>,<end>] (repeatable)")

	// Custom help screen
	fs.Usage = func() {
		fmt.Fprintln(os.Stderr, "Lab Buddy | seq_sim - DNA Sequencing Read Simulator")
		fmt.Fprintln(os.Stderr, "---------------------------------------------")
		fmt.Fprintln(os.Stderr, "Usage: lab_buddy seq_sim [options]")
		fmt.Fprintln(os.Stderr, "\nRequired:")
		fmt.Fprintln(os.Stderr, "  -in_file string           Input FASTA file for sequencing simulation")
	
		fmt.Fprintln(os.Stderr, "\nOptional Output:")
		fmt.Fprintln(os.Stderr, "  -out_file string          Output FASTQ file (default: stdout)")
		fmt.Fprintln(os.Stderr, "  -split_reads              Output paired-end reads into R1 and R2 files")
	
		fmt.Fprintln(os.Stderr, "\nSequencing Parameters:")
		fmt.Fprintln(os.Stderr, "  -read_len int             Fixed read length (default: 150)")
		fmt.Fprintln(os.Stderr, "  -depth int                Target coverage depth (default: 5)")
		fmt.Fprintln(os.Stderr, "  -paired                   Enable paired-end simulation")
		fmt.Fprintln(os.Stderr, "  -frag_len_mean int        Mean fragment length for paired-end (default: 600)")
		fmt.Fprintln(os.Stderr, "  -frag_len_stddev int      Fragment length stddev (default: 150)")
	
		fmt.Fprintln(os.Stderr, "\nLength Distribution:")
		fmt.Fprintln(os.Stderr, "  -read_len_mean int        Mean read length (default: 150)")
		fmt.Fprintln(os.Stderr, "  -read_len_stddev int      Stddev for read length (0 = fixed)")
		fmt.Fprintln(os.Stderr, "  -read_len_min int         Minimum read length (default: 50)")
		fmt.Fprintln(os.Stderr, "  -read_len_max int         Maximum read length (default: 50000)")
	
		fmt.Fprintln(os.Stderr, "\nError Simulation:")
		fmt.Fprintln(os.Stderr, "  -error_rate float         Base substitution rate [0.0–1.0]")
		fmt.Fprintln(os.Stderr, "  -indel_rate float         Indel rate [0.0–1.0]")
		fmt.Fprintln(os.Stderr, "  -ambig_rate float         N-substitution rate [0.0–1.0]")
		fmt.Fprintln(os.Stderr, "  -cluster_bias float       Error multiplier after error (default: 2.0)")
		fmt.Fprintln(os.Stderr, "  -sub_rate_gc_boost float  Substitution rate boost in GC-rich regions (default: 1.5)")
		fmt.Fprintln(os.Stderr, "  -max_indel_len int        Maximum indel length (default: 3)")
		fmt.Fprintln(os.Stderr, "  -homopolymer_multiplier float  Indel boost in homopolymer regions (default: 2.0)")
	
		fmt.Fprintln(os.Stderr, "\nPlatform Presets:")
		fmt.Fprintln(os.Stderr, "  -platform string          Use preset for sequencing platform:")
		fmt.Fprintln(os.Stderr, "                             illumina_hiseq")
		fmt.Fprintln(os.Stderr, "                             illumina_novaseq")
		fmt.Fprintln(os.Stderr, "                             illumina_miseq")
		fmt.Fprintln(os.Stderr, "                             pacbio_hifi")
		fmt.Fprintln(os.Stderr, "                             pacbio_ccs")
		fmt.Fprintln(os.Stderr, "                             ont_minion")
		fmt.Fprintln(os.Stderr, "                             ont_promethion")
			
		fmt.Fprintln(os.Stderr, "\nOther:")
		fmt.Fprintln(os.Stderr, "  -quality_profile string   Quality style: short (Illumina) or long (PacBio)")
		fmt.Fprintln(os.Stderr, "  -log                      Log all simulated error positions")
		fmt.Fprintln(os.Stderr, "  -range <Header>,[start,end]  Limit simulation to a specific region (repeatable)")
	
		fmt.Fprintln(os.Stderr, "\nExample:")
		fmt.Fprintln(os.Stderr, "  lab_buddy seq_sim -in_file genome.fa -depth 10 -platform illumina_miseq")
		fmt.Println()
	}
	
	for _, arg := range args {
		if arg == "-h" || arg == "--help" {
			fs.Usage()
			os.Exit(0)
		}
	}
	
	if len(args) == 0 {
		fs.Usage()
		os.Exit(1)
	}

	err := fs.Parse(args)
	if err != nil {
		fs.Usage()
		os.Exit(1)
	}

	if len(args) == 0 {
		fs.Usage()
		os.Exit(1)
	}

	if *errorRate < 0 || *errorRate > 1 {
		log.Fatal("Error: -error_rate must be between 0.0 and 1.0")
	}
	if *indelRate < 0 || *indelRate > 1 {
		log.Fatal("Error: -indel_rate must be between 0.0 and 1.0")
	}
	if *ambigRate < 0 || *ambigRate > 1 {
		log.Fatal("Error: -ambig_rate must be between 0.0 and 1.0")
	}

	if *platform != "" {
		switch strings.ToLower(*platform) {
		case "illumina_hiseq":
			*readLenMean = 150
			*readLenStdDev = 0
			*readLenMin = 150
			*readLenMax = 150
			*qualityProfile = "short"
			*errorRate = 0.001
			*indelRate = 0.0001
			*ambigRate = 0.0
			*clusterBias = 1.5
			*gcBoost = 1.2
			*maxIndel = 1
			*homoBoost = 1.0
			*paired = true
			*fragLenMean = 400
			*fragLenStddev = 50
			*splitReads = true
	
		case "illumina_novaseq":
			*readLenMean = 250
			*readLenStdDev = 10
			*readLenMin = 200
			*readLenMax = 300
			*qualityProfile = "short"
			*errorRate = 0.002
			*indelRate = 0.0005
			*ambigRate = 0.0005
			*clusterBias = 2.0
			*gcBoost = 1.3
			*maxIndel = 2
			*homoBoost = 1.5
			*paired = true
			*fragLenMean = 600
			*fragLenStddev = 100
			*splitReads = true
	
		case "pacbio_hifi":
			*readLenMean = 15000
			*readLenStdDev = 2000
			*readLenMin = 5000
			*readLenMax = 25000
			*qualityProfile = "long"
			*errorRate = 0.005
			*indelRate = 0.002
			*ambigRate = 0.001
			*clusterBias = 1.2
			*gcBoost = 1.1
			*maxIndel = 3
			*homoBoost = 1.2
			*paired = false
	
		case "ont_minion":
			*readLenMean = 8000
			*readLenStdDev = 2500
			*readLenMin = 1000
			*readLenMax = 20000
			*qualityProfile = "long"
			*errorRate = 0.08
			*indelRate = 0.03
			*ambigRate = 0.005
			*clusterBias = 2.5
			*gcBoost = 1.5
			*maxIndel = 5
			*homoBoost = 3.5
			*paired = false
	
		case "ont_promethion":
			*readLenMean = 12000
			*readLenStdDev = 3000
			*readLenMin = 2000
			*readLenMax = 30000
			*qualityProfile = "long"
			*errorRate = 0.07
			*indelRate = 0.025
			*ambigRate = 0.003
			*clusterBias = 2.2
			*gcBoost = 1.4
			*maxIndel = 6
			*homoBoost = 3.2
			*paired = false
	
		case "illumina_miseq":
			*readLenMean = 250
			*readLenStdDev = 3
			*readLenMin = 243
			*readLenMax = 253
			*qualityProfile = "short"
			*errorRate = 0.002
			*indelRate = 0.0003
			*ambigRate = 0.0002
			*clusterBias = 1.2
			*gcBoost = 1.1
			*maxIndel = 1
			*homoBoost = 1.1
			*paired = true		
			*splitReads = false
	
		case "pacbio_ccs":
			*readLenMean = 15000
			*readLenStdDev = 4000
			*readLenMin = 1000
			*readLenMax = 30000
			*qualityProfile = "long"
			*errorRate = 0.01
			*indelRate = 0.001
			*ambigRate = 0.001
			*clusterBias = 1.5
			*gcBoost = 1.2
			*maxIndel = 2
			*homoBoost = 2.0
			*paired = false
	
		default:
			log.Fatalf("Unknown platform preset: %s", *platform)
		}
	}

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
			log.Printf("Warning: ID %s is not found in FASTA index. Skipping.\n", region.ID)
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
	
		if *paired {
			// PAIR-END MODE
			var w1, w2 io.Writer
	
			if *splitReads {
				// Create R1 and R2 files
				r1Name := strings.TrimSuffix(*outFile, ".fq") + "_R1.fq"
				r2Name := strings.TrimSuffix(*outFile, ".fq") + "_R2.fq"
	
				f1Handle, err := os.Create(r1Name)
				if err != nil {
					log.Fatalf("failed to create R1 output file: %v", err)
				}
				defer f1Handle.Close()
	
				f2Handle, err := os.Create(r2Name)
				if err != nil {
					log.Fatalf("failed to create R2 output file: %v", err)
				}
				defer f2Handle.Close()
	
				w1 = bufio.NewWriter(f1Handle)
				defer w1.(*bufio.Writer).Flush()
				w2 = bufio.NewWriter(f2Handle)
				defer w2.(*bufio.Writer).Flush()
	
			} else {
				// Interleaved mode
				w1 = bufOut
				w2 = bufOut
			}
	
			err := simulateRegionPaired(
				*inFile, index_map, region.ID, start, stop,
				*fragLenMean, *fragLenStddev,
				*readLenMin, *readLenMax,
				*coverageDepth,
				w1, w2,
				*errorRate, *indelRate, *ambigRate,
				*qualityProfile, *logErrors,
				*clusterBias, *gcBoost, *maxIndel, *homoBoost,
			)
	
			if err != nil {
				log.Printf("Paired-end simulation failed for %s [%d-%d]: %v\n", region.ID, start, stop, err)
			}
	
		} else {
			// SINGLE-END MODE
			err := simulateRegion(
				*inFile, index_map, region.ID, start, stop,
				*readLenMean, *readLenStdDev, *readLenMin, *readLenMax,
				*coverageDepth, bufOut,
				*errorRate, *indelRate, *ambigRate,
				*qualityProfile, *logErrors,
				*clusterBias, *gcBoost, *maxIndel, *homoBoost,
			)
	
			if err != nil {
				log.Printf("Simulation failed for %s [%d-%d]: %v\n", region.ID, start, stop, err)
			}
		}
	}
	fmt.Printf("Completed simulation for %d region(s).\n", len(multiSeq))
}
