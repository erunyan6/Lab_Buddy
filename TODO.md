# TODO

## ONGOING

- [ ] **Random DNA Generator: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **ORF Finder: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **FASTA Overview: Finish code annotion and decide on 3bit system**
      Add comments to fully explain the benchmark system.

- [ ] **FASTA 3-Bit Encoder: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **ORF Finder: Explore parallelization opportunities**
      Impliment parallelization of chromosomes or contigs or something to speed up analysis of FASTA files with multiple large sequences.

- [ ] **Sequencing Simulator: Add concurrency**
      Add concurrency to process multiple chromosomes more efficiently

- [ ] **Sequencing Simulator: Add custom help screen**
      There are a lot of options in this tool and organization is DESPARATELY needed

- [ ] **Sequencing Simulator: Increase realism**
      Short reads have increased realism; long reads are untested; goal is to fool antagonistic model (goal: 50% accuracy +/- 3%)

- [ ] **FASTQC_Mimic: Annotate and clean up module**
      Module is messy due to shotgun coding. When able, remove unused functions and annotate code for readability

---

## DONE

- [X] **Sequencing Simulator: Introduce realistic simulator presets**   (7/9/2025)
      Mimic common sequencing techniques with accurate length, error rates, quality scores, adapters, etc

- [X] **Next Step to Annotation: Figure out how to work ORF_Finder output for BLAST**     (7/6/2025)
      Design protein translation tool to configure ORF_Finder output for BLAST and downstream annotation tools.

- [X] **ORF Finder: Perform program accuracy validation**   (7/5/2025)
      Determine accuracy and credibility of ORF_Finder tool by comparing against established systems.

- [X] **ORF Finder: Fix strand and frame logic**      (7/4/2025)
      Program is outputting all frames of a specific strand, even if asked for a specific frame. Also would like to simplify flag set options to reduce needless complexity.

- [X] **Utils: Make common/utils package(s)**   (7/1/2025)
      Consolidate commonly used functions to avoid duplicate code.

- [X] **All tools: Add changelogs** (7/1/2025)
      Add change logs for all tools to display update changes.

- [X] **Random Protein Generator: Create new tool**   (6/30/2025)
      Expand example generation capacity with protein generation options.

- [X] **FASTA Overview: Add protein FASTA capacity**  (6/30/2025)
      Expand FASTA summary capability to protein FASTA files as well as DNA.

- [X] **FASTA Overview: Known error with duplicate header reporting**   (6/30/2025)
      With duplicate headers, information from the later sequence (gc bais, length, etc) overrides the first without warning.

- [X] **Random DNA Generator: Add multi-sequence functionality**  (6/29/2025)
      Add ability to accept list of sequence names, desired size, and GC percentages (if appropriate) to simulate multiple sequences in a single run.

- [X] **Benchmark: Finish code annotation**     (6/29/2025)
      Add comments to fully explain the benchmark system.

- [X] **Organize tools into single directory**  (6/29/2025)
      Move tools to single directory. Refactor import statements to adjust to the new organization.

- [X] **Sanity Check: Finish code annotation**  (6/29/2025)
      Add comments to fully explain the benchmark system.

- [X] **Kmer_Analyzer: Finish code annotation** (6/27/2025)
      *(Starting at line 97)*: Add comments to fully explain the result formatting and output section.

- [X] **Kmer_Analyzer: Streaming mode for large FASTA files**     (6/27/2025)
      *(Lines 48–53)*: Refactor `countKmers` to use a sliding window instead of building the full DNA sequence in memory. Prevents RAM overload on large files.

- [X] **Kmer_Analyzer: Implement frame-specific k-mer analysis**  (6/27/2025)
      Support frame-specific options (e.g., +1, -2 frames).
