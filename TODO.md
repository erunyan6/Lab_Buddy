# TODO

## ONGOING

- [ ] **ORF Finder: Fix strand and frame logic**
      Program is outputting all frames of a specific strand, even if asked for a specific frame. Also would like to simplify flag set options to reduce needless complexity.

- [ ] **Benchmark: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **Sanity Check: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **Random DNA Generator: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **ORF Finder: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **FASTA Overview: Finish code annotion and decide on 3bit system**
      Add comments to fully explain the benchmark system.

- [ ] **FASTA 3-Bit Encoder: Finish code annotation**
      Add comments to fully explain the benchmark system.

- [ ] **Organize tools into single directory**
      Move tools to single directory. Refactor import statements to adjust to the new organization.

- [ ] **ORF Finder: Perform program accuracy validation**
      Determine accuracy and credibility of ORF_Finder tool by comparing against established systems.

- [ ] **Random DNA Generator: Add multi-sequence functionality**
      Add ability to accept list of sequence names, desired size, and GC percentages (if appropriate) to simulate multiple sequences in a single run

---

## DONE

- [X] **Kmer_Analyzer: Finish code annotation** (6/27/2025)
      *(Starting at line 97)*: Add comments to fully explain the result formatting and output section.

- [X] **Kmer_Analyzer: Streaming mode for large FASTA files**     (6/27/2025)
      *(Lines 48–53)*: Refactor `countKmers` to use a sliding window instead of building the full DNA sequence in memory. Prevents RAM overload on large files.

- [X] **Kmer_Analyzer: Implement frame-specific k-mer analysis**  (6/27/2025)
      Support frame-specific options (e.g., +1, -2 frames).
