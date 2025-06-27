# TODO

## Ongoing





---


## DONE

- [X] **Kmer_Analyzer: Finish code annotation** (6/27/2025)
      *(Starting at line 97)*: Add comments to fully explain the result formatting and output section.

- [X] **Kmer_Analyzer: Streaming mode for large FASTA files**     (6/27/2025)
      *(Lines 48–53)*: Refactor `countKmers` to use a sliding window instead of building the full DNA sequence in memory. Prevents RAM overload on large files.

- [X] **Kmer_Analyzer: Implement frame-specific k-mer analysis**  (6/27/2025)
      Support frame-specific options (e.g., +1, -2 frames).
