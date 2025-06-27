# TODO

## kmer_analyzer

- [ ] **Streaming mode for large FASTA files**  
      *(Lines 48–53)*: Refactor `countKmers` to use a sliding window instead of building the full DNA sequence in memory. Prevents RAM overload on large files.

- [ ] **Finish code annotation**  
      *(Starting at line 97)*: Add comments to fully explain the result formatting and output section.

- [ ] **Implement strand-specific k-mer analysis**  
      Support strand-specific options (e.g., +1, -2 frames) to enable more granular biological insights, especially for coding regions.
