# Lab Buddy Main Executable Change Log

| Release Date | Version | Key Updates |
|--------------|---------|-------------|
| July 2025    | v1.9.1  | Removed slice index out of bounds error (bug). Removed bug causing HTML output to be generated when not requested. |
| July 2025    | v1.9.0  | Added FASTA Isolate tool for rapid extraction of specific entries / ranges from FASTA files.  Removed unused 3Bit encoder for the time being. |
| July 2025    | v1.8.2  | Updated FASTQC_Mimic with concurrency in key functions, reducing run time and memory usage by ~ 37%. |
| July 2025    | v1.8.1  | Reduced overhead and addressed bugs associated with FASTQC_Mimic. |
| July 2025    | v1.8.0  | Added FASTQC_Mimic, a Lab_Buddy version of the popular FASTQ file analyzer and graph / report generator. |
| July 2025    | v1.7.0  | Added a sequencing simulator (seq_sim) tool for producing realistic sequencing reads from a DNA FASTA file. |
| July 2025    | v1.6.0  | Added FASTA indexer for enabling efficient and rapid access to FASTA sequences. |
| July 2025    | v1.5.0  | Added chunking helper function to `utils/` directory, increasing upper limit of FASTA input size while further reducing computational usage.|
| July 2025    | v1.4.0  | Added Common package for often-used functions under a new `utils/` directory. |
| July 2025    | v1.3.0  | Added a Lab Buddy Art tool for fun. |
| July 2025    | v1.2.0  | Added main and tool-specific changelogs to support transparent and accurate version tracking. |
| June 2025    | v1.1.0  | Reorganized all tools under a new `tools/` directory for a cleaner home directory. Introduced `config/` directory to store version control information. |
| June 2025    | v1.0.0  | Initial release of the Lab Buddy tool suite — a modular Go application containing lightweight, commonly used bioinformatics utilities under a unified executable. |
