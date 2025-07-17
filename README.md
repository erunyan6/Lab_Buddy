# Lab Buddy

**Lab Buddy** is a modular command-line toolkit written in Go for laboratory data processing, simulation, and genomic utilities. The system is designed to be expandable, with each tool operating as an independent module that parses its own arguments. This allows for clean, scalable growth of laboratory functions under a single executable.

---

## 📦 Current Tools

| Tool | Description |
| ---- | ----------- |
| `sanity_check` | Quick program sanity check |
| `seq_generator` | Random DNA/RNA/Protein sequence generator with configurable length, GC bias, and output options |
| `kmer_analyzer` | Efficient streaming k-mer counter with support for strand-specific analysis, reading frames, relative frequency, and sorting/filtering options |
| `orf_finder` | Open reading frame (ORF) detector supporting custom start codons, strand selection, frame filtering, minimum length, and usable output GFF3 |
| `fasta3bit` | Encoder to compress FASTA into custom 3-bit binary format for future tools |
| `fasta_overview` | Quick FASTA report and sanity check. Can be used on DNA, RNA, or Protein 'FASTA' files |
| `benchmark` | Reports enviroment variables and resource usage (RAM, GC cycles, execution time, etc.) of any other tool |
| `fasta_indexer` | Recreation of commonly used '.fai' index file for efficient FASTA access |
| `lab_buddy_art` | ASCII art of Lab Buddy himself, accompanied by a motivational quote or pun |
| `orf_to_faa` | Lightweight protein translator utilizing ORFs identified by the `orf_finder` tool |
| `seq_sim` | Rapid and memory efficient tool mimicking advanced sequencing platforms with realistic error types and probabilities |
| `fastqc_mimic` | FASTQ format analyzer similar in design and output to a mimimized version of the popular package FASTQC |

---

## 🔧 Build

You can build Lab Buddy from source:

```bash
git clone https://github.com/your-username/lab_buddy_go.git
cd lab_buddy_go
go build -o lab_buddy .

```

As this passion project develops, **built** versions will be uploaded.

## 🚀 Usage

```bash
./lab_buddy <tool_name> [flags]
```

## ⚠️ Disclaimer

This project is under active development and is not yet intended for professional use.

## 👨‍💻 Auther

Built and maintained by Elliott Runyan — feel free to fork, contribute, or suggest new tools!
