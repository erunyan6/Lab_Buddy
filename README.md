# Lab Buddy

**Lab Buddy** is a modular command-line toolkit written in Go for laboratory data processing, simulation, and genomic utilities. The system is designed to be expandable, with each tool operating as an independent module that parses its own arguments. This allows for clean, scalable growth of laboratory functions under a single executable.

---

## 📦 Current Tools

| Tool | Description |
| ---- | ----------- |
| `sanity_check` | Quick program sanity check |
| `ran_dna_gen` | Random DNA sequence generator with configurable length, GC bias, and output options |
| `kmer_analyzer` | Simple k-mer counter with support for relative frequencies, sorting, and filtering options |
| `orf_finder` | Open reading frame (ORF) detector supporting custom start codons, strand selection, frame filtering, minimum length, and multiple output formats (TSV/GFF3) |
| `fasta3bit` | Encoder to compress FASTA into custom 3-bit binary format for future tools |
| `fasta_overview` | Quick FASTA report and sanity check |

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
