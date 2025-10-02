# Genomic-Data-Exploration-Pipeline
This Python pipeline performs genomic data analysis by extracting and manipulating sequences from FASTA, FASTQ, and NCBI Entrez. It calculates GC content, nucleotide counts, complements, transcription, translation, identifies longest sequences, computes quality scores, and visualizes sequence lengths and read quality distributions.

## Features
- **FASTA Processing:** Extracts sequence and header, calculates GC content, generates complements, reverse complements, transcribed mRNA, and translated protein sequences.
- **FASTQ Processing:** Parses reads, calculates per-read Phred quality scores, counts nucleotides, identifies positions of ambiguous bases (N), and visualizes read quality distributions.
- **NCBI Entrez Integration:** Searches and fetches gene sequences, supports multiple queries, and retrieves detailed annotations.
- **Visualization:** Generates histograms for sequence length distribution and boxplots for read quality across positions.
- **File Handling:** Reads from FASTA/FASTQ, writes processed sequences to new files.
- **Extensible:** Modular design allows adding custom sequence analyses or plots.

## Requirements
- Python 3.13.5  
- [Biopython](https://biopython.org/)  
- [Matplotlib](https://matplotlib.org/)  
- [Seaborn](https://seaborn.pydata.org/)  
- [Requests](https://docs.python-requests.org/)  
- [Certifi](https://pypi.org/project/certifi/)

**Install dependencies via pip:**
```bash
pip install biopython matplotlib seaborn requests certifi
```
## Usage
**Update your NCBI email in the script:**
  - python
  Copy code
  Entrez.email = "your.email@example.com"
  Prepare input files: FASTA (Long.txt) or FASTQ (SRR390728_1.fastq).

**Run the script:**
  bash
  Copy code
  python genomic_pipeline.py
  Outputs include processed sequences, GC content statistics, quality score analysis, motif occurrences, and visual plots.

## Outputs
  output.txt – Extracted FASTA sequences.
  Plots – Sequence length distribution and read quality boxplots.
  Terminal/console summary – Sequence statistics, nucleotide frequencies, motif analysis, and read quality metrics.

## License
MIT License
