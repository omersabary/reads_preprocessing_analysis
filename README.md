# DNA Reads Processing and Barcode (Index) Matching Pipeline

This script (written in Python) processes DNA sequencing reads and performs adapter trimming and barcode (index) matching against a known design list (provided in a csv file). 
The output of the script is a CSV file with trimmed (preprocessed) sequences, their counts, and their matched barcode/index (variant IDs).

The script was used in the following paper:
And is based on SOLQC tool (see paper [link](https://academic.oup.com/bioinformatics/article/37/5/720/5896982))

---

## The Pipeline
1. Adapter/primer detection and trimming (forward and reverse orientation, beginning and end of the reads)
2. Barcode (index) extraction and matching using Levenshtein or Hamming distance (user selected metric)
3. Output of a csv file that matches each read with its barcode (index), the script also aggregates reads with exact same sequence and include counts per unique read sequence.  

---

## Requirements

- Python 3.x
- [`tqdm`](https://pypi.org/project/tqdm/)
- [`python-Levenshtein`](https://pypi.org/project/python-Levenshtein/)

Install dependencies using pip (alternatively, you can use conda). 

```bash
pip install tqdm python-Levenshtein
```

---

### Input Files

#### 1. FASTQ File (`--fastq`)
Input file containing sequencing reads.
- **Format**: Standard 4-line FASTQ format per read.


#### 2. Design CSV File (`--design_csv`)
CSV file describing predefined barcode (index) and sequence pairs.
Must contain the following columns:
  - `barcode`: A short identifying sequence (used for matching)
  - `sequence`: The expected full insert sequence


---

### Output 

#### Matched Results CSV (`--output_csv`)
- **Description**: Output CSV file containing matched, trimmed sequences and their associated variant IDs.
- **Format**: CSV with the following columns:
  - `enum`: Enumeration of the read (used just for convenience)
  - `sequence`: The  DNA read  after trimming (the part that was found between the two primers)
  - `count`: Number of times this exact read was observed in the fastq file
  - `variant_id`: Index of the best-matching barcode from the design file (or `-1` if no acceptable match found)
- **Default**: `matched_results.csv`



### Notes
#### If primers (adapters) are not found (in either direction), the read is filtered out and is not represented in the output.

#### Barcodes (Indices) are only compared if a trimmed read is successfully extracted.

#### variant_id is the index of the best-matching barcode, as it appears in the design.csv; -1 if no match found.




