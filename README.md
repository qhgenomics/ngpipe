# ngpipe
 
ngpipe is a pipeline for typing _Neisseria gonorrhoeae_.
It takes a folder of reads or contigs and runs sequence typing for MLST, NG-STAR and NG-MAST typing.
It also determines species using the rplF gene and checks for the presence of rplF and ppng.
It will compile this information into a table and create a minimum spanning tree using MLST.

## Installation

Clone this repository

```git clone https://github.com/mjsull/ngpipe.git```

## Requirements
### Essential
#### Python
#### networkx
#### matplotlib
#### Snakemake - https://github.com/snakemake/snakemake
#### pyngoST
#### Minimap
#### Samtools
#### Shovill
#### pymlst


### Optional
#### Abricate
#### FastQC
#### QUAST

## Usage

#### To run ngpipe on a list of contigs
1. Create text file /path/to/samples.txt
2. Enter a list of samples to process, one on each line.
3. create directory /path/to/contigs
4. For each sample copy a FASTA of contigs to ```/path/to/contigs/<sample>.fasta```
5. __Run__

```snakemake --snakefile ngpipe/Snakefile --config workdir=/path/to/output contig_dir=/path/to/contigs  samples=/path/to/samples.txt  --cores 4```

#### To run ngpipe on a folder of Illumina reads
1. Create text file /path/to/samples.txt
2. Enter a list of samples to process, one on each line.
3. create directory /path/to/reads
4. For each sample copy read pairs to to ```/path/to/reads/<sample>_R1.fastq.gz``` and ```/path/to/reads/<sample>_R2.fastq.gz```
5. __Run__

```snakemake --snakefile ngpipe/Snakefile --config workdir=/path/to/output read_dir=/path/to/reads  samples=/path/to/samples.txt  --cores 4```



## Implementation

```mermaid
---
config:
  theme: redux
  layout: elk
---
flowchart TD
    A(["Start"]) --> B("Reads<br>(fastq)") & D("Contigs (fasta)")
    B --> C["<b>Assembly</b> <br> Shovill/<br>Spades"] & F["<b>Typing</b><br>pyMLST"] & I["<b>QC</b><br> fastQC"] & K["<b>Alignment</b><br>Minimap2"] & P["<b>Alignment</b><br>Minimap2"]
    C --> D
    D --> E["<b>Typing</b><br>pyngoST"] & F & G["<b>AMR prediction</b><br>Abricate"] & H["<b>QC</b><br> QUAST"]
    I --> J["<b>QC</b><br> MultiQC"]
    H --> J
    n6["rplF<br>(fasta)"] --> K
    n7["23S<br>(fasta)"] --> K
    n8["PPNG<br>(fasta)"] --> K
    K --> L("BAM")
    L --> M("Coverage information") & N("Minor alleles in 23s")
    n2["pubmlst"] --> n3["NG-MAST"] & n4["MLST"] & n5["rplF"]
    n1["ngstar.ca"] --> E & F
    n3 --> E & F
    n4 --> E & F
    n5 --> F
    E --> T["Merge typing results"]
    T --> O("Sequence types")
    F --> T
    O --> P & S("Final table")
    P --> Q("Locus coverage information")
    R("Version infomration") --> S
    J --> S
    Q --> S
    M --> S
    N --> S
    n2@{ shape: cyl}
    n1@{ shape: cyl}
    n6@{ shape: cyl}
    n7@{ shape: cyl}
    n8@{ shape: cyl}
    n3@{ shape: cyl}
    n4@{ shape: cyl}
    n5@{ shape: cyl}
    style A fill:#C8E6C9
    style S fill:#C8C8E6
```