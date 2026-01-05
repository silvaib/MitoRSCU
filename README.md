# MitoRSCU
MitoRSCU is a tool designed to conduct Relative Synonymous Codon Usage (RSCU) on mitochondrial protein-coding sequences.
It takes fasta files with protein-coding gene sequences as inputs, and conduct an RSCU analysis and generate a clustered heatmap, 
aiming for a better visualization when working with numerous species.



[![DOI](https://zenodo.org/badge/1116780578.svg)](https://doi.org/10.5281/zenodo.17938892)


## Prerequisites
```bash
python=3.14.1
biopython=1.85
matplotlib=3.10.7
scipy=1.16.3
seaborn=0.13.2
```

## Installation

### Clone the project

```bash
  git clone https://github.com/silvaib/MitoRSCU.git
```

### Create mamba environment
```bash
  cd MitoRSCU
```
```bash
  mamba env create -f environment.yml
```

## Usage
The MitoRSCU.py takes a fasta sequence as a required input. You can pass more than one fasta file using the wildcard * or by directly indicating more than one file. 
The default genetic code is 'Vertebrate Mitochondrial', but any genetic code can be used. 
The genetic code must be declared either by the full name or the respective number. You can use --help to see more information.
### Recommendations
We suggest that you use the taxon name as the fasta file filename, since it will be used as default for naming the sequences. 

Start the env
```bash
  mamba activate mitorscu
```
Basic usage
```bash
  python MitoRSCU.py --input protein_coding_sequences.fasta --genetic_code 'Vertebrate Mitochondrial' --output rscu_analysis_output
  python MitoRSCY.py -i seq1.fasta seq2.fasta seq3.fasta --genetic_code 'Invertebrate Mitochondrial'
  python MitoRSCU.py -i *.fasta -g 2 -o rscu_analysis_output
```
#### Example data
The example dataset includes the 13 mitochondrial protein‑coding genes from several vertebrate model organisms.

```bash
 python MitoRSCU.py --input example_data/*.fasta --genetic_code 2 --output example_data_results
```
##### Clustered Heatmap
![RSCU_plot_clustered]

## License

Distributed under the GPL-3.0 license. See LICENSE.txt for more information.


## Citation

If you used this tool, we ask you to cite: https://doi.org/10.5281/zenodo.17938892
