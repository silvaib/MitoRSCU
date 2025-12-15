MitoRSCU is a tool designed to conduct Relative Synonymous Codon Usage (RSCU) on protein-coding sequences of mitochondrial genomes.
It takes fasta files with protein-coding gene sequences as inputs, and conduct an RSCU analysis and generate a clustered heatmap, 
aiming for a better visualization when working with numerous species and analyzing RSCU.


<!-- Getting Started -->
## Getting Started

<!-- Prerequisites -->
### Prerequisites
python=3.14.1
biopython = 1.85
matplotlib = 3.10.7
scipy = 1.16.3
seaborn = 0.13.2


<!-- Installation -->
### Installation

Clone the project

```bash
  git clone https://github.com/silvaib/MitoRSCU.git
```

Go to the project directory

```bash
  cd MitoRSCU
```

Create conda/mamba environment

Conda
```bash
  conda env create -f environment.yml
```

Mamba
```bash
  mamba env create -f environment.yml
```
Start the env

Conda
```bash
  conda activate mitorscu
```
Mamba
```bash
  mamba activate mitorscu
```


<!-- Usage -->
The MitoRSCU.py takes a fasta sequence as a required input. You can pass more than one fasta file using the wildcard * or by directly indicating more than one file. 
The default genetic code is 'Vertebrate Mitochondrial', but any genetic code can be used. 
The genetic code must be declared either by the full name or the respective number. You can use --help to see more information.

```bash
  python MitoRSCU.py --input protein_coding_sequences.fasta --genetic_code 'Vertebrate Mitochondrial' --output rscu_analysis_output
  python MitoRSCY.py -i seq1.fasta seq2.fasta seq3.fasta --genetic_code 'Invertebrate Mitochondrial'
  python MitoRSCU.py -i *.fasta -g 2 -o rscu_analysis_output
```


<!-- License -->
## License

Distributed under the GPL-3.0 license. See LICENSE.txt for more information.


<!-- Citation -->
If you used this tool, we ask you to cite:
