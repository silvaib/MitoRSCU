from Bio import SeqIO
from Bio.Data import CodonTable
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import collections
import argparse
import logging
import csv
import io
import os


parser = argparse.ArgumentParser(description=r"""

                            ,--.   ,--. ,--.   ,--.           ,------.   ,---.    ,-----. ,--. ,--.
                            |   `.'   | `--' ,-'  '-.  ,---.  |  .--. ' '   .-'  '  .--./ |  | |  |
                            |  |'.'|  | ,--. '-.  .-' | .-. | |  '--'.' `.  `-.  |  |     |  | |  |
                            |  |   |  | |  |   |  |   ' '-' ' |  |\  \  .-'    | '  '--'\ '  '-'  '
                            `--'   `--' `--'   `--'    `---'  `--' '--' `-----'   `-----'  `-----'

MitoRSCU is a tool designed to conduct Relative Synonymous Codon Usage (RSCU) on protein-coding sequences of mitochondrial genomes.
It takes fasta files with protein-coding gene sequences as inputs, and conduct an RSCU analysis and generate a clustered heatmap, 
aiming for a better visualization when working with numerous species and analyzing RSCU.

Basic usage:
python MitoRSCU.py --input protein_coding_sequences.fasta --genetic_code 'Vertebrate Mitochondrial' --output rscu_analysis_output
python MitoRSCY.py -i seq1.fasta seq2.fasta seq3.fasta --genetic_code 'Invertebrate Mitochondrial'
python MitoRSCU.py -i *.fasta -g 2 -o rscu_analysis_output

Test data:
MitoRSCU.py --input example_data/*.fasta --genetic_code 2 --output rscu_analysis_output


""",

                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--input', '-i', nargs='+',
                    help="""Fasta file with protein-coding genes. One species per file. Only DNA sequences are accepted.
                    Multiple files can be inputed by using *.fasta or indicating more than one file in this argument. 
                    We recommend using the taxon name as the fasta filename.""")
parser.add_argument('--output', '-o', default='rscu_analysis_output',
                    help='Output name (default = rscu_analysis_output)')
parser.add_argument('--genetic_code', '-g', default='Vertebrate Mitochondrial',
                    help="""The genetic code must be declared either by the full name or the respective number.
                    (deafult='Vertebrate Mitochondrial')
                    "1": "Standard",
                    "2": "Vertebrate Mitochondrial",
                    "3": "Yeast Mitochondrial",
                    "4": "Mold, Protozoan, and Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
                    "5": "Invertebrate Mitochondrial",
                    "6": "Ciliate, Dasycladacean and Hexamita Nuclear",
                    "9": "Echinoderm and Flatworm Mitochondrial",
                    "10": "Euplotid Nuclear",
                    "11": "Bacterial, Archaeal and Plant Plastid",
                    "12": "Alternative Yeast Nuclear",
                    "13": "Ascidian Mitochondrial",
                    "14": "Alternative Flatworm Mitochondrial",
                    "15": "Blepharisma Macronuclear",
                    "16": "Chlorophycean Mitochondrial",
                    "21": "Trematode Mitochondrial",
                    "22": "Scenedesmus obliquus Mitochondrial",
                    "25": "Candidate Division SR1 and Gracilibacteria",
                    "26": "Pachysolen tannophilus Nuclear",
                    "27": "Karyorelict Nuclear",
                    "28": "Condylostoma Nuclear",
                    "29": "Mesodinium Nuclear",
                    "30": "Peritrich Nuclear",
                    "31": "Blastocrithidia Nuclear",
                    "32": "Balanophoraceae Plastid",
                    "33": "Cephalodiscidae Mitochondrial"
""")
parser.add_argument('-v', '--version',
                    action='version',
                    version='MitoRSCU 0.1.0',
                    help="Show program's version number and exit.")

args = parser.parse_args()


aa_dict = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
           'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
           'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'O': 'Pyl', 'S': 'Ser', 'U': 'Sec',
           'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', 'B': 'Asx', 'Z': 'Glx',
           'X': 'XAA', 'J': 'Xle'}


GENETIC_CODES = {
    "1": "Standard",
    "2": "Vertebrate Mitochondrial",
    "3": "Yeast Mitochondrial",
    "4": "Mold, Protozoan, and Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
    "5": "Invertebrate Mitochondrial",
    "6": "Ciliate, Dasycladacean and Hexamita Nuclear",
    "9": "Echinoderm and Flatworm Mitochondrial",
    "10": "Euplotid Nuclear",
    "11": "Bacterial, Archaeal and Plant Plastid",
    "12": "Alternative Yeast Nuclear",
    "13": "Ascidian Mitochondrial",
    "14": "Alternative Flatworm Mitochondrial",
    "15": "Blepharisma Macronuclear",
    "16": "Chlorophycean Mitochondrial",
    "21": "Trematode Mitochondrial",
    "22": "Scenedesmus obliquus Mitochondrial",
    "25": "Candidate Division SR1 and Gracilibacteria",
    "26": "Pachysolen tannophilus Nuclear",
    "27": "Karyorelict Nuclear",
    "28": "Condylostoma Nuclear",
    "29": "Mesodinium Nuclear",
    "30": "Peritrich Nuclear",
    "31": "Blastocrithidia Nuclear",
    "32": "Balanophoraceae Plastid",
    "33": "Cephalodiscidae Mitochondrial"
}

genetic_code_input = str(args.genetic_code)
genetic_code_name = GENETIC_CODES.get(genetic_code_input, genetic_code_input
                                      )
if not os.path.exists(args.output):
    os.makedirs(args.output)


def setup_logging():
    """
    Setup logging configuration.
    """
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s- %(message)s',
                        handlers=[
                            logging.FileHandler(
                                f"{args.output}/rscu_analysis.log"),
                            logging.StreamHandler()
                        ]
                        )


def fasta_to_single_line(input_file):
    """
    Test if the input fasta file is a single-line fasta.
    If not, it converts it into a single-line.

    Returns:
    A single-line fasta as a StringIO object.
    """
    logging.info(f"Checking if the input is a single-line fasta.")

    if isinstance(input_file, list):
        input_file = input_file[0]
    output_file = io.StringIO()
    with open(input_file) as f_input:
        block = []
        for line in f_input:
            if line.startswith('>'):
                if block:
                    output_file.write(''.join(block) + '\n')
                    block = []
                output_file.write(line)
            else:
                block.append(line.strip())
        if block:
            output_file.write(''.join(block) + '\n')
    output_file.seek(0)
    return output_file


def sequence_preparation(pcgs_sequences, sample_dir, sample_name):
    """
    Reads the input fasta, calculates the length and if the sequence is not a
    multiple of three, adds a - chracter until this condition is satisfied and
    concatenates the sequences in a single line.

    Returns:
    A single line fasta file with concatenated sequences.
    """
    logging.info(
        f"Calculating the sequence(s) length and adjusting it if needed. Concatenating if there is more than one sequence.")

    sequences = []
    concat_file = os.path.join(sample_dir, f"{sample_name}.fasta")
    for record in SeqIO.parse(pcgs_sequences, "fasta"):
        seq = str(record.seq)
        while len(seq) % 3 != 0:
            seq += '-'
        sequences.append(seq)
    concat_seq = ''.join(sequences)
    with open(concat_file, "w") as output_handle:
        output_handle.write(f">{sample_name}\n{concat_seq}")
    type(concat_file)

    logging.info(f"Calculation and concatenation are finished.")

    return concat_file


def rscu_analysis(concatenated_sequences, genetic_code):
    """_summary_

    Args:
        concatenated_sequences (fasta file)
        genetic_code (string)

    Returns:
    Two tables: counts_file and codon_table
    """
    logging.info(
        f"Starting RSCU analysis with the genetic code: {genetic_code_name}")

    nested_rscu = {}
    gen_code = CodonTable.unambiguous_dna_by_name[genetic_code]
    codon_counts = collections.Counter()

    all_codons = gen_code.forward_table.keys()
    for codon in all_codons:
        codon_counts[codon] = 0

    for record in SeqIO.parse(concatenated_sequences, "fasta"):
        codon_counts.update(str(record.seq[i:i+3])
                            for i in range(0, len(record.seq), 3))

    counts_file = io.StringIO()
    writer = csv.writer(counts_file)
    writer.writerow(['Codon', 'Count'])

    for codon, count in codon_counts.items():
        writer.writerow([codon, count])

    codon_table = io.StringIO()
    codon_table_writer = csv.writer(codon_table)
    codon_table_writer.writerow(['Codon', 'Amino_acid', 'RSCU'])

    for codon, aa in gen_code.forward_table.items():
        syn_codons = [c for c, a in gen_code.forward_table.items() if a == aa]
        observed_count = codon_counts[codon]
        expected_count = sum(codon_counts[c]
                             for c in syn_codons) / len(syn_codons)
        if expected_count == 0:
            rscu_value = 0
        else:
            rscu_value = observed_count / expected_count

        aa_rscu = {}
        aa_rscu[aa] = rscu_value
        nested_rscu[codon] = aa_rscu
        codon_table_writer.writerow([codon, aa, rscu_value])

    counts_file.seek(0)
    codon_table.seek(0)

    logging.info(f"RSCU analysis is finished.")

    return counts_file, codon_table


def csv_adjusts(counts, rscu, sample_dir, sample_name):
    logging.info(f"Processing results.")
    df_codon_counts = pd.read_csv(counts)
    df_codon_table = pd.read_csv(rscu)

    # Merge the two DataFrames on the 'Codon' column
    df_merged = pd.merge(df_codon_counts, df_codon_table, on='Codon')
    df_sorted = df_merged.sort_values('Amino_acid')
    df_sorted['Amino_acid'] = df_sorted['Amino_acid'].map(aa_dict)
    df_sorted = df_sorted.reindex(
        columns=['Amino_acid', 'Codon', 'Count', 'RSCU'])

    # Save the sorted dataframe as a csv file
    df_sorted.to_csv(os.path.join(
        sample_dir, f'{sample_name}_RSCU_results.csv'), index=False)


def process_fasta_file(fasta, output_dir, genetic_code):
    """Process a  FASTA file and perform RSCU analysis.

    Args:
        fasta (str): Path to the input FASTA file.
        output_dir (str): Path to the output directory.
        genetic_code (str): Name or number of the genetic code.

    Returns:
        None
    """
    try:
        sample_name = os.path.splitext(os.path.basename(fasta))[0]
        sample_dir = os.path.join(output_dir, sample_name)
        os.makedirs(sample_dir, exist_ok=True)

        if not os.path.exists(fasta):
            raise FileNotFoundError(f"Input file not found: {fasta}")

        logging.info(f"Processing {fasta}...")
        single_line = fasta_to_single_line(fasta)
        concat_seq = sequence_preparation(single_line, sample_dir, sample_name)
        rscu = rscu_analysis(concat_seq, genetic_code)
        csv_adjusts(rscu[0], rscu[1], sample_dir, sample_name)
        logging.info(f"Processing is finished.")

    except Exception as e:
        logging.error(f"Error processing {fasta}: {e}")


def read_and_concat_results(results_dir):
    """
    Read and concatenates all CSV files found in subdirectories of results folder.
    Args:
        results_dir (str): Directory containing subdirectories with CSV files from the individual analysis
        output_file (str): combined_df, species_names
    """
    all_results_data = []
    species_names = []

    for dirpath, _, filenames in os.walk(results_dir):
        for file in filenames:
            if file.endswith('.csv'):
                species_name = os.path.basename(dirpath).replace("_", " ")
                csv_path = os.path.join(dirpath, file)
                data = pd.read_csv(csv_path)
                data['Species'] = species_name
                all_results_data.append(data)
                if species_name not in species_names:
                    species_names.append(species_name)
    if all_results_data:
        combined_df = pd.concat(all_results_data, ignore_index=True)
        combined_df.to_csv(os.path.join(results_dir,
                                        f"combined_RSCU_{args.output}.csv",), index=False)
        cols = ['Species'] + \
            [col for col in combined_df.columns if col != 'Species']
        combined_df = combined_df[cols]
        print(
            f"Successfully concatenated RSCU results")
        return combined_df, species_names
    else:
        print("No CSV files found in the specified directory")
        return pd.DataFrame(), []


def create_heatmap(data):
    """
    Process the combine_df to plot a clustered heatmap using the RSCU values.
    Args:
        data (str): Combined dataframe of RSCU values
        output_file (str): Heatmap plot
    """
    logging.info(f"Getting the necessary information for the Heatmap.")

    results_data = data[0]
    results_data['Amino_acid_Codon'] = results_data['Amino_acid'] + \
        '_' + results_data['Codon']
    heatmap_data = results_data.pivot_table(index='Species', columns=[
        'Amino_acid', 'Codon'], values='RSCU')
    max_rscu = results_data['RSCU'].max()

    logging.info(f"Starting clustered heatmap of the RSCU values...")

    g = sns.clustermap(heatmap_data, cmap="plasma", vmin=0, vmax=max_rscu,
                       linewidths=0.5, linecolor='black',
                       cbar_kws={'label': 'RSCU value', 'shrink': 0.2},
                       tree_kws=dict(linewidths=1.0),
                       figsize=(60, 30), annot=False)
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=90, ha='center', fontsize=25)
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), rotation=0, ha='left', fontsize=30, fontstyle='italic')
    g.ax_heatmap.xaxis.tick_bottom()
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.figure.axes[-1].set_ylabel('RSCU value', size=25)
    g.cax.tick_params(labelsize=20)

    plt.savefig(f"{args.output}/RSCU_plot_clustered.svg",
                format='svg', bbox_inches='tight')
    logging.info(r"""

        RSCU heatmap plot is ...      
        
        ______ _____ _   _ _____  _____ _    _ ______ _____  
        |  ____|_   _| \ | |_   _|/ ____| |  | |  ____|  __ \ 
        | |__    | | |  \| | | | | (___ | |__| | |__  | |  | |
        |  __|   | | | . ` | | |  \___ \|  __  |  __| | |  | |
        | |     _| |_| |\  |_| |_ ____) | |  | | |____| |__| |
        |_|    |_____|_| \_|_____|_____/|_|  |_|______|_____/ 
        
        
        You have succesfully executed the MitoRSCU.
        
        """)


def main():
    """Main function to process all input FASTA files."""
    setup_logging()
    for fasta in args.input:
        process_fasta_file(fasta, args.output, genetic_code_name)
    create_heatmap(read_and_concat_results(args.output))


if __name__ == "__main__":
    main()
