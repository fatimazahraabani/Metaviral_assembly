#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# For parse table
import pandas as pd
# Import click for argument parser
import rich_click as click
from collections import defaultdict
from pathlib import Path
# Import orffinder package for predict orf
from orffinder import orffinder

def fasta2dict(filename):
    """
    This function take a file path (fasta), and return a dictionary of sequence
    """
    with open(filename, "r") as fastaFile:
        return SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))

def search_orf_orffinder(biopython_object, length_min=300):
    """
    This function take a biopython sequence in input and create a dict object which contains all information about
    orf predict by orrfinder package (near that ncbi result).
    The minimum length of orf is by default to 75, you can change this with the length_min=XX option.
    """
    # Create the orffinder object with biopython object
    orf_result = orffinder.getORFs(biopython_object, minimum_length=length_min)
    # If the tools not detect orf, the function return a None object
    if not orf_result:
        return None
    # Create liste of gene include in other gene which must be delete
    liste_to_remove = []

    # Compare orf with previously good orf for see if orf is inclued in the previously good orf

    # For frame "+"
    if '+' in [elt["sense"] for elt in orf_result]:
        orf_result_sens_pos = [elt for elt in orf_result if elt['sense'] == '+']
        # Sort result in function of orf start
        orf_result_sens_pos = sorted(orf_result_sens_pos, key=lambda x: x['start'])
        start_orf_prev = orf_result_sens_pos[0]['start']
        end_orf_prev = orf_result_sens_pos[0]['end']
        for orf in orf_result_sens_pos[1:]:
            if orf['start'] >= start_orf_prev and orf['end'] <= end_orf_prev:
                liste_to_remove.append(orf)  # If orf is inclued in previously orf, he is remove or result
            else:  # If orf is not inclued, we change the previously start and end orf by this one
                start_orf_prev = orf['start']
                end_orf_prev = orf['end']

    # For frame "-"
    if '-' in [elt["sense"] for elt in orf_result]:
        orf_result_sens_neg = [elt for elt in orf_result if elt['sense'] == '-']
        # Sort result in function of orf start
        orf_result_sens_neg = sorted(orf_result_sens_neg, key=lambda x: x['start'])
        start_orf_prev = orf_result_sens_neg[0]['start']
        end_orf_prev = orf_result_sens_neg[0]['end']
        for orf in orf_result_sens_neg[1:]:
            if orf['start'] >= start_orf_prev and orf['end'] <= end_orf_prev:
                liste_to_remove.append(orf)  # If orf is inclued in previously orf, he is remove or result
            else:  # If orf is not inclued, we change the previously start and end orf by this one
                start_orf_prev = orf['start']
                end_orf_prev = orf['end']

    # Remove all orf inclued in an another orf
    liste_final = [orf for orf in orf_result if orf not in liste_to_remove]

    # Create dico result for other function
    dico = {}
    for i, orf in enumerate(liste_final):
        id_orf = f"{biopython_object.id}_g{i + 1}"
        dico[id_orf] = {"start": orf['start'],  # Start of ORF
                        "end": orf['end'],  # End of ORF
                        "end_partial": orf['trailing'],  # True if end of ORF is truncated
                        "partial": orf['trailing'],  # partial is True if start or end is truncated
                        "strand": orf['sense']}
    return dico

@click.command("make_source_file", short_help=f'Install snakevir on HPC cluster',
               context_settings=dict(max_content_width=800))
@click.option('--fasta_path', '-f', type=click.Path(exists=True, resolve_path=True), required=True,
              help="Path of fasta file that contains all contigs to submit to GenBank")
@click.option('--output_path', '-o', type=click.Path(exists=False, resolve_path=True), required=True,
              help="Path of blastn result with the option '6 qseqid sseqid qlen slen "
                   "length qstart qend sstart send qcovhsp pident evalue bitscore'")
def main(fasta_path, output_path):
    """
    """
    fasta = fasta2dict(fasta_path)
    with open(output_path, 'w') as output:
        for seq_id in fasta:
            dico_orf = search_orf_orffinder(fasta[seq_id])
            # Create fasta with all orf find
            if dico_orf != None:
                for orf_id in dico_orf:
                    orf = dico_orf[orf_id]
                    if orf['strand'] == '+':
                        sequence = fasta[seq_id].seq[orf['start'] - 1:orf['end'] - 1]
                    if orf['strand'] == '-':
                        sequence = fasta[seq_id].seq[orf['end'] - 1:orf['start'] - 1].reverse_complement()
                    record = SeqRecord(sequence.translate(), id=orf_id, name=orf_id, description='')
                    SeqIO.write(record, output, "fasta")

if __name__ == '__main__':
    main()