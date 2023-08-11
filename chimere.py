#!/usr/bin/env python3
# For parse table
import pandas as pd
# Import click for argument parser
import rich_click as click
from collections import defaultdict
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.backends.backend_pdf

def parse_multi_hit_blast(path_file):
    """
    """
    data = pd.read_csv(path_file, sep='\t')
    dico_aln = defaultdict(dict)
    for qseqid in data['qseqid']:
        data_qseqid = data.loc[data['qseqid'] == qseqid]
        for sseqid in data_qseqid['sseqid'] :
            data_qseqid_sseqid = data_qseqid.loc[data_qseqid['sseqid'] == sseqid]
            if len(data_qseqid_sseqid) > 2:
                for index in data_qseqid_sseqid.index:
                    contig = data_qseqid_sseqid.loc[index,'qseqid']
                    len_aln = data_qseqid_sseqid.loc[index,'send'] - data_qseqid_sseqid.loc[index,'sstart']
                    # Here we create a dict with length alignement, for the frame (sens) we use this equation for keep only -1 for reverse frame ans 1 for positive frame.
                    dico_aln[contig][abs(len_aln)] = {"sseqid": sseqid.split('|')[-2], "pident":data_qseqid_sseqid.loc[index,'pident'], 'sens' : len_aln/abs(len_aln),'qlen':data_qseqid_sseqid.loc[index,'qlen'], 'start' : data_qseqid_sseqid.loc[index,'qstart'], 'end': data_qseqid_sseqid.loc[index,'qend']}
    return(dico_aln)

def merge_block(dico_aln,contig):
    '''
    '''
    sens = 0
    new_dict_aln = defaultdict(dict)
    for elt in sorted(dico_aln[contig], key=lambda x:dico_aln[contig][x]['start']):
        if dico_aln[contig][elt]["pident"] < 80 :
            continue
        if dico_aln[contig][elt]["sens"] == sens :
            new_len = dico_aln[contig][elt]['end'] - old_start
            new_dict_aln[contig][new_len] = dico_aln[contig][elt].copy()
            new_dict_aln[contig][new_len]['start'] = old_start
            if old_elt in new_dict_aln[contig] :
                if old_elt < new_len :
                    del new_dict_aln[contig][old_elt]
                elif old_elt != new_len :
                    del new_dict_aln[contig][new_len]
        else :
            new_dict_aln[contig][elt] = dico_aln[contig][elt]
            old_start = dico_aln[contig][elt]["start"]
            old_elt = elt
        sens = dico_aln[contig][elt]["sens"]
    return(new_dict_aln)


def trim_contig(dico, contig):
    """
    """
    dico_aln = merge_block(dico, contig)
    max_aln = max(dico_aln[contig])
    max_aln_position = dico_aln[contig][max_aln]

    trim_start, trim_end = 1, dico_aln[contig][max_aln]['qlen']  # At start we keep all contig !
    for len_aln in dico_aln[contig]:
        if max_aln >= 0.33 * dico_aln[contig][max_aln]['qlen']:  # Check if blast hit cover 1/3 of contig
            if dico_aln[contig][max_aln]['sens'] != dico_aln[contig][len_aln]['sens']:
                if dico_aln[contig][len_aln]['start'] < dico_aln[contig][max_aln]['start']:
                    trim_start = dico_aln[contig][max_aln]['start']  # But if they have gap between end and start ?
                if dico_aln[contig][len_aln]['end'] > dico_aln[contig][max_aln]['end']:
                    trim_end = dico_aln[contig][max_aln]['end']  # But if they have gap between end and start ?
    return (trim_start, trim_end)

def fig_blast_hit(dico_aln,contig_test):
    """
    This function trace blast
    """
    max_aln = max(dico_aln[contig_test])
    qlen = int(dico_aln[contig_test][max_aln]['qlen'])
    good_strand = int(dico_aln[contig_test][max_aln]['sens'])
    dict_sseqid = {dico_aln[contig_test][elt]['sseqid']: 0 for elt in dico_aln[contig_test]}
    for elt in dico_aln[contig_test]:
        if float(dico_aln[contig_test][elt]['pident']) > 80:
            dict_sseqid[dico_aln[contig_test][elt]['sseqid']] += int(elt)
    sseqid_long = [ sseqid for sseqid,len in dict_sseqid.items() if len == max(dict_sseqid.values()) ][0]
    features=[GraphicFeature(start=int(dico_aln[contig_test][elt]['start']),
                             end=int(dico_aln[contig_test][elt]['end']),
                             strand=int(dico_aln[contig_test][elt]['sens']),
                             color= f'{"#0131b4" if int(dico_aln[contig_test][elt]["sens"]) == good_strand else "#fb1919"}',
                             label=str(dico_aln[contig_test][elt]['sseqid'])) for elt in dico_aln[contig_test] if (float(dico_aln[contig_test][elt]['pident']) > 80 and str(dico_aln[contig_test][elt]['sseqid']) == sseqid_long)]
                             # label=str(elt)) for elt in dico_aln[contig_test] if float(dico_aln[contig_test][elt]['pident']) > 80]

    record = GraphicRecord(sequence_length=qlen, features=features)
    fig,_ = record.plot(figure_width=15)
    start,end = trim_contig(dico_aln,contig_test)
    fig.axvline(x=start, color='red')
    fig.axvline(x=end, color='red')
    fig.set_title(f'{contig_test} - keep {start} to {end}')
    return(fig)

@click.command("chimere_detection", short_help=f'Install snakevir on HPC cluster',
               context_settings=dict(max_content_width=800))
@click.option('--blast_path', '-b', type=click.Path(exists=True, resolve_path=True), required=True,
              help="Path of fasta file that contains all contigs to submit to GenBank")
@click.option('--output_path', '-o', type=click.Path(exists=False, resolve_path=True), required=True,
              help="Path of blastn result with the option '6 qseqid sseqid qlen slen "
                   "length qstart qend sstart send qcovhsp pident evalue bitscore'")
def main(blast_path, output_path):
    dico_aln = parse_multi_hit_blast(blast_path)
    pdf = matplotlib.backends.backend_pdf.PdfPages(output_path)
    for contig in dico_aln.keys():
        fig = fig_blast_hit(dico_aln, contig)
        pdf.savefig(fig.figure)
    pdf.close()


if __name__ == '__main__':
    main()
