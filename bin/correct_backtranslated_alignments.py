#! /usr/bin/env python

# Author: Sam Nooij
# Date: 2-2-2018
# Email: sam.nooij [at] rivm [dot] nl

# Script to check (reverse translated) multiple alignment
# files and correct the number of gaps to make all sequences
# equally long.

from Bio import SeqIO
from sys import argv
from operator import itemgetter

try:
    input_fasta = argv[1]
    fasta_extensions = ["fasta", "fa", "fas", "fna"]
    assert input_fasta.split('.')[-1] in fasta_extensions
except:
    print("Please provide a fasta file as argument.")
    print("For example: $ python correct_backtranslated_alignments.py alignment.fasta")
    exit()

def collect_sequences(input_fasta):
    """
    Read the fasta file and save all the sequences in a
    list, as tuple of (sequence, length)
    """
    fasta_file = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        fasta_file.append((seq_record, len(seq_record.seq)))
        
    return fasta_file

def correct_lengths(fasta_file):
    """
    Uses a list of seq records + sequence lengths to
    determine the required length and fill up too short
    sequences.
    """
    max_length = max(fasta_file, key=itemgetter(1))[1]
    
    for seq_record in fasta_file:
        if len(seq_record[0].seq) < max_length:
            difference = max_length - len(seq_record[0].seq)
            seq_record[0].seq += difference * "-"
        else:
            pass

    #Remove the lengths from the SeqRecords list
    fasta_file = [item[0] for item in fasta_file]

    return fasta_file


if __name__ == "__main__":
    fasta_file = correct_lengths(collect_sequences(input_fasta))
    
    #Change the name of the output file if desired:
    #output_fasta = input_fasta.split('.')[0] + "_corrected." + input_fasta.split('.')[-1]

    #SeqIO.write(fasta_file, output_fasta, "fasta")
    
    #Else use the same name:
    SeqIO.write(fasta_file, input_fasta, "fasta")
