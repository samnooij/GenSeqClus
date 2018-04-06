"""
Date: 2018-01-04
Author: Sam Nooij

Title: Sapovirus Cluster Analysis

This snakemake workflow analyses clusters in Sapovirus
genomic sequences. The clustering is based on three (3)
different distance measures: similarity (BLAST identity),
composition (NVR, like in Yu et al., 2013), and phylogeny 
(DEmARC, like in Lauber & Gorbalenya, 2012).
Therefore, first distances are calculated by running tBLASTn,
a custom-made NVR script (in R), and phylogeny using a
multiple sequence alignment from MAFFT and Maximum Likelihood
phylogeny from IQTree.
(N.B. Three multiple alignments are generated: full protein
sequences, back-translated coding sequences, and trimmed
protein sequences (with Gblocks, relaxed settings, like
in Talavera & Castresana, 2007).
The resulting distance values are imported into R and 
analysed using custom R code, visualised using RMarkdown.
"""

SAMPLES = ["SaV_genomes", "SaV_VP1"]

localrules: calculate_nvr_distances, cluster_analysis

rule all:
    input:
        expand("results/{sample}-Cluster_Analysis.html", sample = SAMPLES)

### BLAST ###

rule make_blastdb:
    input:
        "data/{sample}_nt.fasta"
    output:
        "tmp/{sample}_nt.nsq"
    shell:
        "makeblastdb -dbtype nucl -in {input} -out tmp/{wildcards.sample}_nt"

rule run_tblastn:
    input:
        fasta="data/{sample}_aa.fasta",
        # tblastn requires amino acid sequences as query
        db="tmp/{sample}_nt.nsq"
    params:
        hsps="-max_hsps 1",
        output_format="-outfmt 6"
    output:
        "tmp/{sample}_tblastn_aa.blast"
    shell:
        """
tblastn -db tmp/{wildcards.sample}_nt \
-query {input.fasta} \
{params.output_format} {params.hsps} \
-out {output}
        """
        

### NVR ###

rule calculate_nvr_distances:
    input:
        "data/{sample}_nt.fasta"
    output:
        "tmp/{sample}_nt.nvr"
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/NVR.R"

### kmer frequencies ###

import os
BASE_PATH = os.path.realpath("")

rule calculate_kmer_frequencies:
    input:
        "data/{sample}_nt.fasta"
    params:
        kmers=[1, 2, 3, 4, 5],
        wd=BASE_PATH
    output:
        [''.join(["tmp/{sample}_nt_", str(kmer), "mers.csv"]) for kmer in [1, 2, 3, 4, 5]]
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/kmer_frequencies.R"

### ML PHYLOGENY ###

rule make_multiple_alignment:
    input:
        "data/{sample}_aa.fasta"
    #make multile alignments of protein sequences
    output:
        "tmp/{sample}-mafft_aa.fas"
    shell:
        "mafft --auto {input} > {output}"

rule download_revtrans:
    input:
        "bin/download_revtrans.sh"
    output:
        "bin/revtrans.py"
    shell:
        "bash {input}"

rule download_gblocks:
    input:
        "bin/download_gblocks.sh"
    output:
        "bin/Gblocks"
    shell:
        "bash {input}"

rule backtranslate:
    input:
        bin="bin/revtrans.py",
        cor="bin/correct_backtranslated_alignments.py",
        nt="data/{sample}_nt.fasta",
        aa="tmp/{sample}-mafft_aa.fas"
    output:
        "tmp/{sample}-mafft-RevTrans_aa.fas"
    conda:
        "envs/cluster_analysis_py27.yaml"
    shell:
        """
        python {input.bin} -match name {input.nt} {input.aa} > {output}
        python {input.cor} {output}
        """

#Sometimes, RevTrans returns a multiple alignment with 
# sequences of differing lengths. The second script adds
# gaps "-" to the end of shorter sequences to fix this.

rule clean_multiple_alignment:
    input:
        bin="bin/Gblocks",
        fas="tmp/{sample}-mafft_aa.fas"
    params:
        relaxed="-b1=37 -b2=37 -b3=10 -b4=5 -b5=h"
    output:
        "tmp/{sample}-mafft_aa.fas-gb"
    shell:
        """
{input.bin} {input.fas} {params.relaxed} || true
        """
#Gblocks returns exit status 1 by default, which stops
# snakemake. Adding "|| true" should reset the exit status to 0.
        
rule realign_cleaned_multiple_alignment:
    input:
        "tmp/{sample}-mafft_aa.fas-gb"
    output:
        "tmp/{sample}-mafft-Gblocks-mafft_aa.fas"
    shell:
        "mafft --auto {input} > {output}"

rule infer_phylogeny_protein:
    input:
        "tmp/{sample}-mafft_aa.fas"
    params:
        test="-m TEST",
        bootstrap="-bb 1000"
    output:
        "tmp/{sample}-mafft_aa.fas.mldist"
    threads: 4
    shell:
        "iqtree {params.test} {params.bootstrap} -nt {threads} -s {input}"

rule infer_phylogeny_cds:
    input:
        "tmp/{sample}-mafft-RevTrans_aa.fas"
    params:
        test="-m TEST",
        bootstrap="-bb 1000"
    output:
        "tmp/{sample}-mafft-RevTrans_aa.fas.mldist"
    threads: 4
    shell:
        "iqtree {params.test} {params.bootstrap} -nt {threads} -s {input}"

### CLUSTER ANALYSIS ###

#Install one required R package that is not in conda:
#R: "install.packages("directlabels", repo="http://r-forge.r-project.org")"

rule cluster_analysis:
    input:
        "tmp/{sample}_tblastn_aa.blast",
        "tmp/{sample}_nt.nvr",
        "tmp/{sample}-mafft_aa.fas.mldist",
        "tmp/{sample}-mafft-RevTrans_aa.fas.mldist",
        "tmp/{sample}-mafft-Gblocks-mafft_aa.fas.mldist",
        "data/{sample}_metadata.tsv"
    output:
        "results/{sample}-Cluster_Analysis.html"
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/Cluster_analysis.Rmd"
