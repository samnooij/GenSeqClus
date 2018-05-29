"""
Date: 2018-01-04
Author: Sam Nooij

Title: Sapovirus Cluster Analysis

This snakemake workflow analyses clusters in Sapovirus
genomic sequences. The clustering is based on three (3)
different distance measures: similarity (BLAST identity),
composition (NVR, like in Yu et al., 2013; and oligomer
frequency-distances: 1-5mers), and phylogeny 
(DEmARC, like in Lauber & Gorbalenya, 2012).
Therefore, distances are calculated by running tBLASTn,
custom-made R scripts (for NVR, and the kmer frequencies), 
and phylogeny using a multiple sequence alignment from MAFFT
 and Maximum Likelihood phylogeny from IQTree.
(N.B. Three multiple alignments are generated: full protein
sequences, back-translated coding sequences, and trimmed
protein sequences (with Gblocks, relaxed settings, like
in Talavera & Castresana, 2007).)
The resulting distance values are imported into R and 
analysed using custom R code, visualised using RMarkdown.
"""

SAMPLES = ["SaV_genomes", "SaV_VP1"]
KMERS = [1, 2, 3, 4, 5]

localrules: calculate_nvr_distances, cluster_analysis

rule all:
    input:
        expand("results/{sample}-Cluster_Analysis.html", sample = SAMPLES)

### BLAST --------------------------------------------------

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
        

### NVR ----------------------------------------------------

rule calculate_nvr:
    input:
        "data/{sample}_nt.fasta"
    output:
        "tmp/{sample}_nt.nvr"
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/NVR.R"

rule convert_nvr_to_distances:
    input:
        "tmp/{sample}_nt.nvr"
    output:
        matrix="tmp/{sample}_nt_nvr_distances.RDS",
        dataframe="tmp/{sample}_nt_nvr_distances.csv"
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/euclidean_distances.R"


### kmer frequencies ---------------------------------------

rule calculate_kmer_frequencies:
    input:
        "data/{sample}_nt.fasta"
    params:
        kmers=KMERS
    output:
        [''.join(["tmp/{sample}_nt_", str(kmer), "mers.csv"]) for kmer in KMERS]        
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/kmer_frequencies.R"

rule convert_kmers_to_distances:
    input:
        [''.join(["tmp/{sample}_nt_", str(kmer), "mers.csv"]) for kmer in KMERS]        
    output:
        matrix=[''.join(["tmp/{sample}_nt_", str(kmer), "mer_distances.RDS"]) for kmer in KMERS],
        dataframe=[''.join(["tmp/{sample}_nt_", str(kmer), "mer_distances.csv"]) for kmer in KMERS]
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/euclidean_distances.R"


### ML PHYLOGENY -------------------------------------------

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

#Sometimes RevTrans returns a multiple alignment with 
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

### CLUSTER ANALYSIS ---------------------------------------

#Install one required R package that is not in conda:
#R: "install.packages("directlabels", repo="http://r-forge.r-project.org")"

rule cluster_analysis:
    input:
        blast="tmp/{sample}_tblastn_aa.blast",
        nvr_mat="tmp/{sample}_nt_nvr_distances.RDS",
        nvr_df="tmp/{sample}_nt.nvr",
        kmer_dfs=expand("{sample}_nt_{kmer}mer_distances.csv", kmer = KMERS, sample = SAMPLES),
        kmer_mats=expand("{sample}_nt_{kmer}mer_distances.RDS", kmer = KMERS, sample = SAMPLES),
        aa_ml="tmp/{sample}-mafft_aa.fas.mldist",
        nt_ml="tmp/{sample}-mafft-RevTrans_aa.fas.mldist",
        aa_ml_clean="tmp/{sample}-mafft-Gblocks-mafft_aa.fas.mldist",
        metadata="data/{sample}_metadata.tsv"
    output:
        blast="tmp/{sample}_blast_clustering_results.csv"
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/cluster_annotated_sequences.R"

rule create_cluster_report:
    input:
        ""
#all the different "_clustering_results.csv" files
    output:
        ".html"
#separate reports? one big report?
    conda:
        "envs/cluster_analysis_py27.yaml"
    script:
        "bin/Cluster_analysis_visualisation.Rmd"
