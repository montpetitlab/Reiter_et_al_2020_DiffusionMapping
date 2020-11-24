RAW_DATADIR = "inputs/raw_data"
SAMPLES = ["2017_AS1_1_1", "2017_AS1_1_2", "2017_AS1_1_3", "2017_AS1_1_4",
           "2017_AS1_4_1", "2017_AS1_4_2", "2017_AS1_4_3", "2017_AS1_4_4",
           "2017_AS2_5_1", "2017_AS2_5_2", "2017_AS2_5_3", "2017_AS2_5_4",
           "2017_AS2_6_1", "2017_AS2_6_2", "2017_AS2_6_3", "2017_AS2_6_4",
           "2017_SNC2_9_1", "2017_SNC2_9_2", "2017_SNC2_9_3", "2017_SNC2_9_4",
           "2017_SNC2_9_5", "2017_SNC2_12_1", "2017_SNC2_12_2", "2017_SNC2_12_3",
           "2017_SNC2_12_4", "2017_SNC2_12_5", "2017_CRN1_13_1", "2017_CRN1_13_2",
           "2017_CRN1_13_3", "2017_CRN1_13_4", "2017_CRN1_13_5", "2017_CRN1_16_1",
           "2017_CRN1_16_2", "2017_CRN1_16_3", "2017_CRN1_16_4", "2017_CRN1_16_5",
           "2017_SMV1_19_1", "2017_SMV1_19_2", "2017_SMV1_19_3", "2017_SMV1_19_4",
           "2017_SMV1_19_5", "2017_SMV1_20_1", "2017_SMV1_20_2", "2017_SMV1_20_3",
           "2017_SMV1_20_4", "2017_SMV1_20_5", "2017_SMV2_22_1", "2017_SMV2_22_2",
           "2017_SMV2_22_3", "2017_SMV2_22_4", "2017_SMV2_22_5", "2017_SMV2_24_1",
           "2017_SMV2_24_2", "2017_SMV2_24_3", "2017_SMV2_24_4", "2017_SMV2_24_5",
           "2017_RRV2_25_1", "2017_RRV2_25_2", "2017_RRV2_25_3", "2017_RRV2_25_4",
           "2017_RRV2_25_5", "2017_RRV2_28_1", "2017_RRV2_28_2", "2017_RRV2_28_3",
           "2017_RRV2_28_4", "2017_RRV2_28_5", "2017_RRV1_29_1", "2017_RRV1_29_2",
           "2017_RRV1_29_3", "2017_RRV1_29_4", "2017_RRV1_29_5", "2017_RRV1_31_1",
           "2017_RRV1_31_2", "2017_RRV1_31_3", "2017_RRV1_31_4", "2017_RRV1_31_5",
           "2017_SNC1_33_1", "2017_SNC1_33_2", "2017_SNC1_33_3", "2017_SNC1_33_4",
           "2017_SNC1_33_5", "2017_SNC1_34_1", "2017_SNC1_34_2", "2017_SNC1_34_3",
           "2017_SNC1_34_4", "2017_SNC1_34_5", "2017_SRH1_37_1", "2017_SRH1_37_2",
           "2017_SRH1_37_3", "2017_SRH1_37_4", "2017_SRH1_37_5", "2017_SRH1_40_1", 
           "2017_SRH1_40_2", "2017_SRH1_40_3", "2017_SRH1_40_4", "2017_SRH1_40_5",
           "2017_AV1_42_1", "2017_AV1_42_2", "2017_AV1_42_3", "2017_AV1_42_4",
           "2017_AV1_42_5", "2017_AV1_44_1", "2017_AV1_44_2", "2017_AV1_44_3",
           "2017_AV1_44_4", "2017_AV1_44_5", "2017_RRV3_46_1", "2017_RRV3_46_2",
           "2017_RRV3_46_3", "2017_RRV3_46_4", "2017_RRV3_46_5", "2017_RRV3_47_1", 
           "2017_RRV3_47_2", "2017_RRV3_47_3", "2017_RRV3_47_4", "2017_RRV3_47_5",
           "2017_AV2_50_1", "2017_AV2_50_2", "2017_AV2_50_3", "2017_AV2_50_4",
           "2017_AV2_50_5", "2017_AV2_51_1", "2017_AV2_51_2", "2017_AV2_51_3",
           "2017_AV2_51_4", "2017_AV2_51_5", "2017_OR1_53_1", "2017_OR1_53_2",
           "2017_OR1_53_3", "2017_OR1_53_4", "2017_OR1_53_5", "2017_OR1_54_1",
           "2017_OR1_54_2", "2017_OR1_54_3", "2017_OR1_54_4", "2017_OR1_54_5",
           "2017_OR2_59_1", "2017_OR2_59_2", "2017_OR2_59_3", "2017_OR2_59_4",
           "2017_OR2_59_5", "2017_OR2_60_1", "2017_OR2_60_2", "2017_OR2_60_3",
           "2017_OR2_60_4", "2017_OR2_60_5"] 

rule all:
    input:
       "outputs/counts/raw_counts.tsv",
       expand("outputs/gather/{sample}_gather.csv", sample = SAMPLES),
       "outputs/counts_all_organisms/raw_counts.tsv"

rule cat_fastq:
    output: 'inputs/cat/{sample}.fq.gz'
    params: indir = RAW_DATADIR
    shell:'''
    cat {params.indir}/{wildcards.sample}_S*_L00*_R1_001.fastq.gz > {output} 
    '''

rule trim_first_seven:
    output: 'outputs/quality-7/{sample}.trim7.fq.gz'
    input: 'inputs/cat/{sample}.fq.gz'
    conda: 'envs/env.yml'
    params: outdir='inputs/cat'
    shell:'''
    gunzip {input}
    fastx_trimmer -f 7 -z -i {params.outdir}/{wildcards.sample}.fq  -o {output}    
    gzip {params.outdir}/{wildcards.sample}.fq
    '''

rule split_TATA:
    output: 
        'outputs/quality-TATA/{sample}TATA1.fq.gz',
        'outputs/quality-TATA/{sample}unmatched.fq'
    input: 
        reads='outputs/quality-7/{sample}.trim7.fq.gz',
        barcode='inputs/tata.txt'
    conda: 'envs/env.yml'
    params:
        outdir = 'outputs/quality-TATA'
    shell:'''
    gunzip {input.reads}
    cat outputs/quality-7/{wildcards.sample}.trim7.fq | fastx_barcode_splitter.pl --bcfile {input.barcode} --bol --prefix {params.outdir}/{wildcards.sample} --suffix .fq
    gzip {params.outdir}/{wildcards.sample}TATA1.fq
    #gzip {params.outdir}/{wildcards.sample}unmatched.fq
    '''

rule trim_next_five:
    output: 
        twelve = 'outputs/quality-12/{sample}.trim12.fq.gz',
        un = 'outputs/quality-TATA/{sample}unmatched.fq.gz'
    input: 'outputs/quality-TATA/{sample}unmatched.fq'
    conda: 'envs/env.yml'
    params: outdir = 'outputs/quality-TATA'
    shell:''' 
    fastx_trimmer -f 5 -z -i {input} -o {output.twelve} 
    gzip {input}
    '''

rule bbduk_qc:
    output: 'outputs/quality-bb/{sample}.trimbb.fq.gz'
    input:
        fq = 'outputs/quality-12/{sample}.trim12.fq.gz',
        polyA = 'inputs/polya.fa',
        truseqr = 'inputs/truseq_rna.fa.gz',
    conda: 'envs/bbmap.yml'
    shell:'''
    bbduk.sh in={input.fq} out={output} ref={input.polyA},{input.truseqr} k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
    '''

rule download_genome:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
    '''

rule decompress_genome:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna'
    shell: "gunzip {input}"

rule download_gtf:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
    '''

rule decompress_gtf:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    shell: "gunzip {input}"

rule star_index_genome:
    input:
        genome = 'inputs/genome/GCF_000146045.2_R64_genomic.fna',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    params: input_dir = 'inputs/genome' 
    output: 'inputs/genome/SAindex'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir {params.input_dir} \
         --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang  99
    '''

rule gunzip_read:
    input: 
        reads = 'outputs/quality-bb/{sample}.trimbb.fq.gz'
    output: 
        reads = 'outputs/quality-bb/{sample}.trimbb.fq'
    shell:'''
    gunzip {input}
    '''

rule star_align:
    #v2.5.2a
    input:
        reads = 'outputs/quality-bb/{sample}.trimbb.fq',
        genome_index = 'inputs/genome/SAindex'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs/star/' + wildcards.sample,
        genome_dir = 'inputs/genome'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 2 --genomeDir {params.genome_dir}      \
        --readFilesIn {input.reads} --outFilterType BySJout  \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
        --alignIntronMax 1000000 --alignMatesGapMax 1000000  \
        --outSAMattributes NH HI NM MD --outSAMtype BAM      \
        SortedByCoordinate --outFileNamePrefix {params.out_prefix}
    '''

rule index_bam:
    input: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam.bai'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools index {input}
    '''
    
rule htseq_count:
    input:
        bam = 'outputs/star/{sample}Aligned.sortedByCoord.out.bam',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    output: "outputs/htseq/{sample}_readcounts.txt"
    conda: "envs/htseq.yml"
    shell:'''
    htseq-count -m intersection-nonempty -s yes -f bam -r pos {input.bam} {input.gtf} > {output}
    '''

rule make_counts:
    input: expand("outputs/htseq/{sample}_readcounts.txt", sample = SAMPLES)
    output: "outputs/counts/raw_counts.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/make_raw_counts.R"
