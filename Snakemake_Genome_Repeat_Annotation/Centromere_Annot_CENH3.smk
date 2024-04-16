###########################################################################################################
rule centromere_regions_2_trim_galore:
    input:
        fq1 = "results/{prefix}/Centromere_CENH3/Chip-Seq_Data/{srarunid}_R1.fastq.gz",
        fq2 = "results/{prefix}/Centromere_CENH3/Chip-Seq_Data/{srarunid}_R2.fastq.gz"
    output:
        fq1 = "results/{prefix}/Centromere_CENH3/Trim_galore/{srarunid}/{srarunid}_R1.fastq.gz",
        fq2 = "results/{prefix}/Centromere_CENH3/Trim_galore/{srarunid}/{srarunid}_R2.fastq.gz"
    threads: 2
    log:
        "logs/{prefix}/Centromere_CENH3/Trim_galore_{srarunid}.log"
    conda:
        "../envs/Centromere_Annot.yaml"
    params:
        workdir = config["workdir"],
        samples = "{srarunid}",
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        trim_galore --paired --fastqc -j {threads} --gzip --max_n 50 \
          --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
          -o {params.outdir} \
          {input.fq1} {input.fq2} > {params.workdir}/{log} 2>&1
        
        mv {params.outdir}/{params.samples}_R1_val_1.fq.gz {output.fq1}
        mv {params.outdir}/{params.samples}_R2_val_2.fq.gz {output.fq2}
        '''

rule centromere_regions_3_bwa_index:
    input:
        genome = config["chr_genome"]
    output:
        multiext("results/{}/Centromere_CENH3/BWA_Alignment/{}".format(config["prefix"], os.path.basename(config["chr_genome"])), ".amb", ".ann", ".bwt", ".pac", ".sa"),
        "results/{}/Centromere_CENH3/BWA_Alignment/{}".format(config["prefix"], os.path.basename(config["chr_genome"]))
    params:
        prefix = os.path.basename(config["chr_genome"]),
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        module load BWA/0.7.17
        ln -s {input.genome} {params.outdir}
        cd {params.outdir}
        bwa index -p {params.prefix} {params.prefix}
        '''

rule centromere_regions_4_bwa_mem:
    input:
        fq1 = "results/{prefix}/Centromere_CENH3/Trim_galore/{srarunid}/{srarunid}_R1.fastq.gz",
        fq2 = "results/{prefix}/Centromere_CENH3/Trim_galore/{srarunid}/{srarunid}_R2.fastq.gz",
        genome = "results/{{prefix}}/Centromere_CENH3/BWA_Alignment/{}".format(os.path.basename(config["chr_genome"])),
        idx = rules.centromere_regions_3_bwa_index.output
    output:
        sam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.sam"),
        bam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.bam"),
        sorted_bam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.sorted.bam"),
        rmdup_bam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.bam"),
        marked_dup_metrics = "results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.marked_dup_metrics.txt",
        sorted_rmdup_bam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.sorted.bam"),
        sorted_rmdup_bam_bai = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.sorted.bam.bai"),
        q20_bam = temp("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.q20.bam"),
        q20_sorted_bam = "results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.q20.sorted.bam",
        flagstat = "results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.ChIP.flagstat"
    threads: 10
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        srarunid = "{srarunid}",
        index = lambda w, input: os.path.splitext(input.idx[0])[0]
    shell:
        '''
        module load SAMtools/1.9
        module load BWA/0.7.17
        bwa mem -t {threads} {input.genome} {input.fq1} {input.fq2} > {output.sam}
        samtools view -@ {threads} -b -o {output.bam} {output.sam}
        samtools sort -@ {threads} -o {output.sorted_bam} -T {params.srarunid}.ChIP {output.bam}

        module load picard/2.23.9
        java -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates I={output.sorted_bam} O={output.rmdup_bam} M={output.marked_dup_metrics} REMOVE_DUPLICATES=true
        samtools sort -@ {threads} -o {output.sorted_rmdup_bam} -T {params.srarunid}.ChIP {output.rmdup_bam}
        samtools index -@ {threads} {output.sorted_rmdup_bam}

        samtools view -@ {threads} -bhq 20 {output.sorted_rmdup_bam} -o {output.q20_bam}
        samtools sort -@ {threads} -o {output.q20_sorted_bam} -T {params.srarunid}.ChIP {output.q20_bam}
        samtools flagstat {output.q20_sorted_bam} > {output.flagstat}
        samtools index -@ {threads} {output.q20_sorted_bam}
        '''

rule centromere_regions_5_SICER2:
    input:
        expand("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.q20.sorted.bam", prefix = config["prefix"], srarunid = config["CENH3_Chip_srarunid"])
    output:
        "results/{}/Centromere_CENH3/SICER2/{}.rmdup.q20.sorted-W{}-G{}-islands-summary".format(config["prefix"], config["CENH3_Chip_srarunid"][0], config["SICER2"]["windows"], config["SICER2"]["gaps"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        species = config["SICER2"]["species"],
        windows = config["SICER2"]["windows"],
        gaps = config["SICER2"]["gaps"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        module load BEDTools/2.27
        ~/miniconda2/envs/SnakemakeV5.15/bin/sicer -t {input[0]} -c {input[1]} -s {params.species} -w {params.windows} -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g {params.gaps} -e 1000 -o {params.outdir} -cpu {threads}
        '''

rule centromere_regions_6_SICER2_filter:
    input:
        "results/{}/Centromere_CENH3/SICER2/{}.rmdup.q20.sorted-W{}-G{}-islands-summary".format(config["prefix"], config["CENH3_Chip_srarunid"][0], config["SICER2"]["windows"], config["SICER2"]["gaps"])
    output:
        "results/{}/Centromere_CENH3/SICER2/{}.rmdup.q20.sorted-W{}-G{}-islands-summary_filter.bed".format(config["prefix"], config["CENH3_Chip_srarunid"][0], config["SICER2"]["windows"], config["SICER2"]["gaps"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        fold_change_control = config["SICER2"]["fold_change_control"],
        false_discovery_rate = config["SICER2"]["false_discovery_rate"]
    log:
        "logs/{}/Centromere_Annot/SICER2.log".format(config["prefix"])
    shell:
        '''
        cat {input} | awk -F"\\t" 'BEGIN{{OFS="\\t"; print"chrom\\tstart\\tend\\tChIP_island_read_count\\tCONTROL_island_read_count\\tp_value\\tfold_change\\tFDR_threshold"}} $7 >= {params.fold_change_control} && $8 < {params.false_discovery_rate} {{print $0}}' > {output}
        '''


rule centromere_regions_7_epic2_effective_genome_fraction:
    input:
        genome = config["chr_genome"]
    output:
        jf = "results/{}/Centromere_CENH3/epic2/{}_genome_mer_counts.jf".format(config["prefix"], config["prefix"]),
        stats = "results/{}/Centromere_CENH3/epic2/{}_genome.stats".format(config["prefix"], config["prefix"])
    threads: 10
    conda:
        "../envs/Centromere_Annot.yaml"
    shell:
        '''
        module load Jellyfish/2.3.0
        jellyfish count -t {threads} -m 150 -s 200M \
          --out-counter-len 1 --counter-len 1 -o {output.jf} {input.genome}
        jellyfish stats {output.jf} > {output.stats}
        '''

rule centromere_regions_8_epic2:
    input:
        bams = expand("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.q20.sorted.bam", prefix = config["prefix"], srarunid = config["CENH3_Chip_srarunid"]),
        genome = config["chr_genome"],
        stats = "results/{}/Centromere_CENH3/epic2/{}_genome.stats".format(config["prefix"], config["prefix"])
    output:
        chr_size = "results/{}/Centromere_CENH3/epic2/{}_genome_chr.size".format(config["prefix"], config["prefix"]),
        epic2 = "results/{}/Centromere_CENH3/epic2/{}_epic2_peaks.bed".format(config["prefix"], config["prefix"]),
        filtered_epic2 = "results/{}/Centromere_CENH3/epic2/{}_epic2_peaks.filtered.bed".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Centromere_Annot.yaml"
    log:
        "logs/{}/Centromere_Annot/epic2.log".format(config["prefix"])
    shell:
        '''
        egf=`cat {input.stats} | \\grep -E "Unique|Total" | sed "s/ //g" | cut -f 2 -d \\: | xargs | awk '{{printf ("%.2f\\n",$1/$2)}}'`
        bioawk -c fastx '{{ print $name, length($seq) }}' < {input.genome} | sort -k1,1V > {output.chr_size}
        epic2 -t {input.bams[0]} -c {input.bams[1]} --chromsizes {output.chr_size} --mapq 20 --effective-genome-fraction ${{egf}} --bin-size 5000 --gaps-allowed 0 -o {output.epic2} > {log} 2>&1
        cat {output.epic2} | awk '{{if($10>2&&$5>250){{print$0}}}}' > {output.filtered_epic2}
        '''


rule centromere_regions_9_Enrichment_calculation:
    """
    https://blog.csdn.net/samhuairen/article/details/107204924
    """
    input:
        expand("results/{prefix}/Centromere_CENH3/BWA_Alignment/{srarunid}.rmdup.q20.sorted.bam", prefix = config["prefix"], srarunid = config["CENH3_Chip_srarunid"])
    output:
        "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph"
    params:
        window = "{window}"
    threads: 5
    shell:
        '''
        module load deepTools/3.5.0
        bamCompare -b1 {input[0]} -b2 {input[1]} -p {threads} -bs {params.window} -o {output} -of bedgraph --minMappingQuality 20 --ignoreDuplicates --normalizeUsing RPKM --scaleFactorsMethod None
        '''

rule centromere_regions_10_genome_size:
    input:
        genome = config["chr_genome"]
    output:
        genome_size = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}_genome_size.txt"
    shell:
        """
        bioawk -c fastx '{{ print $name, "1", length($seq) }}' < {input.genome} | sort -k1,1V | sed '1i chr\\tstart\\tend' > {output.genome_size}
        """

rule centromere_regions_11_Enrichment_calculation_plot:
    input:
        genome_size = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}_genome_size.txt",
        bedgraph = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph"
    output:
        pdf = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph.pdf",
        png = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph.png",
        pdf2 = "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph2.pdf"
    shell:
        '''
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        Rscript workflow/scripts/Centromere_Region_Plot.R {input.bedgraph} {output.pdf} {output.png}
        Rscript workflow/scripts/Centromere_Region_Plot2.R {input.genome_size} {input.bedgraph} {output.pdf2}
        '''

rule centromere_regions_12_Island_identification:
    input:
        "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.ChIP_Input.RPKM.w{window}.bedgraph"
    output:
        "results/{prefix}/Centromere_CENH3/DeepTools_bamCompare/{prefix}.centromere.w{window}.bg"
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"]
    shell:
        '''
        module load BEDTools/2.27
        Enriched_Islands_Ratio="0.2"
        cat {input} | awk -v ratio=${{Enriched_Islands_Ratio}} '{{if($4 > ratio) {{print$0}}}}' | bedtools merge -i - -d 1000000 -o mean -c 4 | grep -v -i "scaff" | \
          sort -k1,1V -k4,4nr | awk -F"\\t" '!a[$1]++' > {output}
        '''