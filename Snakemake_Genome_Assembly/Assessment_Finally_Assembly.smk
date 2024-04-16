def select_Assessment_input_genome():
    return "results/{}/Plus_finally_genome/{}_genome.fasta".format(config["prefix"], config["prefix"])

def select_Assessment_input_Chr_genome():
    return "results/{}/Plus_finally_genome/Chr_{}_genome.fasta".format(config["prefix"], config["prefix"])

def input_genome_name():
    input_file = select_Assessment_input_genome()
    name = os.path.basename(input_file)
    return name

localrules:
    AFA_26_genome_windows,
    AFA_34_bowtie2_stat

##################################################################################################################################
################################################# 1.  #######################################################
##################################################################################################################################
rule AFA_1_assembly_stat:
    input:
        select_Assessment_input_genome()
    output:
        stat = "results/{}/Assessment_Finally_Assembly/{}.stat".format(config["prefix"], input_genome_name()),
        stat2 = "results/{}/Assessment_Finally_Assembly/{}_seqkit.stat".format(config["prefix"], input_genome_name()),
        length = "results/{}/Assessment_Finally_Assembly/{}.length".format(config["prefix"], input_genome_name())
    threads: 10
    params:
        workdir = config["workdir"]
    shell:
        '''
        cd {params.workdir}/workflow/bin
        ./assemblathon_stats.pl {params.workdir}/{input} > {params.workdir}/{output.stat}
        seqkit stats -a {params.workdir}/{input} > {params.workdir}/{output.stat2}
        bioawk -c fastx '{{ print $name, length($seq) }}' < {params.workdir}/{input} | sort -k1,1V > {params.workdir}/{output.length}
        '''

rule AFA_3_BUSCO:
    input:
        genome = select_Assessment_input_genome()
    output:
        summary = report("results/{}/Assessment_Finally_Assembly/BUSCO/{}_{{BUSCODB}}/short_summary.specific.{{BUSCODB}}.{}.txt".format(config["prefix"], os.path.basename(select_Assessment_input_genome()), os.path.basename(select_Assessment_input_genome())), caption="../report/BUSCO.rst", category="9. Assessment_Finally_Assembly", subcategory="BUSCO"),
    threads: 20
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    params:
        workdir = config["workdir"],
        prefix = lambda w, input: os.path.basename(input.genome),
        BUSCODB = "{BUSCODB}",
        outdir = lambda w, output: "/".join(os.path.dirname(output['summary']).split("/")[0:4])
    shell:
        '''
        cd {params.outdir}
        busco -i {params.workdir}/{input.genome} --out_path {params.workdir}/{params.outdir} -o {params.prefix} --lineage_dataset {params.workdir}/resources/busco_downloads/lineages/{params.BUSCODB} -f -m genome --augustus -c {threads} --long --offline
        '''

rule AFA_4_BBMap_stats:
    input:
        genome = select_Assessment_input_genome()
    output:
        stat = "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats.txt".format(config["prefix"], config["prefix"]),
        gc = "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_ACGTN.txt".format(config["prefix"], config["prefix"]),
        gchist = "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_ACGTN_content_histogram.hist".format(config["prefix"], config["prefix"]),
        shist = "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_cumulative_scaffold_length_histogram.hist".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        stats.sh in={input.genome} out={output.stat} gc={output.gc} gchist={output.gchist} shist={output.shist} phs=t extended=t score=t
        sed -i "s/^#//" {output.gc}
        '''

rule AFA_4_BBMap_stats_plot:
    input:
        "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_ACGTN.txt".format(config["prefix"], config["prefix"])
    output:
        "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_ACGTN.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        Rscript workflow/scripts/BBMap_stats_ACGTN.R {input} {output}
        '''

rule AFA_5_seqkit_fx2tab:
    input:
        genome = select_Assessment_input_genome()
    output:
        "results/{}/Assessment_Finally_Assembly/{}_seqkit_fx2tab.txt".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        seqkit fx2tab {input.genome} -l -g -G -n -i -H -B AT -B GC -B N > {output}
        '''

rule AFA_6_Genome_GC_content:
    input:
        genome = select_Assessment_input_genome()
    output:
        size = "results/{}/Assessment_Finally_Assembly/GC_content/{}.size".format(config["prefix"], config["prefix"]),
        windows = "results/{}/Assessment_Finally_Assembly/GC_content/{}_1kb.bed".format(config["prefix"], config["prefix"]),
        windows_nuc = "results/{}/Assessment_Finally_Assembly/GC_content/{}_nuc_1kb.txt".format(config["prefix"], config["prefix"]),
        windows_igv = "results/{}/Assessment_Finally_Assembly/GC_content/{}_nuc_1kb.igv".format(config["prefix"], config["prefix"]),
        windows_tdf = "results/{}/Assessment_Finally_Assembly/GC_content/{}_nuc_1kb.tdf".format(config["prefix"], config["prefix"])
    threads: 5
    params:
        workdir = config["workdir"],
        prefix = config["prefix"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load IGVTools/2.3.98-Java-1.8.0_45
        width=1000
        bioawk -c fastx '{{ print $name, length($seq) }}' < {input.genome} | sort -k1,1V > {output.size}
        bedtools makewindows -g {output.size} -w ${{width}} > {output.windows}
        bedtools nuc -fi {input.genome} -bed {output.windows} > {output.windows_nuc}
        gawk -v w=${{width}} 'BEGIN{{FS="\\t"; OFS="\\t"}}
            {{
            if (FNR>1) {{print $1,$2,$3,"GCpc_"w"bps",$5}}
            }}' {output.windows_nuc} > {output.windows_igv}
        igvtools toTDF -z 5 -f min,max,mean {output.windows_igv} {output.windows_tdf} {input.genome}
        '''

rule AFA_7_Gap_Stats:
    input:
        genome = select_Assessment_input_genome()
    output:
        "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_Gaps.txt".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        python3 {params.workdir}/workflow/scripts/Genome_Gap_Stats.py {input} {output}
        '''

rule AFA_8_telomeric:
    input:
        genome = select_Assessment_input_genome()
    output:
        "results/{}/Assessment_Finally_Assembly/Telomeric/{}_genome_telomeric.txt".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        python {params.workdir}/workflow/scripts/FindTelomeres.py {input} --window 1000 -c 10 > {output}
        '''

rule AFA_11_Gap_Stats_telomeric_plot:
    input:
        genome = select_Assessment_input_genome(),
        telomeric = "results/{}/Assessment_Finally_Assembly/Telomeric/{}_genome_telomeric.txt".format(config["prefix"], config["prefix"]),
        gaps = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_Gaps.txt".format(config["prefix"], config["prefix"])
    output:
        genome_size = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_size.txt".format(config["prefix"], config["prefix"]),
        telomeric = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_telomeric.txt".format(config["prefix"], config["prefix"]),
        plot = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_Gaps.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        bioawk -c fastx '{{ print $name, "1", length($seq) }}' < {input.genome} | sort -k1,1V | \\grep -i -v -E "scaffold|contig" | sed '1i chr\\tstart\\tend' > {output.genome_size}
        cat {input.telomeric} | grep -v -E "^#|Telomeres found|sequences to analyze|^$|forward:|reverse:" | \\grep -i "chr" | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$4,$5,$6,$7}}' > {output.telomeric}
        cat {output.telomeric} | cut -f 1 | sort -V | uniq | grep -v -w -f - {output.genome_size} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"middle","gpos25"}}' >> {output.telomeric}
        Rscript {params.workdir}/workflow/scripts/Genome_Gap_Stats_Plot.R {output.genome_size} {output.telomeric} {input.gaps} {output.plot}
        '''
#################################################################################################################################
########################################################## 2. T2T ###########################################################
#################################################################################################################################
rule AFA_12_meryl_count_genome:
    input:
        genome = select_Assessment_input_genome()
    output:
        merylDB = temp(directory("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_{}kmer.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))),
        repetitive = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))
    threads: 5
    params:
        kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        meryl count k={params.kmer} output {output.merylDB} {input.genome}
        meryl print greater-than distinct=0.9998 {output.merylDB} > {output.repetitive}
        '''

rule AFA_13_winnowmap_mapping_HiFi:
    input:
        repetitive = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        ghifi = "results/{}/Pseudomolecule_Polishing_Gap_Filling/Additional_Polishing_HiFi/{}_g16k_hifi_reads.fasta".format(config["prefix"], config["prefix"])
    output:
        sam = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome.sam".format(config["prefix"], config["prefix"]))
    params:
        kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    threads: 5
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        winnowmap -k {params.kmer} -W {input.repetitive} --MD -ax map-pb {input.genome} {input.ghifi} > {output.sam}
        '''

rule AFA_14_winnowmap_mapping_ONT:
    input:
        repetitive = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        ont = "results/{}/Pseudomolecule_Polishing_Gap_Filling/ONT/{}_R9ONT.fasta.gz".format(config["prefix"], config["prefix"])
    output:
        sam = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_ont2genome.sam".format(config["prefix"], config["prefix"]))
    params:
        kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    threads: 5
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        winnowmap -k {params.kmer} -W {input.repetitive} --MD -ax map-pb {input.genome} {input.ont} > {output.sam}
        '''

rule AFA_15_bwa_index:
    input:
        genome = select_Assessment_input_genome()
    output:
        multiext("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}".format(config["prefix"], input_genome_name()), ".amb", ".ann", ".pac"),
        "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}".format(config["prefix"], input_genome_name())
    params:
        prefix = lambda w, output: ".".join(output[0].split(".")[:-1]),
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load BWA/0.7.17
        bwa index -p {params.prefix} {input.genome}
        cp {input.genome} {params.outdir}
        '''

rule AFA_16_bwa:
    input:
        unpack(get_illumina_pe_fastq),
        genome = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}".format(config["prefix"], input_genome_name()),
        idx = rules.AFA_15_bwa_index.output
    output:
        filter_sorted_bam = temp("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{sample_illumina_pe}-{unit_illumina_pe}.sorted.bam")
    threads: 5
    conda:
        "../envs/HiC_QC.yaml"
    shell:
        '''
        module load BWA/0.7.17
        bwa mem -t {threads} {input.genome} {input.r1} {input.r2} | samblaster | samtools view -@ {threads} -S -b > {output.filter_sorted_bam}
        '''

rule AFA_17_bwa_merge_bam:
    input:
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{u.sample}-{u.unit}.sorted.bam", prefix= config["prefix"], u=units[units['platform'].str.match('ILLUMINA_PE')].itertuples())
    output:
        sam = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_Illumina2genome.sam".format(config["prefix"], config["prefix"]),
        header = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_sam_header.txt".format(config["prefix"], config["prefix"])
    threads: 5
    shell:
        '''
        samtools merge -@ {threads} --output-fmt SAM {output.sam} {input}
        cat {output.sam} | \\grep "^@SQ" > {output.header}
        '''

rule AFA_18_winnowmap_mapping_filter_hifi:
    input:
        header = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_sam_header.txt".format(config["prefix"], config["prefix"]),
        hifi_sam = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome.sam".format(config["prefix"], config["prefix"]),
    output:
        hifi_sam = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome_F256_Q20.sam".format(config["prefix"], config["prefix"]))
    threads: 20
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        samtools view -@ {threads} -F 256 -q 20 -S {input.hifi_sam} | cat {input.header} - > {output.hifi_sam}
        '''

rule AFA_20_winnowmap_mapping2Paf_illumina:
    input:
        illumina_sam = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_Illumina2genome.sam".format(config["prefix"], config["prefix"]),
    output:
        illumina_paf = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_Illumina2genome.paf".format(config["prefix"], config["prefix"])),
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        paftools.js sam2paf {input.illumina_sam} > {output.illumina_paf}
        '''

rule AFA_21_winnowmap_mapping2Paf_hifi:
    input:
        hifi_sam = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome_F256_Q20.sam".format(config["prefix"], config["prefix"])
    output:
        hifi_paf = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome_F256_Q20.paf".format(config["prefix"], config["prefix"]))
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        paftools.js sam2paf {input.hifi_sam} > {output.hifi_paf}
        '''

rule AFA_22_winnowmap_mapping2Paf_ONT:
    input:
        ont_sam = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_ont2genome_F256.sam".format(config["prefix"], config["prefix"])
    output:
        ont_paf = temp("results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_ont2genome_F256.paf".format(config["prefix"], config["prefix"]))
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        paftools.js sam2paf {input.ont_sam} > {output.ont_paf}
        '''

rule AFA_26_genome_size:
    input:
        genome = select_Assessment_input_genome()
    output:
        genome_size = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_genome_size.txt"
    shell:
        """
        bioawk -c fastx '{{ print $name, "1", length($seq) }}' < {input.genome} | sort -k1,1V | sed '1i chr\\tstart\\tend' > {output.genome_size}
        """

rule AFA_27_genome_windows:
    input:
        genome = select_Assessment_input_genome()
    output:
        windows_bed = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_genome_w{windows}.bed"
    params:
        windows = "{windows}"
    shell:
        """
        module load BEDTools/2.27
        samtools faidx {input.genome}
        bedtools makewindows -g {input.genome}.fai -w {params.windows} > {output.windows_bed}
        """

rule AFA_28_windows_statistics:
    input:
        stats_bg = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome.{asstype}.bedgraph",
        windows_bed = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_genome_w{windows}.bed"
    output:
        windows_bedgraph = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome.{asstype}.w{windows}.{statistics}.bedgraph",
    conda:
        "../envs/Genome_Mappability.yaml"
    params:
        statistics = "{statistics}"
    shell:
        """
        module load BEDTools/2.27
        bedtools map -a {input.windows_bed} -b {input.stats_bg} -c 4 -o {params.statistics} | grep -v -P "\\t\\.$" > {output.windows_bedgraph}
        """

rule AFA_29_windows_statistics_all_plot:
    """
    We used these alignments to generate summary statistics (min, max, mean, median) for 10 kbp non-overlapping windows of the \
      assembly for each of: alignment length, read length, identity, and MAPQ.
    """
    input:
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome.{asstype}.w{windows}.{statistics}.bedgraph", prefix = config["prefix"], assdata = config["params"]["assessment_t2t_data"], windows = config["params"]["assessment_t2t_windows"], asstype = config["params"]["assessment_t2t_type"], statistics = config["params"]["assessment_t2t_tatistics"])
    output:
        stat = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_T2T_Alignment.stat".format(config["prefix"], config["prefix"]),
        plot = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_T2T_Alignment.pdf".format(config["prefix"], config["prefix"])
    threads: 10
    shell:
        """
        module load R/3.6.0
        module load GCC/7.2.0-2.29

        for i in {input}
        do
        name=`basename ${{i}}`
        assdata=`echo ${{name}} | cut -f 1 -d \\. | cut -f 2 -d \\_ | sed "s/2genome//"`
        asstype=`echo ${{name}} | cut -f 2 -d \\.`
        windows=`echo ${{name}} | cut -f 3 -d \\.`
        statistics=`echo ${{name}} | cut -f 4 -d \\.`
        cat ${{i}} | awk -F"\\t" -v name="${{name}}" -v assdata="${{assdata}}" -v asstype="${{asstype}}" -v windows="${{windows}}" -v statistics="${{statistics}}" 'BEGIN{{OFS="\\t"}} {{print $0,name,assdata,asstype,windows,statistics}}' >> {output.stat}
        done

        sed -i '1i chr\\tstart\\tend\\tcvalue\\tfilename\\tdatatype\\tstatistics\\tbinsize\\tstatisticstype' {output.stat}
        Rscript workflow/scripts/t2t_mapping_statistics_plot.R {output.stat} {output.plot}
        """

rule AFA_30_winnowmap_mapping_sam2bam:
    input:
        sam = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.sam"
    output:
        bam = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.bam"
    threads: 10
    conda:
        "../envs/Genome_Mappability.yaml"
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} | samtools sort -@ {threads} - > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule AFA_31_winnowmap_mapping_coverage_plot_karyoploteR:
    input:
        bam = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.bam",
        genome_size = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_genome_size.txt"
    output:
        bedgraph = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.w{windows}.coverage.bedgraph",
        bedgraph_plot_line = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.w{windows}.coverage.bedgraph.line.pdf",
        bedgraph_plot_bar = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.w{windows}.coverage.bedgraph.bar.pdf"
    threads: 10
    conda:
        "../envs/Genome_Mappability.yaml"
    params:
        windows = "{windows}"
    shell:
        """
        module load deepTools/3.5.0
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        bamCoverage --bam {input.bam} -o {output.bedgraph} \
            --outFileFormat bedgraph \
            --numberOfProcessors {threads} \
            --binSize {params.windows} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize `cat {input.genome_size} | awk -F'\\t' '{{sum+=$3;}} END{{print sum;}}'`

        Rscript workflow/scripts/genmap_mappability_plot.R {input.genome_size} {output.bedgraph} {output.bedgraph_plot_line} {output.bedgraph_plot_bar}
        """

rule AFA_32_winnowmap_mapping_coverage_plot_tinycov:
    input:
        bam = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.bam",
        windows_bed = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_genome_w{windows}.bed"
    output:
        tinycov_covplot = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot.pdf",
        tinycov_covplot_bg = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot.bedgraph",
        tinycov_covhist = "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covhist.pdf"
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    shell:
        '''
        tinycov covplot --max-depth 50 --out {output.tinycov_covplot} --whitelist {params.chromosome} --bins {input.windows_bed} --text {output.tinycov_covplot_bg} --ploidy 0 {input.bam}
        tinycov covhist --max-depth 50 --out {output.tinycov_covhist} --whitelist {params.chromosome} --bins {input.windows_bed} {input.bam}
        '''

rule AFA_32_winnowmap_mapping_coverage_plot_tinycov_Rplot:
    input:
        "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot.bedgraph"
    output:
        "results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot_Rplot.pdf"
    threads: 1
    shell:
        '''
        # R plot tinycov covplot output file
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        Rscript workflow/scripts/t2t_mapping_tinycov_plot.R {input} 50 {output}
        '''

rule AFA_32_winnowmap_mapping_coverage_plot_list:
    input:
        "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_T2T_Alignment.pdf".format(config["prefix"], config["prefix"]),
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.w{windows}.coverage.bedgraph.line.pdf", prefix = config["prefix"], assdata = config["params"]["assessment_t2t_data"], windows = config["params"]["assessment_t2t_windows"]),
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot.pdf", prefix = config["prefix"], assdata = config["params"]["assessment_t2t_data"], windows = config["params"]["assessment_t2t_windows"]),
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}_w{windows}_tinycov_covplot_Rplot.pdf", prefix = config["prefix"], assdata = config["params"]["assessment_t2t_data"], windows = config["params"]["assessment_t2t_windows"])

##################################################################################################################################
########################################################### 3.  ###########################################################
##################################################################################################################################

rule AFA_32_genome_dict_fai_bowtie2_index:
    input:
        genome = select_Assessment_input_genome()
    output:
        gfai = "{}.fai".format(rules.AFA_3_BUSCO.input.genome),
        gdict = "{}.dict".format(os.path.splitext(rules.AFA_3_BUSCO.input.genome)[0])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        bowtie2-build {input.genome} {input.genome}
        samtools faidx {input.genome}
        picard CreateSequenceDictionary R={input.genome} O={output.gdict}
        '''

rule AFA_33_bowtie2:
    """ Mapping Illumina reads to finally assembly genome """
    input:
        unpack(get_illumina_pe_fastq),
        gfai = "{}.fai".format(rules.AFA_3_BUSCO.input.genome),
        gdict = "{}.dict".format(os.path.splitext(rules.AFA_3_BUSCO.input.genome)[0]),
        genome = select_Assessment_input_genome()
    output:
        #sam = temp("results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}.sam"),
        sam = "results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}.sam",
        log = "results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}_bowtie2.log"
    threads: 10
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        bowtie2 --threads {threads} -x {input.genome} -1 {input.r1} -2 {input.r2} -S {output.sam} > {output.log} 2>&1
        sed -i '1i\\{input.r1} and {input.r2}\\n' {output.log}
        '''

rule AFA_34_bowtie2_stat:
    input:
        expand("results/{prefix}/Assessment_Finally_Assembly/bowtie2/{u.sample}-{u.unit}_bowtie2.log", prefix= config["prefix"], u=units[units['platform'].str.match('ILLUMINA_PE')].itertuples())
    output:
        report("results/{}/Assessment_Finally_Assembly/bowtie2/All_sample_bowtie2.log".format(config["prefix"]), caption="../report/bowtie2.rst", category="9. Assessment_Finally_Assembly", subcategory="bowtie2 mapping Illumina")
    threads: 1
    shell:
        '''
        cat {input} > {output}
        '''

rule AFA_35_bowtie2_sorted_bam:
    """ Mapping Illumina reads to finally assembly genome """
    input:
        "results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}.sam"
    output:
        bam = "results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}.bam",
        sorted_bam = "results/{prefix}/Assessment_Finally_Assembly/bowtie2/{sample_illumina_pe}-{unit_illumina_pe}.sorted.bam"
    threads: 10
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        samtools view -bS -@ {threads} {input} -o {output.bam}
        samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
        samtools index -@ {threads} {output.sorted_bam}
        '''

rule AFA_36_bowtie2_merge_bam:
    input:
        expand("results/{prefix}/Assessment_Finally_Assembly/bowtie2/{u.sample}-{u.unit}.sorted.bam", prefix= config["prefix"], u=units[units['platform'].str.match('ILLUMINA_PE')].itertuples())
    output:
        "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_bowtie2.sorted.bam".format(config["prefix"], config["prefix"])
    threads: 10
    shell:
        '''
        samtools merge -@ {threads} {output} {input}
        samtools index -@ {threads} {output}
        '''

rule AFA_37_bowtie2_coverage_plot:
    input:
        genome = select_Assessment_input_genome(),
        gfai = "{}.fai".format(rules.AFA_3_BUSCO.input.genome),
        gdict = "{}.dict".format(os.path.splitext(rules.AFA_3_BUSCO.input.genome)[0]),
        sort_bam = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_bowtie2.sorted.bam".format(config["prefix"], config["prefix"])
    output:
        physical_coverage = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_physical_coverage.xls".format(config["prefix"], config["prefix"]),
        svg = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_wgscoverageplotter.svg".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"]),
        genome_prefix = lambda w, input: os.path.splitext(input.genome)[0]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## bamtocov
        bamtocov --physical --report {output.physical_coverage} -T {threads} {input.sort_bam}
        ## picard + wgscoverageplotter
        samtools faidx {input.genome}
        java -jar ${{EBROOTPICARD}}/picard.jar CreateSequenceDictionary R={input.genome} O={params.genome_prefix}.dict
        wgscoverageplotter.py --dimension 1500x500 -C -1 --clip --include-contig-regex "Chr.*" -R {input.genome} {input.sort_bam} --percentile median > {output.svg}
        '''

rule AFA_38_bowtie2_coverage_plot_tinycov:
    input:
        sort_bam = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_bowtie2.sorted.bam".format(config["prefix"], config["prefix"])
    output:
        tinycov_covplot = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_tinycov_covplot.pdf".format(config["prefix"], config["prefix"]),
        tinycov_covhist = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_tinycov_covhist.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    shell:
        '''
        ## tinycov
        tinycov covplot --out {output.tinycov_covplot} --whitelist {params.chromosome} --ploidy 0 {input.sort_bam}
        tinycov covhist --out {output.tinycov_covhist} --whitelist {params.chromosome} {input.sort_bam}
        '''

rule AFA_39_PB_coverage_map:
    input:
        genome = select_Assessment_input_genome(),
        subreads_fasta = "results/{}/Pseudomolecule_Polishing_Gap_Filling/Additional_Polishing_HiFi/{}_g16k_hifi_reads.fasta".format(config["prefix"], config["prefix"])
    output:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"]),
        depth = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb.depth".format(config["prefix"], config["prefix"])
    threads: 10
    params:
        bigGenome = config["bigGenome"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        if [ {params.bigGenome} = True ] ; then
            echo "For Big Genome"
            minimap2 -ax map-hifi -t {threads} {input.genome} {input.subreads_fasta} | samtools view -@ {threads} -T {input.genome} -bS - | samtools sort -@ {threads} -o {output.sort_bam} -
        else
            echo "For Normal Genome"
            minimap2 -ax map-hifi -t {threads} {input.genome} {input.subreads_fasta} | samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output.sort_bam} -
        fi
        samtools depth {output.sort_bam} > {output.depth}
        '''

rule AFA_40_PB_coverage:
    input:
        depth = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb.depth".format(config["prefix"], config["prefix"])
    output:
        cov = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_cov.txt".format(config["prefix"], config["prefix"])
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        python3 {params.workdir}/workflow/scripts/minBinDepth.py {input.depth} {output} 5000000
        '''

rule AFA_41_PB_coverage_plot:
    input:
        genome = select_Assessment_input_genome(),
        gfai = "{}.fai".format(rules.AFA_3_BUSCO.input.genome),
        gdict = "{}.dict".format(os.path.splitext(rules.AFA_3_BUSCO.input.genome)[0]),
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"])
    output:
        physical_coverage = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_physical_coverage.xls".format(config["prefix"], config["prefix"]),
        svg = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_wgscoverageplotter.svg".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## bamtocov
        bamtocov --physical --report {output.physical_coverage} -T {threads} {input.sort_bam}
        ## picard + wgscoverageplotter
        wgscoverageplotter.py --dimension 1500x500 -C -1 --clip --include-contig-regex "Ghjin_[AD].*" -R {input.genome} {input.sort_bam} --percentile median > {output.svg}
        '''

rule AFA_42_PB_coverage_plot_tinycov:
    input:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"])
    output:
        tinycov_covplot = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_tinycov_covplot.png".format(config["prefix"], config["prefix"]),
        tinycov_covhist = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_tinycov_covhist.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    shell:
        '''
        ## tinycov
        tinycov covplot --out {output.tinycov_covplot} --whitelist {params.chromosome} --ploidy 0 {input.sort_bam}
        tinycov covhist --out {output.tinycov_covhist} --whitelist {params.chromosome} {input.sort_bam}
        '''

rule AFA_43_PB_coverage_plot_qualimap:
    input:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"]),
        genome = select_Assessment_input_genome()
    output:
        directory("results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_qualimap".format(config["prefix"], config["prefix"]))
    threads: 10
    params:
        prefix = config["prefix"]
    conda: "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        unset DISPLAY
        qualimap bamqc -bam {input.sort_bam} -nt {threads} -outdir {output} -outfile {params.prefix} -c --java-mem-size=100G
        '''

rule AFA_44_ONT_coverage_map:
    input:
        genome = select_Assessment_input_genome(),
        subreads_fasta = "results/{}/Pseudomolecule_Polishing_Gap_Filling/ONT/{}_R9ONT.fasta.gz".format(config["prefix"], config["prefix"])
    output:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"]),
        depth = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT.depth".format(config["prefix"], config["prefix"])
    threads: 10
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        minimap2 -ax map-ont -t {threads} {input.genome} {input.subreads_fasta} | samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output.sort_bam} -
        samtools depth {output.sort_bam} > {output.depth}
        '''

rule AFA_45_ONT_coverage:
    input:
        depth = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT.depth".format(config["prefix"], config["prefix"])
    output:
        cov = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_cov.txt".format(config["prefix"], config["prefix"])
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        python3 {params.workdir}/workflow/scripts/minBinDepth.py {input.depth} {output} 5000000
        '''

rule AFA_46_ONT_coverage_plot:
    input:
        genome = select_Assessment_input_genome(),
        gfai = "{}.fai".format(rules.AFA_3_BUSCO.input.genome),
        gdict = "{}.dict".format(os.path.splitext(rules.AFA_3_BUSCO.input.genome)[0]),
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"])
    output:
        physical_coverage = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_physical_coverage.xls".format(config["prefix"], config["prefix"]),
        svg = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_wgscoverageplotter.svg".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## bamtocov
        bamtocov --physical --report {output.physical_coverage} -T {threads} {input.sort_bam}
        ## picard + wgscoverageplotter
        wgscoverageplotter.py --dimension 1500x500 -C -1 --clip --include-contig-regex "Gh.*" -R {input.genome} {input.sort_bam} --percentile median > {output.svg}
        '''

rule AFA_47_ONT_coverage_plot_tinycov:
    input:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"])
    output:
        tinycov_covplot = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_tinycov_covplot.pdf".format(config["prefix"], config["prefix"]),
        tinycov_covhist = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_tinycov_covhist.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    shell:
        '''
        ## tinycov
        tinycov covplot --out {output.tinycov_covplot} --whitelist {params.chromosome} --ploidy 0 {input.sort_bam}
        tinycov covhist --out {output.tinycov_covhist} --whitelist {params.chromosome} {input.sort_bam}
        '''

rule AFA_48_ONT_coverage_plot_qualimap:
    input:
        sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"]),
        genome = select_Assessment_input_genome()
    output:
        directory("results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_qualimap".format(config["prefix"], config["prefix"]))
    threads: 10
    params:
        prefix = config["prefix"]
    conda: "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        unset DISPLAY
        qualimap bamqc -bam {input.sort_bam} -nt {threads} -outdir {output} -outfile {params.prefix} -c --java-mem-size=1000G
        '''

rule AFA_49_Illumina_PB_ONT_coverage_plot:
    input:
        genome_size = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_size.txt".format(config["prefix"], config["prefix"]),
        centromere = "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_telomeric.txt".format(config["prefix"], config["prefix"]),
        #centromere = config["params"]["centromere"]["centromere_bed"],
        in_sort_bam = "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_bowtie2.sorted.bam".format(config["prefix"], config["prefix"]),
        pb_sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"]),
        ont_sort_bam = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"])
    output:
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_Illumina_PB_ONT_covplot.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        chromosome = ",".join(config["params"]["chromosome"])
    shell:
        '''
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        Rscript workflow/scripts/bam_coverage_plot.R {input.genome_size} {input.centromere} {input.in_sort_bam} {input.pb_sort_bam} {input.ont_sort_bam} {output}
        '''
##################################################################################################################################
######################################################## 5. merqury ##############################################################
##################################################################################################################################
rule AFA_51_merqury_QV_Illumina_Prepare_meryl_dbs_R1:
    input:
        unpack(get_illumina_pe_fastq)
    output:
        r1_meryl = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.R1.meryl".format(config["prefix"], config["prefix"])))
    threads: 1
    params:
        Assessment_Finally_Assembly_QV_kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## Build k-mer dbs with meryl
        meryl k={params.Assessment_Finally_Assembly_QV_kmer} count output {output.r1_meryl} {input.r1}
        '''

rule AFA_52_merqury_QV_Illumina_Prepare_meryl_dbs_R2:
    input:
        unpack(get_illumina_pe_fastq)
    output:
        r2_meryl = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.R2.meryl".format(config["prefix"], config["prefix"])))
    threads: 1
    params:
        Assessment_Finally_Assembly_QV_kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## Build k-mer dbs with meryl
        meryl k={params.Assessment_Finally_Assembly_QV_kmer} count output {output.r2_meryl} {input.r2}
        '''

rule AFA_53_2_merqury_QV_Illumina_Prepare_meryl_dbs_Merge:
    input:
        r_meryl = expand("results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{{prefix}}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.R{group}.meryl", group = [1,2])
    output:
        reads_meryl = directory("results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{{prefix}}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.k{}.meryl".format(config['params']['Assessment_Finally_Assembly_QV_kmer']))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        # 2. Merge
        meryl union-sum output {output.reads_meryl} {input.r_meryl}
        '''

rule AFA_54_2_merqury_QV_Illumina:
    input:
        reads_meryl = "results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{{prefix}}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.k{}.meryl".format(config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome()
    output:
        temp(directory("results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.QV/{}.meryl".format(config["prefix"], input_genome_name().split(".")[0]))),
        "results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.QV/{}.qv".format(config["prefix"], config["prefix"]),
        "results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.QV/{}.{}.qv".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]),
        "results/{{prefix}}/Assessment_Finally_Assembly/Merqury_QV/{}.{{sample_illumina_pe}}-{{unit_illumina_pe}}.Illumina.QV/{}.completeness.stats".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        cd {params.outdir}
        merqury.sh {params.workdir}/{input.reads_meryl} {params.workdir}/{input.genome} {params.prefix}
        '''

rule AFA_54_3_merqury_QV_Illumina_list:
    input:
        expand("results/{prefix}/Assessment_Finally_Assembly/Merqury_QV/{prefix}.{u.sample}-{u.unit}.Illumina.QV/{prefix}.qv", prefix= config["prefix"], u=units[units['platform'].str.match('ILLUMINA_PE')].itertuples(), group = [1,2])
    output:
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/Merqury_QV_Illumina_{}.txt".format(config["prefix"], config["prefix"])
    threads: 1
    shell:
        '''
        cat {input} > {output}
        '''

rule AFA_53_merqury_QV_Illumina_Prepare_meryl_dbs_Merge:
    input:
        r_meryl = expand("results/{prefix}/Assessment_Finally_Assembly/Merqury_QV/{prefix}.{u.sample}-{u.unit}.R{group}.meryl", prefix= config["prefix"], u=units[units['platform'].str.match('ILLUMINA_PE')].itertuples(), group = [1,2])
    output:
        reads_meryl = directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        # 2. Merge
        meryl union-sum output {output.reads_meryl} {input.r_meryl}
        '''

rule AFA_54_merqury_QV_Illumina:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome()
    output:
        temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.completeness.stats".format(config["prefix"], config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        cd {params.outdir}
        merqury.sh {params.workdir}/{input.reads_meryl} {params.workdir}/{input.genome} {params.prefix}
        '''

# HiFi QV
rule AFA_55_merqury_QV_HiFi_Prepare_meryl_dbs:
    input:
        hifi_fq = "results/{}/RawData/HiFi_reads/{}_All_hifi.fastq".format(config["prefix"], config["prefix"])
    output:
        hifi_meryl = directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))
    threads: 10
    params:
        Assessment_Finally_Assembly_QV_kmer = config['params']['Assessment_Finally_Assembly_QV_kmer']
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        ## Build k-mer dbs with meryl
        meryl k={params.Assessment_Finally_Assembly_QV_kmer} count output {output.hifi_meryl} {input.hifi_fq}
        '''

rule AFA_56_merqury_QV_HiFi:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome()
    output:
        temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.completeness.stats".format(config["prefix"], config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        cd {params.outdir}
        merqury.sh {params.workdir}/{input.reads_meryl} {params.workdir}/{input.genome} {params.prefix}
        '''

# Illumina + FiHi QV
rule AFA_57_merqury_QV_hybrid_gt1:
    input:
        reads_meryl_hifi = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        reads_meryl_Illumina = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer'])
    output:
        reads_meryl_hifi_gt1 = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.gt1.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))),
        reads_meryl_Illumina_gt1 = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.gt1.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer'])))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        meryl greater-than 1 {input.reads_meryl_hifi} output {output.reads_meryl_hifi_gt1}
        meryl greater-than 1 {input.reads_meryl_Illumina} output {output.reads_meryl_Illumina_gt1}
        '''

rule AFA_58_merqury_QV_hybrid_matching:
    input:
        reads_meryl_hifi_gt1 = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.gt1.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        reads_meryl_Illumina_gt1 = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.gt1.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer'])
    output:
        reads_meryl_hifi_gt1_add4 = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.gt1.add4.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))),
        reads_meryl_Illumina_gt1_div3 = temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.gt1.div3.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer'])))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        meryl divide-round 3 {input.reads_meryl_Illumina_gt1} output {output.reads_meryl_Illumina_gt1_div3}
        meryl increase 4 {input.reads_meryl_hifi_gt1} output {output.reads_meryl_hifi_gt1_add4}
        '''

rule AFA_59_merqury_QV_hybrid_union:
    input:
        reads_meryl_hifi_gt1_add4 = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.gt1.add4.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        reads_meryl_Illumina_gt1_div3 = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.gt1.div3.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer'])
    output:
        reads_meryl_hybrid = directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        meryl union-max {input.reads_meryl_Illumina_gt1_div3} {input.reads_meryl_hifi_gt1_add4} output {output}
        '''

rule AFA_60_merqury_QV_hybrid:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome()
    output:
        temp(directory("results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]),
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.completeness.stats".format(config["prefix"], config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        cd {params.outdir}
        merqury.sh {params.workdir}/{input.reads_meryl} {params.workdir}/{input.genome} {params.prefix}
        '''

rule AFA_61_merqury_QV_plot:
    input:
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.{}.qv".format(config["prefix"], config["prefix"], config["prefix"], input_genome_name().split(".")[0]) if config["module"]["Assessment_Finally_Assembly_Illumina"] and config["module"]["Assessment_Finally_Assembly_HiFi"] else []
    output:
        stat = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}_Merqury.QV.stat.txt".format(config["prefix"], config["prefix"], config["prefix"]),
        pdf = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}_Merqury.QV.pdf".format(config["prefix"], config["prefix"], config["prefix"])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        for i in {input}
        do
        type=`dirname ${{i}} | awk -F '/' '{{print $NF}}'`
        cat ${{i}} | awk -F"\\t" -v type=${{type}} 'BEGIN{{OFS="\\t"}} $1 !~ /Scaffold/ || $1 !~ /utg/ {{print $1,$4, type}}' >> {output.stat}
        done
        Rscript workflow/scripts/merqury_QV_plot.R {output.stat} {output.pdf}
        '''

##################################################################################################################################
######################################################### 6. Merfin ##############################################################
##################################################################################################################################

rule AFA_62_peak_estimate:
    input:
        Illumina_meryl_gt1 = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.gt1.meryl".format(config["prefix"], config["prefix"], config["params"]["Assessment_Finally_Assembly_QV_kmer"]),
        genome = select_Assessment_input_genome()
    output:
        reads_hist = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_{}kmer.reads.hist".format(config["prefix"], config["prefix"], config["params"]["Assessment_Finally_Assembly_QV_kmer"]),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        Assessment_Finally_Assembly_QV_kmer = config["params"]["Assessment_Finally_Assembly_QV_kmer"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        ploidy = config["Survey"]["ploidy"]
    conda:
        "../envs/Scaffolds_3_Additional_Polishing_HiFi.yaml"
    shell:
        '''
        meryl histogram {input.Illumina_meryl_gt1} > {output.reads_hist}
        genomescope2 -i {output.reads_hist} -k {params.Assessment_Finally_Assembly_QV_kmer} -o {params.outdir} -p {params.ploidy} --fitted_hist > {output.peak_value}
        '''

rule AFA_63_merfin_QV_Illumina:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        hist = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV.hist".format(config["prefix"], config["prefix"]),
        log = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV.hist.log".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | cut -f 2 -d :`
        workflow/bin/merfin/build/bin/merfin -hist -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} -output {output.hist} 2> {output.log}
        '''

rule AFA_64_merfin_kmer_completeness_Illumina:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV.completeness".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | cut -f 2 -d :`
        workflow/bin/merfin/build/bin/merfin -completeness -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} 2> {output}
        cat {output} | grep "COMPLETENESS:"
        '''

rule AFA_65_merfin_QV_HiFi:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        hist = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV.hist".format(config["prefix"], config["prefix"]),
        log = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV.hist.log".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | cut -f 2 -d :`
        workflow/bin/merfin/build/bin/merfin -hist -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} -output {output.hist} 2> {output.log}
        '''

rule AFA_66_merfin_kmer_completeness_HiFi:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV.completeness".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | awk -F":" '{{print $2*4}}'`
        workflow/bin/merfin/build/bin/merfin -completeness -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} 2> {output}
        cat {output} | grep "COMPLETENESS:"
        '''

rule AFA_67_merfin_QV_hybrid:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        hist = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV.hist".format(config["prefix"], config["prefix"]),
        log = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV.hist.log".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | cut -f 2 -d :`
        workflow/bin/merfin/build/bin/merfin -hist -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} -output {output.hist} 2> {output.log}
        '''

rule AFA_68_merfin_kmer_completeness_hybrid:
    input:
        reads_meryl = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.k{}.meryl".format(config["prefix"], config["prefix"], config['params']['Assessment_Finally_Assembly_QV_kmer']),
        genome = select_Assessment_input_genome(),
        peak_value = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/{}_peak.txt".format(config["prefix"], config["prefix"]),
        #lookup_table = "results/{}/Assessment_Finally_Assembly/Merfin_QV/GenomeScope/lookup_table.txt".format(config["prefix"], config["prefix"])
    output:
        #temp(directory("results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV/{}.meryl".format(config["prefix"], config["prefix"], input_genome_name().split(".")[0]))),
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV.completeness".format(config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load GCC/7.2.0-2.29
        peak=`cat {input.peak_value} | sed "s/ /\\n/g" | \\grep "kcov" | cut -f 2 -d :`
        workflow/bin/merfin/build/bin/merfin -completeness -sequence {input.genome} -readmers {input.reads_meryl} -peak ${{peak}} 2> {output}
        cat {output} | grep "COMPLETENESS:"
        '''

rule AFA_69_merfin_QV_plot:
    input:
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV.hist.log".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV.hist.log".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV.hist.log".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] and config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
    output:
        stat = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}_merfin.QV.stat.txt".format(config["prefix"], config["prefix"]),
        pdf = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}_merfin.QV.pdf".format(config["prefix"], config["prefix"])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        for i in {input}
        do
        cat ${{i}} | awk -F"\\t" -v type=${{i##*/}} 'BEGIN{{OFS="\\t"}} $1 ~ /^Chr/ {{print $1,$5, type}}' | sort -k1,1V >> {output.stat}
        done
        Rscript workflow/scripts/merqury_QV_plot.R {output.stat} {output.pdf}
        '''

rule AFA_70_merfin_completeness_plot:
    input:
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.Illumina.QV.completeness".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.HiFi.QV.completeness".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.hybrid.QV.completeness".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] and config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
    output:
        stat = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.completeness.stat.txt".format(config["prefix"], config["prefix"]),
        pdf = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.completeness.pdf".format(config["prefix"], config["prefix"])
    params:
        prefix = config["prefix"]
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        for i in {input}
        do
        cat ${{i}} | \\grep "COMPLETENESS" | cut -f 2 -d : | sed "s/\\s\\+/${{i##*/}}\\t/g" | sed "s/$/\\t{params.prefix}/" >> {output.stat}
        done
        Rscript workflow/scripts/merqury_QV_plot.R {output.stat} {output.pdf}
        '''

rule AFA_71_Merqury_merfin_QV_plot:
    input:
        Merqury = "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}_Merqury.QV.stat.txt".format(config["prefix"], config["prefix"], config["prefix"]),
        merfin = "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}_merfin.QV.stat.txt".format(config["prefix"], config["prefix"])
    output:
        qv_stat = "results/{}/Assessment_Finally_Assembly/{}_Merqury_Merfin_QV.txt".format(config["prefix"], config["prefix"]),
        qv_plot = "results/{}/Assessment_Finally_Assembly/{}_Merqury_Merfin_QV.pdf".format(config["prefix"], config["prefix"])
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    shell:
        '''
        module load R/3.6.0
        cat {input.Merqury} | sed "s/$/_Merqury/" >> {output.qv_stat}
        cat {input.merfin} | sed "s/$/_merfin/" >> {output.qv_stat}
        Rscript workflow/scripts/merqury_QV_plot.R {output.qv_stat} {output.qv_plot}
        '''
##################################################################################################################################
############################################################ 8. BioNano  #####################################################
##################################################################################################################################

rule AFA_78_NGSgenome_fa2cmap:
    input:
        genome = select_Assessment_input_genome()
    output:
        genome_cmap = "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/{}_{}_0kb_0labels.cmap".format(config["prefix"], os.path.basename(select_Assessment_input_genome()).split(".")[0], config["BioNano"]["enzyme"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        enzyme = config["BioNano"]["enzyme"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/BioNano.yaml"
    shell:
        '''
        echo "#########################################################################################"
        echo "1. fa to cmap"
        echo "#########################################################################################"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/HybridScaffold/12162019/scripts
        perl fa2cmap_multi_color.pl -i {params.workdir}/{input.genome} -e {params.enzyme} 1 -o {params.workdir}/{params.outdir}
        '''

rule AFA_79_BioNano_DeNovo_cogtigs_mapping_NGSgenome:
    input:
        DeNovoCmap = "results/{}/BioNano/DeNovo/Without_reference/contigs/{}_No_Rcmap_refineFinal1/{}_NO_RCMAP_REFINEFINAL1.cmap".format(config["prefix"], config["prefix"], config["prefix"].upper()),
        genome_cmap = "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/{}_{}_0kb_0labels.cmap".format(config["prefix"], os.path.basename(select_Assessment_input_genome()).split(".")[0], config["BioNano"]["enzyme"])
    output:
        alignmentstatistics = "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/BioNano_DeNovo_cogtigs_mapping/{}_alignmentstatistics.out".format(config["prefix"], os.path.basename(select_Assessment_input_genome()).split(".")[0])
    threads: 5
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        deNovo_optArguments_XML = config["BioNano"]["deNovo_optArguments_XML_DLE"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] == "CTTAAG" 
                             else config["BioNano"]["deNovo_optArguments_XML"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] != "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype_DLE"] if config["BioNano"]["assembly_type"] == "haplotype" and config["BioNano"]["enzyme"] == "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype"]
    conda:
        "../envs/BioNano.yaml"
    shell:
        '''
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019
        python runCharacterize.py -o {params.workdir}/{params.outdir} \
            -t {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/RefAligner \
            -q {params.workdir}/{input.DeNovoCmap} \
            -r {params.workdir}/{input.genome_cmap} \
            -p {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019/ \
            -a {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/{params.deNovo_optArguments_XML} \
            -n {threads} 1>{params.workdir}/{output.alignmentstatistics}
        '''

rule AFA_79_BioNano_bnx_molecule_mapping_BioNano_DeNovo_cogtigs:
    input:
        filterbnx_mg = "results/{}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx".format(config["prefix"]),
        DeNovoCmap = "results/{}/BioNano/DeNovo/Without_reference/contigs/{}_No_Rcmap_refineFinal1/{}_NO_RCMAP_REFINEFINAL1.cmap".format(config["prefix"], config["prefix"], config["prefix"].upper())
    output:
        "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/BioNano_bnx_molecule_mapping_BioNano_DeNovo_cogtigs/alignments.tar.gz".format(config["prefix"])
    threads: 5
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        deNovo_optArguments_XML = config["BioNano"]["deNovo_optArguments_XML_DLE"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] == "CTTAAG" 
                             else config["BioNano"]["deNovo_optArguments_XML"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] != "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype_DLE"] if config["BioNano"]["assembly_type"] == "haplotype" and config["BioNano"]["enzyme"] == "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype"]
    conda:
        "../envs/BioNano.yaml"
    shell:
        '''
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019
        python align_bnx_to_cmap.py \
            --mol {params.workdir}/{input.filterbnx_mg} \
            --pipeline {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019/ \
            --ra {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/ \
            --nthreads {threads} \
            --output {params.workdir}/{params.outdir} \
            --ref {params.workdir}/{input.DeNovoCmap} \
            --optArgs {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/{params.deNovo_optArguments_XML} \
            --prefix {params.prefix}
        '''


rule AFA_80_BioNano_bnx_mapping_genome:
    input:
        filterbnx_mg = "results/{prefix}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx",
        genome_cmap = "results/{{prefix}}/Assessment_Finally_Assembly/BioNano_Mapping/{{Chromosome}}_{}_0kb_0labels.cmap".format(config["BioNano"]["enzyme"])
    output:
        "results/{prefix}/Assessment_Finally_Assembly/BioNano_Mapping/BioNano_bnx_molecule_mapping_NGS_genome/{Chromosome}/{Chromosome}_pipelineReport.txt"
    threads: 5
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        enzyme = config["BioNano"]["enzyme"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        deNovo_optArguments_XML = config["BioNano"]["deNovo_optArguments_XML_DLE"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] == "CTTAAG" 
                             else config["BioNano"]["deNovo_optArguments_XML"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] != "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype_DLE"] if config["BioNano"]["assembly_type"] == "haplotype" and config["BioNano"]["enzyme"] == "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype"]
    conda:
        "../envs/BioNano.yaml"
    shell:
        '''
        echo "#########################################################################################"
        echo "2. Align BNX to Reference"
        echo "#########################################################################################"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019
        python align_bnx_to_cmap.py --mol {params.workdir}/{input.filterbnx_mg} --ref {params.workdir}/{input.genome_cmap} \
          --prefix {params.prefix} \
          --pipeline {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019 \
          --optArgs {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/{params.deNovo_optArguments_XML} \
          --ra {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel \
          --nthreads {threads} --output {params.workdir}/{params.outdir} --snrFilter 0 --color 1 > {params.workdir}/{output} 2>&1
        '''

##################################################################################################################################
########################################################### 9. list ###############################################################
##################################################################################################################################

rule AFA_81_list:
    input:
        ## Just Genome fasta file
        "results/{}/Assessment_Finally_Assembly/{}.length".format(config["prefix"], input_genome_name()),
        "results/{}/Assessment_Finally_Assembly/QUAST".format(config["prefix"]),
        expand("results/{prefix}/Assessment_Finally_Assembly/BUSCO/{genome}_{BUSCODB}/short_summary.specific.{BUSCODB}.{genome}.txt", prefix = config["prefix"], genome = os.path.basename(select_Assessment_input_genome()), BUSCODB = config["BUSCO_databsse"]),
        "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats_ACGTN.pdf".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/{}_BBMap_stats.txt".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/{}_seqkit_fx2tab.txt".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/GC_content/{}_nuc_1kb.igv".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Telomeric/{}_genome_telomeric.txt".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Telomeric/Telomere_Identification/{}_telomeric_repeat_windows.svg".format(config["prefix"], config["prefix"]),
        "results/{}/Assessment_Finally_Assembly/Telomeric/vgp_telomeric/{}.summary".format(config["prefix"], input_genome_name().split(".")[0]),
        "results/{}/Assessment_Finally_Assembly/Gap_Stats/{}_genome_Gaps.pdf".format(config["prefix"], config["prefix"]),
        ## T2T
        "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_T2T_Alignment.pdf".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_T2T"] else [],
        expand("results/{prefix}/Assessment_Finally_Assembly/T2T_Alignment/{prefix}_{assdata}2genome_F256.w{windows}.coverage.bedgraph.line.pdf", prefix = config["prefix"], assdata = config["params"]["assessment_t2t_data"], windows = config["params"]["assessment_t2t_windows"]) if config["module"]["Assessment_Finally_Assembly_T2T"] else [],
        ## normal
        "results/{}/Assessment_Finally_Assembly/bowtie2/All_sample_bowtie2.log".format(config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_wgscoverageplotter.svg".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/bowtie2/{}_Illumina_tinycov_covplot.pdf".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_cov.txt".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_wgscoverageplotter.svg".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_tinycov_covplot.png".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pb_qualimap".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_cov.txt".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_ONT"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_wgscoverageplotter.svg".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_ONT"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_tinycov_covplot.pdf".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_ONT"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONT_qualimap".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_ONT"] else [],
        "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_Illumina_PB_ONT_covplot.pdf".format(config["prefix"], config["prefix"]),
        ## Merqury
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.Illumina.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] else [],
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.HiFi.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}.hybrid.QV/{}.qv".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] and config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merqury_QV/{}_Merqury.QV.pdf".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] or config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        ## Merfin
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}_merfin.QV.pdf".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] or config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        "results/{}/Assessment_Finally_Assembly/Merfin_QV/{}.completeness.pdf".format(config["prefix"], config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] or config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        ## Merqury + Merfin
        "results/{}/Assessment_Finally_Assembly/{}_Merqury_Merfin_QV.pdf".format(config["prefix"], config["prefix"]) if config["module"]["Assessment_Finally_Assembly_Illumina"] or config["module"]["Assessment_Finally_Assembly_HiFi"] else [],
        ## BioNano
        "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/BioNano_DeNovo_cogtigs_mapping/{}_alignmentstatistics.out".format(config["prefix"], os.path.basename(select_Assessment_input_genome()).split(".")[0]) if config["module"]["Assessment_Finally_Assembly_BioNano"] else [],
        "results/{}/Assessment_Finally_Assembly/BioNano_Mapping/BioNano_bnx_molecule_mapping_BioNano_DeNovo_cogtigs/alignments.tar.gz".format(config["prefix"]) if config["module"]["Assessment_Finally_Assembly_BioNano"] else []
