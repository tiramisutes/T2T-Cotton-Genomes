###################################################################################################################################
# ----- 1. Trimming and Quality control ----- #
###################################################################################################################################
rule ATAC_1_cp_fastq_pe:
    input:
        get_lanes
    output:
        fastq1=temp("results/{prefix}/temp_fastq/{sample}-{lane}.1.fastq.gz"),
        fastq2=temp("results/{prefix}/temp_fastq/{sample}-{lane}.2.fastq.gz")
    shell:
        """
        ln -s {input[0]} {output.fastq1}
        ln -s {input[1]} {output.fastq2}
        """

rule ATAC_2_mergeFastq_pe:
    input:
        fw = lambda w: expand("results/{prefix}/temp_fastq/{lane.sample}-{lane.lane}.1.fastq.gz", prefix = config["prefix"], lane=units.loc[w.sample].itertuples()),
        rv = lambda w: expand("results/{prefix}/temp_fastq/{lane.sample}-{lane.lane}.2.fastq.gz", prefix = config["prefix"], lane=units.loc[w.sample].itertuples())
    output:
        fastq1 = temp("results/{prefix}/temp_fastq/{sample}.1.fastq.gz"),
        fastq2 = temp("results/{prefix}/temp_fastq/{sample}.2.fastq.gz")
    shell:
        """
        cat {input.fw} > {output.fastq1}
        cat {input.rv} > {output.fastq2}
        """

rule ATAC_3_fastp_pe_trimming:
    input:
        fw = "results/{prefix}/temp_fastq/{sample}.1.fastq.gz",
        rv = "results/{prefix}/temp_fastq/{sample}.2.fastq.gz"
    output:
        fastq1 = "results/{prefix}/Step1_fastq_trimmed/{sample}.1.fastq.gz",
        fastq2 = "results/{prefix}/Step1_fastq_trimmed/{sample}.2.fastq.gz",
        json = "results/{prefix}/Step1_fastq_trimmed/{sample}.json",
        html = "results/{prefix}/Step1_fastq_trimmed/{sample}.html"
    threads:
        10
    conda:
        "../envs/ATAC_seq.yaml"
    params:
        workdir = config["workdir"],
        fastp_params = config["params"]["fastp"]["pe"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        """
        fastp -i {input.fw} \
            -I {input.rv} \
            -o {output.fastq1} \
            -O {output.fastq2} \
            -w {threads} \
            {params.fastp_params} \
            -j {output.json} -h {output.html}
        """

###################################################################################################################################
# ----- 2. Aligning with BOWTIE2 and processing with SAMTOOLS ----- #
###################################################################################################################################
rule ATAC_4_bowtie2_index:
    input:
        genome = config["DATA"]["genome_fasta"],
        genome_gtf = config["DATA"]["genome_gtf"]
    output:
        indexes = "results/{}/Input_Genome_Data/{}.fasta".format(config["prefix"], config["DATA"]["genome_name"]),
        autosome_name_file = "results/{}/Input_Genome_Data/{}_autosome_name.txt".format(config["prefix"], config["DATA"]["genome_name"]),
        genome_gtf = "results/{}/Input_Genome_Data/{}.gz".format(config["prefix"], os.path.basename(config["DATA"]["genome_gtf"])),
        genome_size = "results/{}/Input_Genome_Data/{}.size".format(config["prefix"], config["DATA"]["genome_name"]),
        genome_size2 = "results/{}/Input_Genome_Data/{}_chr.size".format(config["prefix"], config["DATA"]["genome_name"])
    params:
        bowtie2 = config["software"]["bowtie2"],
        autosome_name = config["DATA"]["chromosome"],
        genome_gtf_name = lambda w, input: os.path.basename(input["genome_gtf"]),
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 2
    shell:
        '''
        {params.bowtie2}
        cp {input.genome} {output.indexes}
        bowtie2-build {output.indexes} {output.indexes}
        echo {params.autosome_name} | sed "s/ /\\n/g" > {output.autosome_name_file}
        ## 
        cp {input.genome_gtf} {params.outdir}
        gzip {params.outdir}/{params.genome_gtf_name}
        ## genome size
        bioawk -c fastx '{{ print $name, length($seq) }}' < {input.genome} | sort -k1,1V > {output.genome_size}
        bioawk -c fastx '{{ print $name, "1", length($seq) }}' < {input.genome} | sort -k1,1V | grep -i -v "scaffold" | sed '1i chr\\tstart\\tend' > {output.genome_size2}
        '''

rule ATAC_5_bowtie2_align_remove:
    input:
        fastq1 = "results/{prefix}/Step1_fastq_trimmed/{sample}.1.fastq.gz",
        fastq2 = "results/{prefix}/Step1_fastq_trimmed/{sample}.2.fastq.gz",
        indexes = "results/{}/Input_Genome_Data/{}.fasta".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        sortbam   = temp("results/{prefix}/Step2_bowtie2_align/{sample}.sorted.bam"),
        sortbamindex   = temp("results/{prefix}/Step2_bowtie2_align/{sample}.sorted.bam.bai"),
        markedbam = temp("results/{prefix}/Step2_bowtie2_align/{sample}.marked.bam"),
        markedbamindex = temp("results/{prefix}/Step2_bowtie2_align/{sample}.marked.bai"),
        bam   = "results/{prefix}/Step2_bowtie2_align/{sample}.bam",
        index = "results/{prefix}/Step2_bowtie2_align/{sample}.bam.bai"
    threads:
        10
    params:
        bowtie2SF = config["software"]["bowtie2"],
        samblasterSF = config["software"]["samblaster"],
        bowtie2_global = config["params"]["bowtie2"]["global"],
        bowtie2_pe     = config["params"]["bowtie2"]["pe"],
        samblaster   = config["params"]["samblaster"],
        samtools_mem = config["params"]["samtools"]["memory"],
        chrM_name = config["DATA"]["chrM_name"]
    log:
       align   = "results/{prefix}/Step2_bowtie2_align/{sample}_bowtie2_align.log",
       rm_dups = "results/{prefix}/Step2_bowtie2_align/{sample}_rm_dup.log"
    shell:
        """
        {params.bowtie2SF}
        {params.samblasterSF}
        # -q2 removes multimapping, -F 4 removes unmapped, -f 2 removes unpaired reads.
        # fgrep to remove reads mapping to chromosome M 
        bowtie2 -p {threads} --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x {input.indexes} -1 {input.fastq1} -2 {input.fastq2}  2> {log.align} \
            | samtools view -@ {threads} -bS - \
            | samtools sort -@ {threads} - -o {output.sortbam}
        samtools index -@ {threads} {output.sortbam}

        module load picard/2.23.9
        java -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates QUIET=true INPUT={output.sortbam} OUTPUT={output.markedbam} METRICS_FILE={log.rm_dups} REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
        samtools view -@ {threads} -h -b -f 2 -F 1548 -q 30 {output.markedbam} | samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule ATAC_6_bam2bigwig_noSubstract:
    input: 
        "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output:  
        "results/{prefix}/Step3_Sample_bigwig/{sample}.bw"
    params: 
        deeptools = config["software"]["deeptools"],
        params = config["params"]["bam2bigwig"]
    threads: 
        10
    shell:
        """
        {params.deeptools}
        bamCoverage -b {input} -o {output} --numberOfProcessors {threads}
        """

###################################################################################################################################
# ----- 3. Peak-calling with MACS2 ----- #
###################################################################################################################################
rule ATAC_7_Shift_Read_Coordinates:
    """
    alignmentSieve: If the --shift or --ATACshift options are used, then only properly-paired reads will be used.
    """
    input:
        bam = "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output:
        bam_temp = temp("results/{prefix}/Step2_bowtie2_align/{sample}.dt.bam"),
        shifted_bam = "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bam",
        shifted_bam_index = "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bam.bai"
    threads: 2
    conda:
         "../envs/ATAC_deeptools.yaml"
    shell:
        """
        alignmentSieve -b {input.bam} -p {threads} --samFlagInclude 2 --ATACshift -o {output.bam_temp}
        samtools sort -@ {threads} -o {output.shifted_bam} {output.bam_temp}
        samtools index -@ {threads} {output.shifted_bam}
        """

rule ATAC_8_Shift_Read_bam2bed:
    input:
        shifted_bam = "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bam",
        #shifted_bam_index = "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bam.bai"
    output: 
        "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bed"
    threads: 10
    params:
        deeptools = config["software"]["deeptools"]
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        macs2 randsample -i {input.shifted_bam} -f BAMPE -p 100 -o {output}
        """

rule ATAC_9_call_peaks_macs2_Individual:
    input:
        "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bed"
    output: 
        bed = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.{peak_shapes}Peak",
        peaks = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.xls",
        summits = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_summits.bed",
        bdgt = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_treat_pileup.bdg",
        bdgtc = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_control_lambda.bdg"
    params:
        prefix = "{sample}_macs2",
        peak_shapes = "{peak_shapes}",
        macs2_pvalue = config["params"]["macs2_pvalue"],
        macs2_smooth_win = config["params"]["macs2_smooth_win"],
        macs2_pvalue_broad = config["params"]["macs2_pvalue_broad"],
        genome_size = config["params"]["egenome_size"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        if [ "{params.peak_shapes}" == 'narrow' ]
        then
            echo -e "####################################################################### \\nnarrow peaks is call. \\n#######################################################################"
            macs2 callpeak -t {input} -f BED -n {params.prefix} -g {params.genome_size} -p {params.macs2_pvalue} --outdir {params.outdir} \
                --shift -`echo "{params.macs2_smooth_win}/2" | bc` --extsize {params.macs2_smooth_win} --nomodel -B --SPMR --keep-dup all --call-summits
        elif [ "{params.peak_shapes}" == 'broad' ]
        then
            echo -e "####################################################################### \\broad peaks is call. \\n#######################################################################"
            macs2 callpeak -t {input} -f BED -n {params.prefix} -g {params.genome_size} -p {params.macs2_pvalue} --outdir {params.outdir} \
                --broad --shift -`echo "{params.macs2_smooth_win}/2" | bc` --extsize {params.macs2_smooth_win} --nomodel -B --SPMR --keep-dup all --call-summits
        else
            echo -e "####################################################################### \\nnarrow or broad should be set on macs2_peak_shapes in config.yaml file. \\n#######################################################################"
        fi
        """

rule ATAC_Step1_call_peaks_macs2_list:
    input:
        expand("results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.{peak_shapes}Peak", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------- Biological Replicates ------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if config["params"]["Biological_Replicates"]:
    rule ATAC_11_call_peaks_macs2_RepsJoin:
        input:
            lambda w: expand("results/{prefix}/Step2_bowtie2_align/{condition}.shifted.bed", prefix= config["prefix"], condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
        output: 
            bed = "results/{prefix}/Step4_call_peaks/MACS2_Condition_{peak_shapes}Peak/{condition}_macs2_peaks.{peak_shapes}Peak",
            peaks = "results/{prefix}/Step4_call_peaks/MACS2_Condition_{peak_shapes}Peak/{condition}_macs2_peaks.xls",
            summits = "results/{prefix}/Step4_call_peaks/MACS2_Condition_{peak_shapes}Peak/{condition}_macs2_summits.bed",
            bdgt = "results/{prefix}/Step4_call_peaks/MACS2_Condition_{peak_shapes}Peak/{condition}_macs2_treat_pileup.bdg",
            bdgtc = "results/{prefix}/Step4_call_peaks/MACS2_Condition_{peak_shapes}Peak/{condition}_macs2_control_lambda.bdg"
        params:
            prefix = "{condition}_macs2",
            peak_shapes = "{peak_shapes}",
            macs2_smooth_win = config["params"]["macs2_smooth_win"],
            macs2_pvalue = config["params"]["macs2_pvalue"],
            macs2_pvalue_broad = config["params"]["macs2_pvalue_broad"],
            genome_size = config["params"]["egenome_size"],
            outdir = lambda w, output: os.path.dirname(output[0])
        conda:
            "../envs/ATAC_seq.yaml"
        shell:
            """
            if [ "{params.peak_shapes}" == 'narrow' ]
            then
                echo -e "####################################################################### \\nnarrow peaks is call. \\n#######################################################################"
                macs2 callpeak -t {input} -f BED -n {params.prefix} -g {params.genome_size} -p {params.macs2_pvalue} --outdir {params.outdir} \
                    --shift -`echo "{params.macs2_smooth_win}/2" | bc` --extsize {params.macs2_smooth_win} --nomodel -B --SPMR --keep-dup all --call-summits
            elif [ "{params.peak_shapes}" == 'broad' ]
            then
                echo -e "####################################################################### \\broad peaks is call. \\n#######################################################################"
                macs2 callpeak -t {input} -f BED -n {params.prefix} -g {params.genome_size} -p {params.macs2_pvalue} --outdir {params.outdir} \
                    --broad --shift -`echo "{params.macs2_smooth_win}/2" | bc` --extsize {params.macs2_smooth_win} --nomodel -B --SPMR --keep-dup all --call-summits
            else
                echo -e "####################################################################### \\nnarrow or broad should be set on macs2_peak_shapes in config.yaml file. \\n#######################################################################"
            fi
            """

rule ATAC_13_Downstream_Analysis_Peaks:
    input:
        "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak_Overlap_Condition_{peak_shapes}Peak_Overlap/sets/11_{condition}_idr_{condition}_macs2_peaks.bed" 
          if config["params"]["Biological_Replicates"] else "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.{peak_shapes}Peak"
    output:
        "results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{condition}_macs2_peaks.{peak_shapes}Peak" 
          if config["params"]["Biological_Replicates"] else "results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{sample}_macs2_peaks.{peak_shapes}Peak"
    shell:
        '''
        cp {input} {output}
        '''

rule ATAC_Step2_Downstream_Analysis_Peaks_list:
    input:
        expand("results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{condition}_macs2_peaks.{peak_shapes}Peak", prefix= config["prefix"], condition=ALL_CONDITIONS, peak_shapes = config["params"]["macs2_peak_shapes"]) 
          if config["params"]["Biological_Replicates"] else expand("results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{sample}_macs2_peaks.{peak_shapes}Peak", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"])

###################################################################################################################################
# ----- 4. Repetition and Peak Quality Control ----- #
###################################################################################################################################
rule ATAC_14_reads_fastqc:
    input:
        expand("results/{prefix}/Step1_fastq_trimmed/{sample}.{group}.fastq.gz", prefix = config["prefix"], group=[1, 2], sample = ALL_SAMPLES)
    output:
        expand("results/{prefix}/Step5_QCs/FastQC/{sample}.{group}_fastqc.zip", prefix = config["prefix"], group=[1, 2], sample = ALL_SAMPLES)
    params:
        FastQC = config["software"]["FastQC"],
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 10
    shell:
        '''
        {params.FastQC}
        for i in {input}
        do
        fastqc -o {params.outdir} -t {threads} ${{i}}
        done
        '''

# ------- Aligment QC ------- #
rule ATAC_15_mapping_qc:
    input:
        "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output:
        #stat = "results/{prefix}/Step5_QCs/Distribution_Mapped_Reads/{sample}_Distribution_Mapped_Reads.stat",
        pdf = "results/{prefix}/Step5_QCs/Distribution_Mapped_Reads/{sample}_Distribution_Mapped_Reads.pdf"
    params:
        R = config["software"]["R"]
    threads: 1
    shell:
        '''
        {params.R}
        Rscript workflow/scripts/AlignemrQC.R {input} {output.pdf}
        '''

# ------- InsertSize calculation ------- #
rule ATAC_16_insert_size:
    input:
        "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output:
        txt="results/{prefix}/Step5_QCs/insert_size/{sample}.isize.txt",
        pdf="results/{prefix}/Step5_QCs/insert_size/{sample}.isize.pdf"
    params:
        picardSF = config["software"]["picard"],
        R = config["software"]["R"],
        picardPA = "VALIDATION_STRINGENCY=LENIENT METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
    shell:
        """
        {params.picardSF}
        {params.R}
        java -jar ${{EBROOTPICARD}}/picard.jar CollectInsertSizeMetrics {params.picardPA} INPUT={input} OUTPUT={output.txt} HISTOGRAM_FILE={output.pdf}
        """

# ------- Calculate FRiP score, The higher the score, the better. ------- #
rule ATAC_17_calculate_FRiP:
    input:
        shifted_bam = "results/{prefix}/Step2_bowtie2_align/{sample}.shifted.bam",
        peaks = "results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.{peak_shapes}Peak"
    output:
        peaks_count = temp("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{peak_shapes}Peak.count_mqc.tsv"),
        peaks_FRiP = temp("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.tsv")
    threads: 5
    params:
        deeptools = config["software"]["deeptools"],
        sample = "{sample}"
    shell:
        '''
        {params.deeptools}
        cat {input.peaks} | wc -l | awk -v OFS='\\t' '{{ print "{params.sample}", $1 }}' > {output.peaks_count}
        python workflow/scripts/calculate_FRiP_score.py -b {input.shifted_bam} -p {input.peaks} -t {threads} -o {output.peaks_FRiP}
        '''

rule ATAC_18_calculate_FRiP_score_plot:
    input: 
        peaks_count = expand("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{{peak_shapes}}Peak.count_mqc.tsv", prefix= config["prefix"], sample = ALL_SAMPLES),
        peaks_FRiP = expand("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{{peak_shapes}}Peak.FRiP_mqc.tsv", prefix= config["prefix"], sample = ALL_SAMPLES),
        samples = config["samples"]
    output: 
        peaks_count = "results/{prefix}/Step5_QCs/FRiP_score/{prefix}_All_sample_macs2_peaks.{peak_shapes}Peak.count_mqc.tsv",
        peaks_count_plot = "results/{prefix}/Step5_QCs/FRiP_score/{prefix}_All_sample_macs2_peaks.{peak_shapes}Peak.count_mqc.pdf",
        stat = "results/{prefix}/Step5_QCs/FRiP_score/{prefix}_All_sample_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.tsv",
        stat_plot = "results/{prefix}/Step5_QCs/FRiP_score/{prefix}_All_sample_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.pdf"
    params:
        R = config["software"]["R"]
    shell:
        """
        {params.R}
        echo -e "Sample\\tPeaks_Count" | cat - {input.peaks_count} > {output.peaks_count}
        Rscript workflow/scripts/Peaks_Count_Plot.R {output.peaks_count} {output.peaks_count_plot}
        echo -e "Sample\\tPeaks\\tFRiP" | cat - {input.peaks_FRiP} > {output.stat}
        Rscript workflow/scripts/calculate_FRiP_score_plot.R {output.stat} {input.samples} {output.stat_plot}
        """


rule ATAC_19_replicates_correlation:
    input:
        lambda w: expand("results/{prefix}/Step3_Sample_bigwig/{condition}.bw", prefix= config["prefix"], condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME)
    output:
        "results/{prefix}/Step5_QCs/RepeatabilityCorrelate/{condition}_Correlate.txt"
    params:
        deeptools = config["software"]["deeptools"],
        indir = lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/ATAC_seq.yaml"
    threads: 5
    shell:
        """
        {params.deeptools}
        wigCorrelate {input[0]} {input[1]} > {output}
        """

rule ATAC_20_replicates_correlation_summary:
    input:
        expand("results/{prefix}/Step5_QCs/RepeatabilityCorrelate/{condition}_Correlate.txt", prefix= config["prefix"], condition=ALL_CONDITIONS)
    output:
        "results/{prefix}/Step5_QCs/RepeatabilityCorrelate/{prefix}_All_Sample_Correlate.txt"
    shell:
        """
        cat {input} | sed '1i sample_rep1\\tsample_rep2\\treplicates_correlation' > {output}
        """

rule ATAC_21_BAM_multiBamSummary:
    input:
        expand("results/{prefix}/Step2_bowtie2_align/{sample}.bam", prefix= config["prefix"], sample = ALL_SAMPLES)
    output:
        readCounts_npz = "results/{}/Step5_QCs/RepeatabilityTest/{}_readCounts.npz".format(config["prefix"], config["prefix"]),
        readCounts_tab = "results/{}/Step5_QCs/RepeatabilityTest/{}_readCounts.tab".format(config["prefix"], config["prefix"])
    params:
        deeptools = config["software"]["deeptools"],
        indir = lambda w, input: os.path.dirname(input[0])
    threads: 5
    shell:
        """
        {params.deeptools}
        input_sample_name=`echo {input} | sed "s/ /\\n/g" | awk -F"/" '{{print $NF}}' | awk -F"." '{{print $1}}' | xargs`
        multiBamSummary bins \
            --bamfiles {input} \
            --minMappingQuality 30 \
            --labels ${{input_sample_name}} \
            --numberOfProcessors {threads} \
            -out {output.readCounts_npz} --outRawCounts {output.readCounts_tab}
        """

rule ATAC_22_BAM_plotCorrelation:
    input:
        readCounts_npz = "results/{}/Step5_QCs/RepeatabilityTest/{}_readCounts.npz".format(config["prefix"], config["prefix"])
    output:
        plotCorrelation_heatmap_SpearmanCorr_pdf = "results/{}/Step5_QCs/RepeatabilityTest/{}_heatmap_SpearmanCorr_readCounts.pdf".format(config["prefix"], config["prefix"]),
        plotCorrelation_SpearmanCorr_tab = "results/{}/Step5_QCs/RepeatabilityTest/{}_SpearmanCorr_readCounts.tab".format(config["prefix"], config["prefix"]),
        plotCorrelation_heatmap_PearsonCorr_pdf = "results/{}/Step5_QCs/RepeatabilityTest/{}_heatmap_PearsonCorr_readCounts.pdf".format(config["prefix"], config["prefix"]),
        plotCorrelation_PearsonCorr_tab = "results/{}/Step5_QCs/RepeatabilityTest/{}_PearsonCorr_readCounts.tab".format(config["prefix"], config["prefix"]),
        plotCorrelation_scatterplot_SpearmanCorr_pdf = "results/{}/Step5_QCs/RepeatabilityTest/{}_scatterplot_SpearmanCorr_readCounts.pdf".format(config["prefix"], config["prefix"]),
        plotCorrelation_scatterplot_PearsonCorr_pdf = "results/{}/Step5_QCs/RepeatabilityTest/{}_scatterplot_PearsonCorr_readCounts.pdf".format(config["prefix"], config["prefix"]),
        PCA_readCounts_pdf = "results/{}/Step5_QCs/RepeatabilityTest/{}_PCA_readCounts.pdf".format(config["prefix"], config["prefix"])
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        plotCorrelation \
            -in {input} \
            --corMethod spearman --skipZeros \
            --plotTitle "Spearman Correlation of Read Counts" \
            --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
            -o {output.plotCorrelation_heatmap_SpearmanCorr_pdf} \
            --outFileCorMatrix {output.plotCorrelation_SpearmanCorr_tab}

        plotCorrelation \
            -in {input} \
            --corMethod pearson --skipZeros \
            --plotTitle "Pearson Correlation of Read Counts" \
            --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
            -o {output.plotCorrelation_heatmap_PearsonCorr_pdf} \
            --outFileCorMatrix {output.plotCorrelation_PearsonCorr_tab}

        plotCorrelation \
            -in {input} \
            --corMethod spearman --skipZeros \
            --plotTitle "Spearman Correlation of Average Scores Per Transcript" \
            --whatToPlot scatterplot \
            -o {output.plotCorrelation_scatterplot_SpearmanCorr_pdf}
        
        plotCorrelation \
            -in {input} \
            --corMethod pearson --skipZeros \
            --plotTitle "Pearson Correlation of Average Scores Per Transcript" \
            --whatToPlot scatterplot \
            -o {output.plotCorrelation_scatterplot_PearsonCorr_pdf}
        
        plotPCA -in {input} \
            -o {output.PCA_readCounts_pdf} \
            -T "PCA of read counts"
        """

rule ATAC_23_Count_plotCorrelation_PCA_R:
    input:
        readCounts_tab = "results/{}/Step5_QCs/RepeatabilityTest/{}_readCounts.tab".format(config["prefix"], config["prefix"]),
        samples = config["samples"]
    output:
        readCounts_tsv = "results/{}/Step5_QCs/RepeatabilityTest/{}_readCounts.tab.tsv".format(config["prefix"], config["prefix"]),
        readCounts_plot = "results/{}/Step5_QCs/RepeatabilityTest/{}_Correlation_PCA.pdf".format(config["prefix"], config["prefix"]),
        readCounts_plot2 = "results/{}/Step5_QCs/RepeatabilityTest/{}_Spearman_Correlation.tiff".format(config["prefix"], config["prefix"]),
        readCounts_plot3 = "results/{}/Step5_QCs/RepeatabilityTest/{}_Spearman_Correlation.pdf".format(config["prefix"], config["prefix"])
    params:
        workdir = config["workdir"],
        R = config["software"]["R"],
        gcc = config["software"]["gcc"]
    shell:
        '''
        {params.R}
        {params.gcc}
        cat {input.readCounts_tab} | sed "s/'//g" | sed "s/#//g" > {output.readCounts_tsv}
        Rscript workflow/scripts/plot-cor-pca.R {output.readCounts_tsv} {input.samples} {output.readCounts_plot} {output.readCounts_plot2} {output.readCounts_plot3}
        '''

# ------- Deeptools quality control ------- #

rule ATAC_24_plotFingerprint:
    input: 
        "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output: 
        qualMetrics = "results/{prefix}/Step5_QCs/fingerPrint/{sample}.qualityMetrics.tsv",
        raw_counts  = "results/{prefix}/Step5_QCs/fingerPrint/{sample}.rawcounts.tsv",
        plot        = "results/{prefix}/Step5_QCs/fingerPrint/{sample}.plot.pdf"
    params:
        read_exten = "--extendReads",
        deeptools = config["software"]["deeptools"]
    threads:
        10
    shell:
        """
        {params.deeptools}
        plotFingerprint --bamfiles {input} \
            -p {threads} \
            --skipZeros \
            --outQualityMetrics {output.qualMetrics} \
            --outRawCounts {output.raw_counts} \
            --plotFile {output.plot}
        """

# ------- Run ataqv on BAM file and corresponding peaks ------- #
rule ATAC_25_ataqv_json:
    input:
        peaks_bed = "results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{condition}_macs2_peaks.{peak_shapes}Peak" 
          if config["params"]["Biological_Replicates"] else "results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{sample}_macs2_peaks.{peak_shapes}Peak",
        bam = lambda w: expand("results/{prefix}/Step2_bowtie2_align/{condition}.bam", prefix = config["prefix"], condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME) if config["params"]["Biological_Replicates"] else "results/{prefix}/Step2_bowtie2_align/{sample}.bam",
        autosome_name_file = "results/{}/Input_Genome_Data/{}_autosome_name.txt".format(config["prefix"], config["DATA"]["genome_name"]),
        tss_bed = "results/{}/Input_Genome_Data/{}_TSS.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        ataqv_json = "results/{prefix}/Step5_QCs/ataqv/{condition}_{peak_shapes}peaks.ataqv.json" 
          if config["params"]["Biological_Replicates"] else "results/{prefix}/Step5_QCs/ataqv/{sample}_{peak_shapes}peaks.ataqv.json"
    threads: 10
    params:
        prefix = config["prefix"],
        chrM_name = config["DATA"]["chrM_name"]
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        ataqv --threads {threads} --autosomal-reference-file {input.autosome_name_file} --mitochondrial-reference-name {params.chrM_name} \
          --peak-file {input.peaks_bed} --tss-file {input.tss_bed} --ignore-read-groups \
          --metrics-file {output.ataqv_json} {params.prefix} {input.bam}
        """

rule ATAC_26_ataqv_json_to_html:
    input:
        expand("results/{prefix}/Step5_QCs/ataqv/{condition}_{{peak_shapes}}peaks.ataqv.json", prefix= config["prefix"], condition=ALL_CONDITIONS) 
          if config["params"]["Biological_Replicates"] else expand("results/{prefix}/Step5_QCs/ataqv/{sample}_{peak_shapes}peaks.ataqv.json", prefix= config["prefix"], sample = ALL_SAMPLES)
    output:
        "results/{prefix}/Step5_QCs/ataqv/{prefix}_All_Sample_ATAC_{peak_shapes}peaks_QC_html/index.html"
    params:
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        mkarv --force {params.outdir} {input}
        """

# ------- macs peak ------- #
rule ATAC_27_macs_peak_condition_QC:
    input:
        expand("results/{prefix}/Step4_call_peaks/MACS2_Downstream_Analysis_Peaks/{condition}_macs2_peaks.{{peak_shapes}}Peak", prefix= config["prefix"], condition=ALL_CONDITIONS)
    output:
        summary = "results/{prefix}/Step5_QCs/{prefix}_All_Condition_macs_{peak_shapes}peak.summary.txt",
        plot = "results/{prefix}/Step5_QCs/{prefix}_All_Condition_macs_{peak_shapes}peak.plots.pdf"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        R = config["software"]["R"],
        gcc = config["software"]["gcc"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        indir = lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        {params.R}
        {params.gcc}
        peaks_bed=`echo {input} | sed "s/ /,/g"`
        peaks_bed_name=`echo {input} | sed "s|{params.indir}/||g" | sed "s/ /,/g"`
        Rscript workflow/scripts/plot_macs_qc.R -i ${{peaks_bed}} -s ${{peaks_bed_name}} -o {params.outdir} -p $(basename {output.summary} ".summary.txt")
        """

rule ATAC_28_macs_peak_sample_QC:
    input:
        expand("results/{prefix}/Step4_call_peaks/MACS2_Sample_{{peak_shapes}}Peak/{sample}_macs2_peaks.{{peak_shapes}}Peak", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"])
    output:
        summary = "results/{prefix}/Step5_QCs/{prefix}_All_Sample_macs_{peak_shapes}peak.summary.txt",
        plot = "results/{prefix}/Step5_QCs/{prefix}_All_Sample_macs_{peak_shapes}peak.plots.pdf"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        R = config["software"]["R"],
        gcc = config["software"]["gcc"],
        outdir = lambda w, output: os.path.dirname(output[0]),
        indir = lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/ATAC_seq.yaml"
    shell:
        """
        {params.R}
        {params.gcc}
        peaks_bed=`echo {input} | sed "s/ /,/g"`
        peaks_bed_name=`echo {input} | sed "s|{params.indir}/||g" | sed "s/ /,/g"`
        Rscript workflow/scripts/plot_macs_qc.R -i ${{peaks_bed}} -s ${{peaks_bed_name}} -o {params.outdir} -p $(basename {output.summary} ".summary.txt")
        """

# ------- SPOT (signal-to-noise ratio) ------- #
rule ATAC_29_signal_to_noise_SPOT_prepare:
    input:
        indexes = "results/{}/Input_Genome_Data/{}.fasta".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mappable_temp = temp("results/{}/Input_Genome_Data/{}.fasta.K{}.mappable_only.bed".format(config["prefix"], config["DATA"]["genome_name"], config["params"]["hotspot_kmer_size"])),
        mappable_sort = "results/{}/Input_Genome_Data/{}.fasta.K{}.mappable_only_sorted.bed".format(config["prefix"], config["DATA"]["genome_name"], config["params"]["hotspot_kmer_size"]),
        chromInfo = "results/{}/Input_Genome_Data/{}.chromInfo.bed".format(config["prefix"], config["DATA"]["genome_name"])
    params:
        workdir = config["workdir"],
        hotspot_kmer_size = config["params"]["hotspot_kmer_size"],
        indir = lambda w, input: os.path.dirname(input["indexes"])
    threads: 1
    shell:
        '''
        module load Bowtie/1.3.1
        module load BEDTools/2.27
        module load GCC/7.2.0-2.29
        module load GSL/2.4
        module load R/3.6.0
        module load bedops/2.4.39
        cd {params.indir}
        {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/hotspot-deploy/bin/enumerateUniquelyMappableSpace `basename {input.indexes}` {params.hotspot_kmer_size}
        sort-bed {params.workdir}/{output.mappable_temp} > {params.workdir}/{output.mappable_sort}
        {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/hotspot-deploy/bin/writeChromInfoBed.pl `basename {input.indexes}` > {params.workdir}/{output.chromInfo}
        '''

rule ATAC_32_signal_to_noise_SPOT_run:
    input:
        indexes = "results/{}/Input_Genome_Data/{}.fasta".format(config["prefix"], config["DATA"]["genome_name"]),
        mappable = "results/{}/Input_Genome_Data/{}.fasta.K{}.mappable_only_sorted.bed".format(config["prefix"], config["DATA"]["genome_name"], config["params"]["hotspot_kmer_size"]),
        chromInfo = "results/{}/Input_Genome_Data/{}.chromInfo.bed".format(config["prefix"], config["DATA"]["genome_name"]),
        bam = "results/{prefix}/Step2_bowtie2_align/{sample}.bam"
    output:
        tokens = "results/{prefix}/Step5_QCs/Signal_to_Noise_SPOT/{sample}/{sample}.tokens.txt",
        runhotspot = "results/{prefix}/Step5_QCs/Signal_to_Noise_SPOT/{sample}/{sample}.runhotspot",
        spot = "results/{prefix}/Step5_QCs/Signal_to_Noise_SPOT/{sample}/{sample}.spot.out"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        hotspot_kmer_size = config["params"]["hotspot_kmer_size"],
        chromosome = config["DATA"]["chromosome"],
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 1
    shell:
        '''
        module load Bowtie/1.3.1
        module load BEDTools/2.27
        module load GCC/7.2.0-2.29
        module load GSL/2.4
        module load R/3.6.0
        module load bedops/2.4.39
        module load Python/2.7.15
        cd {params.outdir}
        chromosome=`echo {params.chromosome} | awk '{{print $NF}}'`
        cat {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/pipeline-scripts/test/runall.tokens.txt \
          | sed "s|^_TAGS_ = .*|_TAGS_ = {params.workdir}/{input.bam}|" \
          | sed "s/^_GENOME_ = .*/_GENOME_ = {params.prefix}/" | sed "s/^_K_ = .*/_K_ = {params.hotspot_kmer_size}/" \
          | sed "s|^_CHROM_FILE_ = .*|_CHROM_FILE_ = {params.workdir}/{input.chromInfo}|" \
          | sed "s|^_MAPPABLE_FILE_ = .*|_MAPPABLE_FILE_ = {params.workdir}/{input.mappable}|" \
          | sed "s/^_DUPOK_ = .*/_DUPOK_ = F/" | sed "s/^_FDRS_ = .*/_FDRS_ = \\"N\\"/" | sed "s/^_CHKCHR_ = .*/_CHKCHR_ = ${{chromosome}}/"\
          | sed "s|^_OUTDIR_ = .*|_OUTDIR_ = {params.workdir}/{params.outdir}|" \
          | sed "s|^_RANDIR_ = .*|_RANDIR_ = {params.workdir}/{params.outdir}|" \
          | sed "s/^_OMIT_REGIONS_: .*/_OMIT_REGIONS_: /" \
          | sed "s|^_HOTSPOT_ = .*|_HOTSPOT_ = {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/hotspot-deploy/bin/hotspot|" \
          | sed "s|^_PKFIND_BIN_ = .*|_PKFIND_BIN_ = {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/hotspot-deploy/bin/wavePeaks|" > {params.workdir}/{output.tokens}
        cat {params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/pipeline-scripts/test/runhotspot \
          | sed "s|^scriptTokBin=.*|scriptTokBin={params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/ScriptTokenizer/src/script-tokenizer.py|" \
          | sed "s|^pipeDir=.*|pipeDir={params.workdir}/workflow/bin/hotspot-4.1.1/hotspot-distr/pipeline-scripts|" \
          | sed "s|^tokenFile=.*|tokenFile={params.workdir}/{output.tokens}|" > {params.workdir}/{output.runhotspot}
        chmod +x {params.workdir}/{output.runhotspot}
        bash {params.workdir}/{output.runhotspot}
        '''

rule ATAC_33_signal_to_noise_SPOT_plot:
    input:
        spot = expand("results/{prefix}/Step5_QCs/Signal_to_Noise_SPOT/{sample}/{sample}.spot.out", prefix= config["prefix"], sample = ALL_SAMPLES),
        samples = config["samples"]
    output:
        spots = "results/{}/Step5_QCs/Signal_to_Noise_SPOT/{}_All_Sample_SPOT.txt".format(config["prefix"], config["prefix"]),
        spots_plot = "results/{}/Step5_QCs/Signal_to_Noise_SPOT/{}_All_Sample_SPOT.pdf".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        R = config["software"]["R"]
    shell:
        '''
        {params.R}
        for i in {input.spot}
        do
        sample=`basename ${{i}} | sed "s/\\.spot\\.out//"`
        cat ${{i}} | sed '1d' | awk -v sample="${{sample}}" '{{print sample"\\t"$NF}}' >> {output.spots}
        done
        sed -i '1i Sample\\tSPOT' {output.spots}
        Rscript workflow/scripts/calculate_SPOT_score_plot.R {output.spots} {input.samples} {output.spots_plot}
        '''

# ---------------- MultiQC report ----------------- #
rule ATAC_34_multiQC_inputs:
    input:
        expand("results/{prefix}/Step2_bowtie2_align/{sample}_bowtie2_align.log", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step2_bowtie2_align/{sample}_rm_dup.log", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step5_QCs/FastQC/{sample}.{group}_fastqc.zip", prefix= config["prefix"], sample = ALL_SAMPLES, group=[1, 2]),
        expand("results/{prefix}/Step5_QCs/insert_size/{sample}.isize.txt", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.{peak_shapes}Peak", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step4_call_peaks/MACS2_Sample_{peak_shapes}Peak/{sample}_macs2_peaks.xls", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{peak_shapes}Peak.count_mqc.tsv", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.tsv", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step5_QCs/RepeatabilityCorrelate/{condition}_Correlate.txt", prefix= config["prefix"], condition=ALL_CONDITIONS),
        expand("results/{prefix}/Step5_QCs/fingerPrint/{sample}.qualityMetrics.tsv", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step5_QCs/fingerPrint/{sample}.rawcounts.tsv", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step5_QCs/{prefix}_All_Condition_macs_{peak_shapes}peak.summary.txt", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"]) if config["params"]["Biological_Replicates"] else [],
        expand("results/{prefix}/Step5_QCs/{prefix}_All_Sample_macs_{peak_shapes}peak.summary.txt", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step5_QCs/Signal_to_Noise_SPOT/{sample}/{sample}.spot.out", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotProfile.tab", prefix= config["prefix"], sample = ALL_SAMPLES)
    output: 
        files = "results/{}/Step5_QCs/multiQC/{}_multiQC_inputs.txt".format(config["prefix"], config["prefix"])
    run:
        with open(output.files, 'w') as outfile:
            for fname in input:
                outfile.write(fname + "\n")

rule ATAC_35_multiQC:
    input:
        "results/{}/Step5_QCs/multiQC/{}_multiQC_inputs.txt".format(config["prefix"], config["prefix"])
    output: 
        report("results/{}/Step5_QCs/multiQC/{}_multiQC_report.html".format(config["prefix"], config["prefix"]), caption="../report/multiQC.rst", category="1. 测序数据质控")
    params:
        name = lambda w, output: "".join(os.path.basename(output[0]).split(".")[0]),
        outdir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/ATAC_seq_multiqc.yaml"
    shell:
        """
        multiqc -o {params.outdir} -l {input} -f -v -n {params.name}
        """

rule ATAC_Step3_Repetition_Peak_Quality_list:
    input:
        expand("results/{prefix}/Step5_QCs/Distribution_Mapped_Reads/{sample}_Distribution_Mapped_Reads.pdf", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step5_QCs/FRiP_score/{sample}_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.tsv", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        expand("results/{prefix}/Step5_QCs/FRiP_score/{prefix}_All_sample_macs2_peaks.{peak_shapes}Peak.FRiP_mqc.pdf", prefix= config["prefix"], sample = ALL_SAMPLES, peak_shapes = config["params"]["macs2_peak_shapes"]),
        "results/{}/Step5_QCs/RepeatabilityCorrelate/{}_All_Sample_Correlate.txt".format(config["prefix"], config["prefix"]),
        "results/{}/Step5_QCs/RepeatabilityTest/{}_PCA_readCounts.pdf".format(config["prefix"], config["prefix"]),
        "results/{}/Step5_QCs/RepeatabilityTest/{}_Correlation_PCA.pdf".format(config["prefix"], config["prefix"]),
        expand("results/{prefix}/Step5_QCs/ataqv/{prefix}_All_Sample_ATAC_{peak_shapes}peaks_QC_html/index.html", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"]),
        "results/{}/Step5_QCs/multiQC/{}_multiQC_report.html".format(config["prefix"], config["prefix"]),
        "results/{}/Step5_QCs/Signal_to_Noise_SPOT/{}_All_Sample_SPOT.pdf".format(config["prefix"], config["prefix"]),
        expand("results/{prefix}/Step5_QCs/macs_peak_intersect/{prefix}_All_Condition_macs_{peak_shapes}peak_intersect.pdf", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"]), # if SAMPLESSIZE > 1 else [],
        expand("results/{prefix}/Step5_QCs/macs_peak_intersect/{prefix}_All_Samples_macs_{peak_shapes}peak_intersect.pdf", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"])

###################################################################################################################################
# ----- 5. Visualization ----- #
###################################################################################################################################
rule ATAC_39_MAKE_GENE_TSS_TTS_BED:
    input:
        genome_gtf = config["DATA"]["genome_gtf"]
    output:
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"]),
        tss_bed = "results/{}/Input_Genome_Data/{}_TSS.bed".format(config["prefix"], config["DATA"]["genome_name"]),
        tts_bed = "results/{}/Input_Genome_Data/{}_TTS.bed".format(config["prefix"], config["DATA"]["genome_name"])
    shell:
        """
        perl workflow/scripts/gtf2bed.pl {input.genome_gtf} > {output.gene_bed}
        cat {output.gene_bed} | awk -v FS='\\t' -v OFS='\\t' '{{ if($6=="+") $3=$2+1; else $2=$3-1; print $1, $2, $3, $4, $5, $6;}}' > {output.tss_bed}
        cat {output.gene_bed} | awk -v FS='\\t' -v OFS='\\t' '{{ if($6=="+") {{print $1, $3, $3+1, $4, $5, $6;}} else {{print $1, $2-1, $2, $4, $5, $6;}}}}' > {output.tts_bed}
        """

rule ATAC_40_computeMatrix_body:
    input:
        bw ="results/{prefix}/Step3_Sample_bigwig/{sample}.bw",
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_computeMatrix.mat.gz",
        mat_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_computeMatrix.vals.mat.tab"
    threads: 5
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        computeMatrix scale-regions \
            --regionsFileName {input.gene_bed} \
            --scoreFileName {input.bw} \
            --outFileName {output.mat_gz} \
            --outFileNameMatrix {output.mat_tab} \
            --regionBodyLength 1000 \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --smartLabels \
            --numberOfProcessors {threads}
        """

rule ATAC_41_plotProfile_plotHeatmap_body:
    input:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_computeMatrix.mat.gz"
    output:
        plotHeatmap = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotHeatmap.pdf",
        plotHeatmap_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotHeatmap.tab",
        plotProfile = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotProfile.pdf",
        plotProfile_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotProfile.tab"
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile} --outFileNameData {output.plotProfile_tab} --plotHeight 8 --plotWidth 12
        plotHeatmap --matrixFile {input.mat_gz} --outFileName {output.plotHeatmap} --outFileNameMatrix {output.plotHeatmap_tab}
        """

rule ATAC_42_computeMatrix_body_whole:
    input:
        bw =expand("results/{prefix}/Step3_Sample_bigwig/{sample}.bw", prefix= config["prefix"], sample = ALL_SAMPLES),
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mat_gz = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_computeMatrix.mat.gz".format(config["prefix"], config["prefix"]),
        mat_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_computeMatrix.vals.mat.tab".format(config["prefix"], config["prefix"])
    threads: 10
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        computeMatrix scale-regions \
            --regionsFileName {input.gene_bed} \
            --scoreFileName {input.bw} \
            --outFileName {output.mat_gz} \
            --outFileNameMatrix {output.mat_tab} \
            --regionBodyLength 1000 \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --smartLabels \
            --numberOfProcessors {threads}
        """

rule ATAC_43_plotProfile_plotHeatmap_body_whole:
    input:
        mat_gz = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_computeMatrix.mat.gz".format(config["prefix"], config["prefix"])
    output:
        plotHeatmap = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotHeatmap.pdf".format(config["prefix"], config["prefix"]),
        plotHeatmap_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotHeatmap.tab".format(config["prefix"], config["prefix"]),
        plotProfile = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile.pdf".format(config["prefix"], config["prefix"]),
        plotProfile_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile.tab".format(config["prefix"], config["prefix"]),
        plotProfile2 = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile2.pdf".format(config["prefix"], config["prefix"]),
        plotProfile_tab2 = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile2.tab".format(config["prefix"], config["prefix"])
    params:
        deeptools = config["software"]["deeptools"]
    threads: 10
    shell:
        """
        {params.deeptools}
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile} --outFileNameData {output.plotProfile_tab} --plotHeight 12 --plotWidth 12 --numPlotsPerRow 5
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile2} --outFileNameData {output.plotProfile_tab2} --plotHeight 12 --plotWidth 12 --perGroup
        plotHeatmap --matrixFile {input.mat_gz} --outFileName {output.plotHeatmap} --outFileNameMatrix {output.plotHeatmap_tab}
        """

rule ATAC_44_plotProfile_body_whole_Rplot:
    input:
        plotProfile = expand("results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotProfile.tab", prefix= config["prefix"], sample = ALL_SAMPLES),
        samples = config["samples"]
    output:
        plot_data = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile_Rplot.txt".format(config["prefix"], config["prefix"]),
        plot = "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile_Rplot.pdf".format(config["prefix"], config["prefix"])
    params:
        R = config["software"]["R"],
        gcc = config["software"]["gcc"]
    shell:
        """
        {params.R}
        {params.gcc}
        cat {input.plotProfile} | sort | uniq | \\grep -E "^bin labels" > {output.plot_data}2
        cat {input.plotProfile} | sort | uniq | \\grep -v -E "^bin" > {output.plot_data}3
        cat {output.plot_data}2 {output.plot_data}3 > {output.plot_data}
        rm -rf {output.plot_data}2 {output.plot_data}3
        Rscript workflow/scripts/plotProfile_body_plot.R {input.samples} {output.plot_data} {output.plot}
        """

rule ATAC_45_computeMatrix_TSS:
    input:
        bw ="results/{prefix}/Step3_Sample_bigwig/{sample}.bw",
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_computeMatrix.mat.gz",
        mat_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_computeMatrix.vals.mat.tab"
    threads: 5
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        computeMatrix reference-point \
            --referencePoint TSS \
            --regionsFileName {input.gene_bed} \
            --scoreFileName {input.bw} \
            --outFileName {output.mat_gz} \
            --outFileNameMatrix {output.mat_tab} \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --smartLabels \
            --numberOfProcessors {threads}
        """

rule ATAC_46_plotProfile_plotHeatmap_TSS:
    input:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_computeMatrix.mat.gz"
    output:
        plotHeatmap = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_plotHeatmap.pdf",
        plotHeatmap_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_plotHeatmap.tab",
        plotProfile = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_plotProfile.pdf",
        plotProfile_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_plotProfile.tab"
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile} --outFileNameData {output.plotProfile_tab}
        plotHeatmap --matrixFile {input.mat_gz} --outFileName {output.plotHeatmap} --outFileNameMatrix {output.plotHeatmap_tab}
        """

rule ATAC_45_computeMatrix_TSS_whole:
    input:
        bw = expand("results/{prefix}/Step3_Sample_bigwig/{sample}.bw", prefix= config["prefix"], sample = ALL_SAMPLES),
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mat_gz = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_computeMatrix.mat.gz".format(config["prefix"], config["prefix"]),
        mat_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_computeMatrix.vals.mat.tab".format(config["prefix"], config["prefix"])
    threads: 5
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        computeMatrix reference-point \
            --referencePoint TSS \
            --regionsFileName {input.gene_bed} \
            --scoreFileName {input.bw} \
            --outFileName {output.mat_gz} \
            --outFileNameMatrix {output.mat_tab} \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --smartLabels \
            --numberOfProcessors {threads}
        """

rule ATAC_46_plotProfile_plotHeatmap_TSS_whole:
    input:
        mat_gz = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_computeMatrix.mat.gz".format(config["prefix"], config["prefix"])
    output:
        plotHeatmap = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_plotHeatmap.pdf".format(config["prefix"], config["prefix"]),
        plotHeatmap_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_plotHeatmap.tab".format(config["prefix"], config["prefix"]),
        plotProfile = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_plotProfile.pdf".format(config["prefix"], config["prefix"]),
        plotProfile_tab = "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_plotProfile.tab".format(config["prefix"], config["prefix"])
    params:
        deeptools = config["software"]["deeptools"]
    threads: 5
    shell:
        """
        {params.deeptools}
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile} --outFileNameData {output.plotProfile_tab} --plotHeight 12 --plotWidth 12 --numPlotsPerRow 6
        plotHeatmap --matrixFile {input.mat_gz} --outFileName {output.plotHeatmap} --outFileNameMatrix {output.plotHeatmap_tab}
        """

rule ATAC_47_computeMatrix_TES:
    input:
        bw ="results/{prefix}/Step3_Sample_bigwig/{sample}.bw",
        gene_bed = "results/{}/Input_Genome_Data/{}_gene.bed".format(config["prefix"], config["DATA"]["genome_name"])
    output:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_computeMatrix.mat.gz",
        mat_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_computeMatrix.vals.mat.tab"
    threads: 5
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        computeMatrix reference-point \
            --referencePoint TES \
            --regionsFileName {input.gene_bed} \
            --scoreFileName {input.bw} \
            --outFileName {output.mat_gz} \
            --outFileNameMatrix {output.mat_tab} \
            --beforeRegionStartLength 3000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --smartLabels \
            --numberOfProcessors {threads}
        """

rule ATAC_48_plotProfile_plotHeatmap_TES:
    input:
        mat_gz = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_computeMatrix.mat.gz"
    output:
        plotHeatmap = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_plotHeatmap.pdf",
        plotHeatmap_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_plotHeatmap.tab",
        plotProfile = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_plotProfile.pdf",
        plotProfile_tab = "results/{prefix}/Step3_Sample_bigwig/{sample}_TES_plotProfile.tab"
    params:
        deeptools = config["software"]["deeptools"]
    shell:
        """
        {params.deeptools}
        plotProfile --matrixFile {input.mat_gz} --outFileName {output.plotProfile} --outFileNameData {output.plotProfile_tab}
        plotHeatmap --matrixFile {input.mat_gz} --outFileName {output.plotHeatmap} --outFileNameMatrix {output.plotHeatmap_tab}
        """

rule ATAC_Step4_TSSlist:
    input:
        expand("results/{prefix}/Step3_Sample_bigwig/{sample}_body_plotProfile.pdf", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step3_Sample_bigwig/{sample}_TSS_plotProfile.pdf", prefix= config["prefix"], sample = ALL_SAMPLES),
        expand("results/{prefix}/Step3_Sample_bigwig/{sample}_TES_plotProfile.pdf", prefix= config["prefix"], sample = ALL_SAMPLES),
        "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile.pdf".format(config["prefix"], config["prefix"]),
        "results/{}/Step3_Sample_bigwig/{}_All_Samples_body_plotProfile_Rplot.pdf".format(config["prefix"], config["prefix"]),
        "results/{}/Step3_Sample_bigwig/{}_All_Samples_TSS_plotHeatmap.pdf".format(config["prefix"], config["prefix"])

###################################################################################################################################
# ----- 10. Comparing peak files ----- #
###################################################################################################################################
rule ATAC_73_Diff_Peaks_DiffBind_whole:
    input:
        bams = expand("results/{prefix}/Step2_bowtie2_align/{sample}.bam", prefix= config["prefix"], sample = ALL_SAMPLES),
        narrowPeaks = expand("results/{prefix}/Step4_call_peaks/MACS2_Sample_{{peak_shapes}}Peak/{sample}_macs2_peaks.{{peak_shapes}}Peak", prefix= config["prefix"], sample = ALL_SAMPLES)
    output:
        metadata = "results/{prefix}/Step10_Comparing_peak/{prefix}_whole_DiffBind_metadata_{peak_shapes}Peak.txt",
        DiffBind = directory("results/{prefix}/Step10_Comparing_peak/DiffBind_{peak_shapes}/{prefix}_whole")
    params:
        workdir = config["workdir"],
        R42 = config["software"]["R42"],
        gcc = config["software"]["gcc"],
        prefix = config["prefix"],
        peak_shapes = "{peak_shapes}",
        bamdir = lambda w, input: os.path.dirname(input['bams'][0]),
        peakdir = lambda w, input: os.path.dirname(input['narrowPeaks'][0])
    threads: 20
    shell:
        """
        {params.R42}
        {params.gcc}
        mkdir -p {output.DiffBind}
        echo -e "SampleID\\tTissue\\tFactor\\tCondition\\tTreatment\\tReplicate\\tbamReads\\tPeaks\\tPeakCaller" >> {output.metadata}
        for i in {input.bams}
        do
        sample=$(basename ${{i}} ".bam")
        condition=`echo ${{sample}} | cut -f 1 -d -`
        replicate=`echo ${{sample}} | cut -f 2 -d -`
        echo -e "${{sample}}\\tHypocotyl\\t{params.prefix}\\t${{condition}}\\tTime\\t${{replicate}}\\t{params.workdir}/${{i}}\\t{params.workdir}/{params.peakdir}/${{sample}}_macs2_peaks.{params.peak_shapes}Peak\\t{params.peak_shapes}Peak" >> {output.metadata}
        done
        Rscript workflow/scripts/DiffBind.R {output.metadata} {output.DiffBind} {params.prefix}
        """

rule ATAC_Diff_Peaks_DiffBind_list:
    input:
        expand("results/{prefix}/Step10_Comparing_peak/DiffBind_{peak_shapes}/{prefix}_whole", prefix= config["prefix"], peak_shapes = config["params"]["macs2_peak_shapes"])
