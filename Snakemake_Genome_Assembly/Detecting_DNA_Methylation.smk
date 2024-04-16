rule DDM_4_Primrose_bam:
    input:
        "results/{prefix}/RawData/HiFi_reads/{sample_pacbioraw}-{unit_pacbioraw}.ccs.bam"
    output:
        "results/{prefix}/Detecting_DNA_Methylation/Primrose/{sample_pacbioraw}-{unit_pacbioraw}.5mc.hifi_reads.bam"
    conda: 
        "../envs/Detecting_DNA_Methylation_Primrose.yaml"
    shell:
        '''
        primrose {input} {output}
        '''

rule DDM_5_hifi_bam_align:
    input:
        hifi_bam_5mc = expand("results/{prefix}/Detecting_DNA_Methylation/Primrose/{u.sample}-{u.unit}.5mc.hifi_reads.bam", prefix = config["prefix"], u=units[units['platform'].str.match('PACBIORAW')].itertuples()),
        genome = "results/{prefix}/Plus_finally_genome/{prefix}_genome.fasta"
    output:
        subreadset = "results/{prefix}/Detecting_DNA_Methylation/Primrose/{prefix}_subreadset.xml",
        align_bam = "results/{prefix}/Detecting_DNA_Methylation/Primrose/{prefix}.5mc.hifi_reads_sorted.bam"
    params:
        prefix = config["prefix"]
    conda: 
        "../envs/Detecting_DNA_Methylation_Primrose.yaml"
    threads: 15
    shell:
        '''
        module load SMRTLink/10.2.0.133434
        dataset create --type SubreadSet --name {params.prefix} {output.subreadset} {input.hifi_bam_5mc}
        pbmm2 align -j {threads} --preset CCS --sort --sample {params.prefix} {input.genome} {output.subreadset} {output.align_bam}
        samtools index {output.align_bam}
        '''

rule DDM_6_Primrose_bam_to_cpg:
    input:
        hifi_bam_5mc = "results/{}/Detecting_DNA_Methylation/Primrose/{}.5mc.hifi_reads_sorted.bam".format(config["prefix"], config["prefix"]),
        genome = "results/{}/Plus_finally_genome/{}_genome.fasta".format(config["prefix"], config["prefix"])
    output:
        "results/{}/Detecting_DNA_Methylation/Primrose/{}_CpG_5mc_hifi_reads.combined.reference.bed".format(config["prefix"], config["prefix"])
    conda: 
        "../envs/Detecting_DNA_Methylation_Primrose.yaml"
    params:
        workdir = config["workdir"],
        outprefix = lambda w, output: os.path.basename(output[0]).split(".")[0],
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 10
    shell:
        '''
        python workflow/bin/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b {input.hifi_bam_5mc} -f {input.genome} -o {params.outdir}/{params.outprefix} -t {threads} \
          --pileup_mode model --model_dir {params.workdir}/workflow/bin/pb-CpG-tools/pileup_calling_model --modsites reference
        '''

rule DDM_7_Primrose_5mc_density:
    input:
        cpg5mc = "results/{}/Detecting_DNA_Methylation/Primrose/{}_5mc_hifi_reads.combined.reference.bed".format(config["prefix"], config["prefix"]),
        genome = "results/{}/Plus_finally_genome/{}_genome.fasta".format(config["prefix"], config["prefix"]),
        genome_gene_gff = config["params"]["genome_gene_gff"],
        genome_repeat_gff = config["params"]["genome_repeat_gff"]
    output:
        chromsize = "results/{}/Detecting_DNA_Methylation/Primrose/{}.genome.size".format(config["prefix"], config["prefix"]),
        windows_bed = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        windows_GC = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.GC.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        windows_AT = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.AT.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        gene_bed = "results/{}/Detecting_DNA_Methylation/Primrose/{}.gene.bed".format(config["prefix"], config["prefix"]),
        gene_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.genes.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        repeat_bed = "results/{}/Detecting_DNA_Methylation/Primrose/{}.repeat.bed".format(config["prefix"], config["prefix"]),
        repeat_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.repeat.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        cpg5mc_bed = "results/{}/Detecting_DNA_Methylation/Primrose/{}.5mc.bed".format(config["prefix"], config["prefix"]),
        cpg5mc_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.5mc.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"])
    params:
        cpg5mc_density_windows = config["params"]["cpg5mc_density_windows"],
        prefix = config["prefix"]
    conda: 
        "../envs/Detecting_DNA_Methylation_Primrose.yaml"
    threads: 1
    shell:
        '''
        #######################################
        ## GC
        #######################################
        bioawk -c fastx '{{ print $name, length($seq) }}' < {input.genome} | grep -v "Scaffold" > {output.chromsize}
        /public/home/cotton/software/CscoreTool/get_windows.bed {output.chromsize} {output.windows_bed} {params.cpg5mc_density_windows}
        bedtools nuc -fi {input.genome} -bed {output.windows_bed} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$5}}' | grep -v "^#" > {output.windows_GC}
        bedtools nuc -fi {input.genome} -bed {output.windows_bed} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$4}}' | grep -v "^#" > {output.windows_AT}

        #######################################
        ## gene density
        #######################################
        \\grep -w "gene" {input.genome_gene_gff} | gff2bed > {output.gene_bed}

        for i in `cat {output.windows_bed} | cut -f 1 | sort -V | uniq | xargs`
        do
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.windows_bed} >> {params.prefix}_gdtmp1
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.gene_bed} >> {params.prefix}_gdtmp2
        done

        sort-bed {params.prefix}_gdtmp1 > {params.prefix}_gdtmp1.txt
        sort-bed {params.prefix}_gdtmp2 > {params.prefix}_gdtmp2.txt
        bedmap --echo --count {params.prefix}_gdtmp1.txt {params.prefix}_gdtmp2.txt > {params.prefix}_gdtmp3

        cat {params.prefix}_gdtmp3 | sort -k1,1V -k2,2n | sed "s/|/\\t/g" > {output.gene_density}
        \\rm {params.prefix}_gdtmp1 {params.prefix}_gdtmp2 {params.prefix}_gdtmp1.txt {params.prefix}_gdtmp2.txt {params.prefix}_gdtmp3

        #######################################
        ## repeat element density
        #######################################
        cat {input.genome_repeat_gff} | gff2bed > {output.repeat_bed}

        for i in `cat {output.windows_bed} | cut -f 1 | sort -V | uniq | xargs`
        do
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.windows_bed} >> {params.prefix}_redtmp1
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.repeat_bed} >> {params.prefix}_redtmp2
        done

        sort-bed {params.prefix}_redtmp1 > {params.prefix}_redtmp1.txt
        sort-bed {params.prefix}_redtmp2 > {params.prefix}_redtmp2.txt
        bedmap --echo --count {params.prefix}_redtmp1.txt {params.prefix}_redtmp2.txt > {params.prefix}_redtmp3

        cat {params.prefix}_redtmp3 | sort -k1,1V -k2,2n | sed "s/|/\\t/g" > {output.repeat_density}
        \\rm {params.prefix}_redtmp1 {params.prefix}_redtmp2 {params.prefix}_redtmp1.txt {params.prefix}_redtmp2.txt {params.prefix}_redtmp3

        #######################################
        ## 5mc_density
        #######################################
        cat {input.cpg5mc} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3}}' > {output.cpg5mc_bed}

        for i in `cat {output.windows_bed} | cut -f 1 | sort -V | uniq | xargs`
        do
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.windows_bed} >> {params.prefix}_5mcdtmp1
        awk -F"\\t" -v chr=${{i}} 'BEGIN{{OFS="\\t"}} $1==chr {{print $0}}' {output.cpg5mc_bed} >> {params.prefix}_5mcdtmp2
        done

        sort-bed {params.prefix}_5mcdtmp1 > {params.prefix}_5mcdtmp1.txt
        sort-bed {params.prefix}_5mcdtmp2 > {params.prefix}_5mcdtmp2.txt
        bedmap --echo --count {params.prefix}_5mcdtmp1.txt {params.prefix}_5mcdtmp2.txt > {params.prefix}_5mcdtmp3

        cat {params.prefix}_5mcdtmp3 | sort -k1,1V -k2,2n | sed "s/|/\\t/g" > {output.cpg5mc_density}
        \\rm {params.prefix}_5mcdtmp1 {params.prefix}_5mcdtmp2 {params.prefix}_5mcdtmp1.txt {params.prefix}_5mcdtmp2.txt {params.prefix}_5mcdtmp3
        '''

rule DDM_8_Primrose_5mc_density_plot:
    input:
        windows_GC = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.GC.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        windows_AT = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.AT.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        gene_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.genes.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        repeat_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.repeat.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        cpg5mc_density = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.5mc.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        cenh3_density = config["params"]["cenh3_density"]
    output:
        genome_feature = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.GC.genes.repeat.5mc.bed".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"]),
        genome_feature_plot = "results/{}/Detecting_DNA_Methylation/Primrose/{}.{}.GC.genes.repeat.5mc.pdf".format(config["prefix"], config["prefix"], config["params"]["cpg5mc_density_windows"])
    threads: 1
    shell:
        '''
        module load R/3.6.0
        echo -e "chr\\tstart\\tend\\tdensity\\ttype" > {output.genome_feature}
        for i in {input}
        do
        cat ${{i}} | awk -F"\\t" -v type=${{i##*/}} 'BEGIN{{OFS="\\t"}} {{print $0,type}}' >> {output.genome_feature}
        done
        Rscript workflow/scripts/CpG5mc_density_plot.R {output.genome_feature} {output.genome_feature_plot}
        '''

rule DDM_9_5mc_list:
    input:
        "results/{}/Detecting_DNA_Methylation/Primrose/{}.5mc.hifi_reads.bw".format(config["prefix"], config["prefix"])
