def select_input_genome():
    genome = "results/{}/Plus_finally_genome/Chr_{}_genome.fasta".format(config["prefix"], config["prefix"])
    return genome

def input_genome_name():
    input_file = select_input_genome()
    name = os.path.basename(input_file)
    return name, name.split(".")[0]

localrules:
    HiCHP_3_genome_HiCPro_Step_to_Step_run_command,
    HiCHP_4_genome_HiCPro_Step_to_Step_run_step1,
    HiCHP_5_genome_HiCPro_Step_to_Step_run_step2

#################################################################################################################################################
rule HiCHP_1_genome_HiCPro_prepare:
    input:
        genome = select_input_genome(),
        HiCPro = "logs/global/HiCPro_install_help.log"
    output:
        bed = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/annotation/{}_genome_{}.bed".format(config["prefix"], config["prefix"], config["HiC"]["enzyme"].split(",")[1]),
        chromsize = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/annotation/{}_genome.chrom.size".format(config["prefix"], config["prefix"]),
        bowtie2index = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/annotation/{}.1.bt2l".format(config["prefix"], input_genome_name()[1]) if config["bigGenome"] else "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/annotation/{}.1.bt2".format(config["prefix"], input_genome_name()[1]),
        config = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/{}_config-hicpro.txt".format(config["prefix"], config["prefix"])
    threads: 10
    conda:
        "../bin/HiC-Pro-3.1.0/environment.yml"
    params:
        workdir = config["workdir"],
        enzyme = config["HiC"]["enzyme"],
        HiCPro_enzyme = config["HiC"]["HiCPro_enzyme"],
        HiCPro_BIN_SIZE = config["HiC"]["HiCPro_BIN_SIZE"],
        prefix = lambda w, input: os.path.basename(input[0]),
        outdir = lambda w, output: os.path.dirname(output["bed"]),
        Bowtie2 = config["software"]["Bowtie2"]
    log:
        "logs/{}/Assessment_Finally_Assembly/HiC_Heatmap/assembly_HiCPro_prepare.log".format(config["prefix"])
    shell:
        '''
        #{params.Bowtie2}
        ## Annotation Files
        ## BED
        cd {params.workdir}
        enzyme_seq=`echo {params.enzyme} | cut -f 1 -d ,`
        {params.workdir}/workflow/bin/HiC-Pro-3.1.0/bin/utils/digest_genome.py -r ${{enzyme_seq}} -o {output.bed} {input.genome} >> {params.workdir}/{log} 2>&1
        ## chromosome size
        seqkit fx2tab -nl {input.genome} | awk '{{print $1"\\t"$2}}' > {output.chromsize}
        # bowtie2 index
        cd {params.workdir}/{params.outdir}
        genome_name="{params.prefix}"
        bowtie2-build --threads {threads} {params.workdir}/{input.genome} ${{genome_name%.*}} >> {params.workdir}/{log} 2>&1

        ## config-hicpro
        cd {params.workdir}
        cat {params.workdir}/workflow/bin/HiC-Pro-3.1.0/config-hicpro.txt | \
          sed "s|N_CPU = .*|N_CPU = 20|g" | \
          sed "s|BOWTIE2_IDX_PATH = .*|BOWTIE2_IDX_PATH = {params.workdir}/{params.outdir}|g" | \
          sed "s|REFERENCE_GENOME = .*|REFERENCE_GENOME = ${{genome_name%.*}}|g" | \
          sed "s|GENOME_SIZE = .*|GENOME_SIZE = {params.workdir}/{output.chromsize}|g" | \
          sed "s|GENOME_FRAGMENT = .*|GENOME_FRAGMENT = {params.workdir}/{output.bed}|g" | \
          sed "s|LIGATION_SITE = .*|LIGATION_SITE = {params.HiCPro_enzyme}|g" | \
          sed "s/BIN_SIZE = .*/BIN_SIZE = {params.HiCPro_BIN_SIZE}/g" > {output.config}
        '''

# rule HiCHP_2_genome_HiCPro_One_Step_run:
#     input:
#         config = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/{}_config-hicpro.txt".format(config["prefix"], config["prefix"]),
#         r1 = expand("results/{prefix}/HiC/Trim_galore/{u.sample}-{u.unit}/{u.sample}-{u.unit}_R1.fastq.gz", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples()),
#         r2 = expand("results/{prefix}/HiC/Trim_galore/{u.sample}-{u.unit}/{u.sample}-{u.unit}_R2.fastq.gz", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples())
#     output:
#         validPairs = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/hic_results/data/{}/{}.allValidPairs".format(config["prefix"], config["prefix"], config["prefix"])
#     threads: 1
#     conda:
#         "../bin/HiC-Pro-3.1.0/environment.yml"
#     params:
#         workdir = config["workdir"],
#         prefix = config["prefix"],
#         Bowtie2 = config["software"]["Bowtie2"],
#         indir = lambda w, input: os.path.dirname(input["config"]),
#         outdir = lambda w, output: "/".join(os.path.dirname(output[0]).split("/")[0:6])
#     log:
#         "logs/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC_Heatmap_assembly_HiCPro_run.log".format(config["prefix"])
#     shell:
#         '''
#         #{params.Bowtie2}
#         rm -rf {params.outdir}
#         mkdir -p {params.indir}/HiC_Trim_fqs/{params.prefix}
#         for i in {input.r1} {input.r2}
#         do
#         ln -s {params.workdir}/${{i}} {params.workdir}/{params.indir}/HiC_Trim_fqs/{params.prefix}
#         done

#         {params.workdir}/workflow/bin/HiC-Pro-3.1.0/bin/HiC-Pro \
#           -i {params.workdir}/{params.indir}/HiC_Trim_fqs/ \
#           -o {params.outdir} -c {input.config} > {log} 2>&1
#         '''

rule HiCHP_3_genome_HiCPro_Step_to_Step_run_command:
    input:
        config = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/{}_config-hicpro.txt".format(config["prefix"], config["prefix"]),
        r1 = expand("results/{prefix}/HiC/Trim_galore/{u.sample}-{u.unit}/{u.sample}-{u.unit}_R1.fastq.gz", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples()),
        r2 = expand("results/{prefix}/HiC/Trim_galore/{u.sample}-{u.unit}/{u.sample}-{u.unit}_R2.fastq.gz", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples())
    output:
        step1_sh = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/HiCPro_step1_.sh".format(config["prefix"]),
        step2_sh = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/HiCPro_step2_.sh".format(config["prefix"])
    threads: 20
    conda:
        "../bin/HiC-Pro-3.1.0/environment.yml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        Bowtie2 = config["software"]["Bowtie2"],
        indir = lambda w, input: os.path.dirname(input["config"]),
        outdir = lambda w, output: "/".join(os.path.dirname(output[0]).split("/")[0:6])
    shell:
        '''
        #{params.Bowtie2}
        rm -rf {params.outdir}
        mkdir -p {params.indir}/HiC_Trim_fqs/{params.prefix}
        for i in {input.r1} {input.r2}
        do
        ln -s {params.workdir}/${{i}} {params.workdir}/{params.indir}/HiC_Trim_fqs/{params.prefix}
        done

        {params.workdir}/workflow/bin/HiC-Pro-3.1.0/bin/HiC-Pro \
          -i {params.workdir}/{params.indir}/HiC_Trim_fqs/ \
          -o {params.workdir}/{params.outdir} -c {params.workdir}/{input.config} -p {threads}
        '''

rule HiCHP_4_genome_HiCPro_Step_to_Step_run_step1:
    input:
        "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/HiCPro_step1_.sh".format(config["prefix"])
    output:
        step1_command = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/My_HiCPro_step1_.sh".format(config["prefix"]),
        step1_results_validPairs = expand("results/{prefix}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/hic_results/data/{prefix}/{u.sample}-{u.unit}_{genomeName}.bwt2pairs.validPairs", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples(), genomeName = input_genome_name()[1]),
        step1_results_bam = expand("results/{prefix}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/bowtie_results/bwt2/{prefix}/{u.sample}-{u.unit}_{genomeName}.bwt2pairs.bam", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples(), genomeName = input_genome_name()[1])
    threads: 1
    conda:
        "../bin/HiC-Pro-3.1.0/environment.yml"
    params:
        workdir = config["workdir"],
        indir = lambda w, input: os.path.dirname(input[0])
    shell:
        '''
        cd {params.indir}
        cat {params.workdir}/{input} | grep -v -E "#BSUB -M|#BSUB -W|#BSUB -N" | sed "s/#BSUB -u.*/#BSUB -R span[hosts=1]/" | sed "s/#BSUB -q.*/#BSUB -q normal/" > {params.workdir}/{output.step1_command}
        bsub < {params.workdir}/{output.step1_command}
        '''

rule HiCHP_5_genome_HiCPro_Step_to_Step_run_step2:
    input:
        step1_results = expand("results/{prefix}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/hic_results/data/{prefix}/{u.sample}-{u.unit}_{genomeName}.bwt2pairs.validPairs", prefix = config["prefix"], u=units[units['platform'].str.match('HiC')].itertuples(), genomeName = input_genome_name()[1]),
        step2_sh = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/HiCPro_step2_.sh".format(config["prefix"])
    output:
        step2_sh = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/My_HiCPro_step2_.sh".format(config["prefix"]),
        HiCPro = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/hic_results/data/{}/{}.allValidPairs".format(config["prefix"], config["prefix"], config["prefix"])
    threads: 1
    conda:
        "../bin/HiC-Pro-3.1.0/environment.yml"
    params:
        workdir = config["workdir"],
        indir = lambda w, input: os.path.dirname(input['step2_sh'])
    shell:
        '''
        cd {params.indir}
        cat {params.workdir}/{input.step2_sh} | grep -v -E "#BSUB -M|#BSUB -W|#BSUB -N" | sed "s/#BSUB -u.*/#BSUB -R span[hosts=1]/" | sed "s/#BSUB -q.*/#BSUB -q smp/" > {params.workdir}/{output.step2_sh}
        bsub < {params.workdir}/{output.step2_sh}
        '''

rule HiCHP_6_genome_HiCPro_Plot:
    input:
        config = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/{}_config-hicpro.txt".format(config["prefix"], config["prefix"]),
        HiCPro = "results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/results/hic_results/data/{}/{}.allValidPairs".format(config["prefix"], config["prefix"], config["prefix"])
    output:
        directory("results/{}/Assessment_Finally_Assembly/HiC_Heatmap/HiC-Pro/Interation_Heatmap_Plot".format(config["prefix"]))
    threads: 10
    conda:
        "../bin/HiC-Pro-3.1.0/environment.yml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        chromosome = config["params"]["chromosome"][:-1],
        indir = lambda w,input: "/".join(os.path.dirname(input["HiCPro"]).split("/")[0:7])
    shell:
        '''
        module load HiCExplorer/3.7.2
        export NUMEXPR_MAX_THREADS={threads}
        mkdir -p {output}
        for i in `cat {input.config} | \\grep "BIN_SIZE" | cut -f 2 -d = | sed "s/^ //" | sed "s/ /\\n/g" | sed '1,2d' | tac | xargs`
        do
        hicConvertFormat \
            -m {params.indir}/matrix/{params.prefix}/iced/${{i}}/{params.prefix}_${{i}}_iced.matrix \
            -o {output}/{params.prefix}_${{i}}_iced.matrix \
            --bedFileHicpro {params.indir}/matrix/{params.prefix}/raw/${{i}}/{params.prefix}_${{i}}_abs.bed \
            --inputFormat hicpro --outputFormat h5
        # Whole chromosomes
        hicPlotMatrix --matrix {output}/{params.prefix}_${{i}}_iced.matrix.h5 -out {output}/{params.prefix}_${{i}}_iced.matrix.pdf --chromosomeOrder {params.chromosome} --log1p --colorMap gist_heat_r --clearMaskedBins

        # Single chromosome
        for j in {params.chromosome}
        do
        hicPlotMatrix --matrix {output}/{params.prefix}_${{i}}_iced.matrix.h5 -out {output}/{params.prefix}_${{j}}_${{i}}_iced.matrix.pdf --region ${{j}} --log1p --colorMap gist_heat_r --clearMaskedBins
        done
        done
        '''
