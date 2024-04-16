###################################################################################################################################
###################################################################################################################################
rule MGCS_1_SyRI_split_genome:
    input:
        expand("resources/{prefix}/{MGCSs}/Chr_{MGCSs}.fasta", prefix = config["prefix"], MGCSs = list(np.unique(np.array(Multiple_Genomic_Chromosome_SyRI_species))))
    output:
        temp("results/{}/Multiple_Genomic_Chromosome_SyRI/All_{}_Genome.fasta".format(config["prefix"], config["prefix"]))
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Structural_Variation.yaml"
    shell:
        '''
        cat {input} > {output}
        '''

rule MGCS_2_SyRI_homo_genome:
    input:
        "results/{}/Multiple_Genomic_Chromosome_SyRI/All_{}_Genome.fasta".format(config["prefix"], config["prefix"])
    output:
        expand("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchr}.fa", prefix = config["prefix"], MGCSchr = list(np.unique(np.array(Multiple_Genomic_Chromosome_SyRI))))
    threads: 1
    params:
        workdir = config["workdir"],
        ucsc = config["software"]["ucsc"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        {params.ucsc}
        faSplit byname {params.workdir}/{input} {params.workdir}/{params.outdir}/
        '''

rule MGCS_3_SyRI_minimap2_split:
    input:
        pchrA = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrA}.fa",
        pchrB = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrB}.fa",
        gnome_format = config["Multiple_Genomic_Chromosome_SyRI"]
    output:
        "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.sam"
    threads: 5
    params:
        workdir = config["workdir"],
        chrA = "{MGCSchrA}",
        chrAdir = lambda w, input: os.path.dirname(input['pchrA']),
        chrAname = lambda w, input: os.path.basename(input['pchrA']),
        chrB = "{MGCSchrB}",
        chrBdir = lambda w, input: os.path.dirname(input['pchrB']),
        chrBname = lambda w, input: os.path.basename(input['pchrB'])
    conda: "../envs/Structural_Variation.yaml"
    shell:
        '''
        set +eu
        source ~/.bashrc
        source activate
        conda activate seqkit
        cat {input.gnome_format} | \\grep -P "{params.chrA}\\t.*{params.chrB}$" | \\grep "-"
        if [ "$?" -eq 0 ]
        then
            cat {input.gnome_format} | \\grep -P "{params.chrA}\\t" | cut -f 1 | \\grep "-"
            if [ "$?" -eq 0 ]
            then
                # chrA
                seqkit seq -t DNA -r -p {input.pchrA} > {params.chrAdir}/RC_{params.chrA}-{params.chrB}_{params.chrAname}
                minimap2 -ax asm5 -t {threads} --eqx {params.chrAdir}/RC_{params.chrA}-{params.chrB}_{params.chrAname} {input.pchrB} > {output}
            else
                echo "Inverse and complementary of {input.pchrB}"
            fi

            cat {input.gnome_format} | \\grep "{params.chrB}$" | cut -f 2 | \\grep "-"
            if [ "$?" -eq 0 ]
            then
                # chrB
                seqkit seq -t DNA -r -p {input.pchrB} > {params.chrBdir}/RC_{params.chrA}-{params.chrB}_{params.chrBname}
                minimap2 -ax asm5 -t {threads} --eqx {input.pchrA} {params.chrBdir}/RC_{params.chrA}-{params.chrB}_{params.chrBname} > {output}
            else
                echo "Inverse and complementary of {input.pchrA}"
            fi
        else
            minimap2 -ax asm5 -t {threads} --eqx {input.pchrA} {input.pchrB} > {output}
        fi
        '''

rule MGCS_4_SyRI_minimap2_split_paf_plot:
    input:
        "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.sam"
    output:
        paf = temp("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.paf"),
        pdf = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/paf_plot/{MGCSchrA}-{MGCSchrB}.png"
    threads: 1
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Structural_Variation.yaml"
    shell:
        '''
        paftools.js sam2paf {input} > {output.paf}
        Rscript workflow/scripts/paf_plot.R {output.paf} {output.pdf}
        '''

rule MGCS_5_SyRI_minimap2_sam2bam_split:
    input:
        "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.sam"
    output:
        "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.bam"
    threads: 2
    params:
        workdir = config["workdir"]
    conda:
        "../envs/Structural_Variation.yaml"
    shell:
        '''
        samtools sort -O BAM {input} > {output}
        samtools index {output}
        '''

rule MGCS_6_SyRI_Annotation_minimap2_split:
    input:
        pchrA = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrA}.fa",
        pchrB = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrB}.fa",
        bam = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}.sam"
    output:
        syri_minimap2_summary = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}_minimap2_syri.summary",
        syri_minimap2_out = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}_minimap2_syri.out",
        syri_minimap2_vcf = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}_minimap2_syri.vcf"
    params:
        workdir = config["workdir"],
        chr = "{MGCSchrA}-{MGCSchrB}",
        outpre = lambda w, output: "_".join(os.path.basename(output[0]).split("_")[:-1]),
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 5
    conda:
        "../envs/Structural_Variation_SyRI.yaml"
    shell:
        '''
        python3 {params.workdir}/workflow/bin/syri-v1.6/bin/syri \
          -c {params.workdir}/{input.bam} \
          -r {params.workdir}/{input.pchrA} -q {params.workdir}/{input.pchrB} \
          --dir {params.outdir} \
          --prefix {params.chr}_minimap2_ --nc {threads} -F S
        '''

rule MGCS_7_SyRI_Run_Plot_pairs_prepare_genome_split:
    input:
        pchrA = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrA}.fa",
        pchrB = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/Homologous_Chromosome/{MGCSchrB}.fa"
    output:
        chrs = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{MGCSchrA}-{MGCSchrB}_chrs.txt"
    params:
        workdir = config["workdir"],
        chrA = "{MGCSchrA}",
        chrB = "{MGCSchrB}",
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 1
    shell:
        '''
        echo -e "{params.workdir}/{input.pchrA}\\t{params.chrA}\\tlw:10" >> {output.chrs}
        echo -e "{params.workdir}/{input.pchrB}\\t{params.chrB}\\tlw:10" >> {output.chrs}
        '''

rule MGCS_8_SyRI_Run_Plot_minimap2_split:
    input:
        syri_minimap2 = expand("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{chr.chrA}-{chr.chrB}_minimap2_syri.out", prefix = config["prefix"], chr = Multiple_Genomic_Chromosome_SyRI.itertuples()),
        chrs = expand("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{chr.chrA}-{chr.chrB}_chrs.txt", prefix = config["prefix"], chr = Multiple_Genomic_Chromosome_SyRI.itertuples()),
        genome_chr = "config/{prefix}_Multiple_Genomic_Chromosome_SyRI_Homologous_Chromosome_Group.txt"
    output:
        chrs = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{commonChromosome}_chrs.txt",
        syri_minimap2 = "results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{commonChromosome}_minimap2_syri.pdf"
    params:
        workdir = config["workdir"],
        commonChromosome = "{commonChromosome}",
        syridir = lambda w, input: os.path.dirname(input['syri_minimap2'][0]),
        chrsdir = lambda w, input: os.path.dirname(input['chrs'][0]),
        outdir = lambda w, output: os.path.dirname(output[0])
    threads: 1
    conda:
        "../envs/Structural_Variation_SyRI.yaml"
    shell:
        '''
        chrs=`cat {input.genome_chr} | \\grep -w "^{params.commonChromosome}" | sed "s/{params.commonChromosome}\\t//" | sed "s/\\t/|/g"`
        sr=`echo {input.syri_minimap2} | sed "s/ /\\n/g" | sed "s/^/--sr /g" | \\grep -w -E "${{chrs}}"`
        cat {input.chrs} | uniq | \\grep -w -E "${{chrs}}" > {output.chrs}
        plotsr ${{sr}} \
            --genomes {output.chrs} \
            -o {output.syri_minimap2} \
            -S 0.5 -W 8 -H 8 -f 8 -b pdf
        '''

###################################################################################################################################################
############################################################### END ###############################################################################
###################################################################################################################################################
def common_chromosome_list():
    df = pd.read_table("config/{}_Multiple_Genomic_Chromosome_SyRI_Homologous_Chromosome_Group.txt".format(config["prefix"]), header = None)
    return df[0].tolist()

rule MGCS_9_list_SyRI:
    input:
        expand("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/paf_plot/{chr.chrA}-{chr.chrB}.png", prefix = config["prefix"], chr = Multiple_Genomic_Chromosome_SyRI.itertuples()),
        expand("results/{prefix}/Multiple_Genomic_Chromosome_SyRI/{commonChromosome}_minimap2_syri.pdf", prefix = config["prefix"], commonChromosome = common_chromosome_list())
