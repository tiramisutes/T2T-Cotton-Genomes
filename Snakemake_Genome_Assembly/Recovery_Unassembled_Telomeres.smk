def select_input_Chr_genome():
    return "results/{}/Plus_finally_genome/{}_genome.fasta".format(config["prefix"], config["prefix"])

def input_genome_name():
    input_file = select_input_Chr_genome()
    name = os.path.basename(input_file)
    return name

#######################################################################################################################################################################
########################################################### 1:  ########################################################################
#######################################################################################################################################################################
rule RUT_0_convert_hifi_reads2fastq:
    input:
        "/data/Cotton/Data/PacBio/m64446e_220710_022523.hifi_reads.bam"
    output:
        "/data/Cotton/Data/PacBio/m64446e_220710_022523.hifi_reads.fastq.gz"
    conda:
        "../envs/PacBio.yaml"
    shell:
        '''
        output_name=`echo {output} | sed "s/\\.gz//"`
        bamtools convert -format fastq -in {input} -out ${{output_name}}
        gzip ${{output_name}}
        '''

rule RUT_1_teloclip_hifi_minimap2:
    input:
        genome = select_input_Chr_genome(),
        hifi_reads = "/data/Cotton/Data/PacBio/m64446e_220710_022523.hifi_reads.fastq.gz"
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/m64446e_220710_022523".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    threads: 5
    shell:
        '''
        module load minimap2/2.23
        minimap2 -ax map-hifi -t {threads} {input.genome} {input.hifi_reads} | \
          samtools view -h -F 0x2308 | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''

rule RUT_2_teloclip_ont_minimap2:
    input:
        unpack(get_ONT_ONTSupPassFastq),
        genome = select_input_Chr_genome()
    output:
        directory("results/{prefix}/Recovery_Unassembled_Telomeres/teloclip_ONT_minimap2/{sample_ONTSupPassFastq}-{unit_ONTSupPassFastq}")
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    threads: 10
    shell:
        '''
        module load minimap2/2.23
        minimap2 -ax map-ont -t {threads} {input.genome} {input.ONTSupPassFastq} | \
          samtools view -h -F 0x2308 | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''


rule RUT_3_teloclip_hifi_minimap2:
    input:
        genome = select_input_Chr_genome(),
        pb2genome = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_pbSorted.bam".format(config["prefix"], config["prefix"]),
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/teloclip_hifi_minimap2/pbSorted_bam".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    shell:
        '''
        samtools view -h -F 0x2308 {input.pb2genome} | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''

rule RUT_4_teloclip_hifi_T2Twinnowmap:
    input:
        genome = select_input_Chr_genome(),
        pb2genome = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_hifi2genome.sam".format(config["prefix"], config["prefix"])
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/teloclip_hifi_T2Twinnowmap/hifi2genome_sam".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    shell:
        '''
        samtools view -h -F 0x2308 {input.pb2genome} | awk '!/SA:/ {{print $0;}}' | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''

rule RUT_5_teloclip_ont_minimap2:
    input:
        genome = select_input_Chr_genome(),
        ont2genome = "results/{}/Assessment_Finally_Assembly/Subreads_Coverage/{}_ONTSorted.bam".format(config["prefix"], config["prefix"])
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/teloclip_ONT_minimap2/ONTSorted_bam".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    shell:
        '''
        samtools view -h -F 0x2308 {input.ont2genome} | awk '!/SA:/ {{print $0;}}' | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''

rule RUT_6_teloclip_ont_T2Twinnowmap:
    input:
        genome = select_input_Chr_genome(),
        ont2genome = "results/{}/Assessment_Finally_Assembly/T2T_Alignment/{}_ont2genome.sam".format(config["prefix"], config["prefix"])
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/teloclip_ONT_T2Twinnowmap/ont2genome_sam".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    shell:
        '''
        samtools view -h -F 0x2308 {input.ont2genome} | awk '!/SA:/ {{print $0;}}' | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''

rule RUT_7_teloclip_genome:
    input:
        genome = select_input_Chr_genome(),
        asm_genome = "/public/home/Cotton/Jin668_CM.bp.p_ctg.fasta"
    output:
        directory("results/{}/Recovery_Unassembled_Telomeres/asm_genome/CM.bp.p_ctg".format(config["prefix"]))
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        telomeres_motifs = config["params"]["telomeres_motifs"]
    conda:
        "../envs/Recovery_Unassembled_Telomeres.yaml"
    threads: 10
    shell:
        '''
        module load minimap2/2.23
        minimap2 -ax asm5 -t {threads} {input.genome} {input.asm_genome} | \
          samtools view -h -F 0x2308 | \
          teloclip --ref {input.genome}.fai --motifs {params.telomeres_motifs} | \
          teloclip-extract --refIdx {input.genome}.fai --prefix {params.prefix} --extractReads --extractDir {output}
        '''
#######################################################################################################################################################################
##################################################################### List ############################################################################################
#######################################################################################################################################################################
rule RUT_8_teloclip_list:
    input:
        "results/{}/Recovery_Unassembled_Telomeres/teloclip_hifi_minimap2/pbSorted_bam".format(config["prefix"]),
        "results/{}/Recovery_Unassembled_Telomeres/teloclip_ONT_minimap2/ONTSorted_bam".format(config["prefix"]),
        "results/{}/Recovery_Unassembled_Telomeres/teloclip_hifi_T2Twinnowmap/hifi2genome_sam".format(config["prefix"]),
        "results/{}/Recovery_Unassembled_Telomeres/teloclip_ONT_T2Twinnowmap/ont2genome_sam".format(config["prefix"]),
        "results/{}/Recovery_Unassembled_Telomeres/asm_genome/CM.bp.p_ctg".format(config["prefix"])
