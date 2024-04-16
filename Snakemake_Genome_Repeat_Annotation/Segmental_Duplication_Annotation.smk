####################################################################
####################### TRF MASKER #################################
####################################################################
rule SDA_1_run_trf:
    input:
        fasta = config["chr_genome"]
    output:
        dat = temp("results/{}/Segmental_Duplication_Annotation/{}_TRF.dat".format(config["prefix"], config["prefix"]))
    threads: 1
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        """
        trf {input.fasta} 2 7 7 80 10 50 15 -l 25 -h -ngs > {output.dat}
        """

rule SDA_2_trf_bed:
    """
    workflow/scripts/trf_to_bed.py
    """
    input:
        dat = "results/{}/Segmental_Duplication_Annotation/{}_TRF.dat".format(config["prefix"], config["prefix"])
    output:
        bed = "results/{}/Segmental_Duplication_Annotation/{}_TRF.bed".format(config["prefix"], config["prefix"])
    threads: 1
    shell:
        """
        python workflow/scripts/trf_to_bed.py --infiles {input.dat} --outfile {output.bed}
        """

####################################################################
#################### RepeatMasker bed ##############################
####################################################################
rule SDA_3_RepeatMasker2Bed:
    input:
        soft_out = "results/{}/RepeatMasker/softmasked/{}.out".format(config["prefix"], input_genome_name())
    output:
        bed = temp("results/{}/Segmental_Duplication_Annotation/{}_rm.bed".format(config["prefix"], input_genome_name()))
    shell:
        """
        sed -i "s/?//g" {input.soft_out}
        python3 workflow/scripts/RM2Bed.py -d $(dirname {output.bed}) {input.soft_out}
        """

rule SDA_4_RepeatMasker2Bed_sort:
    input:
        "results/{}/Segmental_Duplication_Annotation/{}_rm.bed".format(config["prefix"], input_genome_name())
    output:
        "results/{}/Segmental_Duplication_Annotation/{}_sRMs.bed".format(config["prefix"], config["prefix"])
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        """
        bedtools sort -i {input} > {output}
        """

###################################################################
########### Make masked version of the fasta for sedef ############
###################################################################
rule SDA_5_pro_masked_fasta:
    input:
        fasta = config["chr_genome"], 
        trf = "results/{}/Segmental_Duplication_Annotation/{}_TRF.bed".format(config["prefix"], config["prefix"]),
        rm = "results/{}/Segmental_Duplication_Annotation/{}_sRMs.bed".format(config["prefix"], config["prefix"])
    output:
        bed = temp("results/{}/Segmental_Duplication_Annotation/{}.tmp.msk.bed".format(config["prefix"], config["prefix"])),
        bed2 = temp("results/{}/Segmental_Duplication_Annotation/{}.tmp2.msk.bed".format(config["prefix"], config["prefix"])),
        maksed_fasta_pro = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta".format(config["prefix"], config["prefix"]),
        maksed_fasta_pro_fai = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta.fai".format(config["prefix"], config["prefix"]),
        hard_maksed_fasta = "results/{}/Segmental_Duplication_Annotation/{}_hard_masked.fasta".format(config["prefix"], config["prefix"])
    threads:1
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        """
        # it is very common to find a 27 base pair gap between alpha sat annotations similarly it is commone to find a 49 bp gap in anotations in HSAT arrays;
        # therefor I merge some of the sat features from rm;
        # additionally there are some huge HSAT arrays that are not annotated as HSAT by rm and instead are Simple_repeat;


        # merge large simple_repeats
        \\grep Simple_repeat {input.rm} | \
                awk '$3-$2 > 1000 {{print $0}}' | \
                bedtools merge -d 100 -i - | \
                grep -v -i "scaffold" | \
                bedtools slop -b 100 -g {input.fasta}.fai -i - >> {output.bed}

        # combine custome merges with trf and rm
        cat {input.trf} {input.rm} | cut -f 1-3 | grep -v -i "scaffold" >> {output.bed}

        # make large merge where large entries are allowed to merge together further
        cut -f 1-3 {output.bed} | bedtools sort -i - | bedtools merge -i - | \
                awk '$3-$2 > 2000 {{print $0}}' | \
                bedtools merge -d 100 -i - > {output.bed2}

        cut -f 1-3 {output.bed} {output.bed2} | bedtools sort -i - | bedtools merge -i - | \
            seqtk seq -l 50 -M /dev/stdin {input.fasta} | sed "s/_//g" > {output.maksed_fasta_pro}

        samtools faidx {output.maksed_fasta_pro}

        cut -f 1-3 {output.bed} {output.bed2} | bedtools sort -i - | bedtools merge -i - | \
            seqtk seq -l 50 -n N -M /dev/stdin {input.fasta} | sed "s/_//g" > {output.hard_maksed_fasta}
        
        samtools faidx {output.hard_maksed_fasta}
        """

######################################################################################################################################################
#################################################### sedef https://github.com/vpc-ccg/sedef ##########################################################
######################################################################################################################################################
rule SDA_6_run_sedef:
    input:
        maksed_fasta_pro = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta".format(config["prefix"], config["prefix"])
    output:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}_sedef.bed".format(config["prefix"], config["prefix"])
    threads: 10
    params:
        workdir = config["workdir"]
    shell:
        """
        module load GCC/9.4.0
        module load parallel/20180222
        export PATH="$PATH:{params.workdir}/workflow/bin/sedef-master"
        bash workflow/bin/sedef-master/sedef.sh -f -o $(dirname {output.bed}) -j {threads} {input.maksed_fasta_pro}
        mv $(dirname {output.bed})/final.bed {output.bed}
        """

# unnessisary can calculate from uppercase matches.
rule SDA_7_count_sat_sedef:
    input:
        rm = "results/{}/Segmental_Duplication_Annotation/{}_sRMs.bed".format(config["prefix"], config["prefix"]),
        sedef = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}_sedef.bed".format(config["prefix"], config["prefix"])
    output:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sat.count.bed".format(config["prefix"], config["prefix"]),
        tmp = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sat.count.tmp.bed".format(config["prefix"], config["prefix"])
    threads: 1
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        """
        cat {input.rm} | \\grep Satellite | cut -f 1-3 | bedtools sort -i - | bedtools merge -i - | sed "s/_//g" | bedtools coverage -header -a {input.sedef} -b - > {output.tmp}
        \\grep "^#" {input.sedef} > {output.bed}
        cat {output.tmp} >> {output.bed}
        sed -i '1{{s/$/\tcount_ovls\tsat_bases\ttotal_bases\tsat_coverage/}}' {output.bed}
        """

rule SDA_8_genome_centromere:
    input:
        cen = "results/{}/Centromere_CENH3/{}.centromere.bed".format(config["prefix"], config["prefix"])
    output:
        cen = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.centromere.bed".format(config["prefix"], config["prefix"])
    shell:
        """
        cat {input.cen} | sed "s/_//g" > {output.cen}
        """

from datetime import date
today = date.today()
DATE =  today.strftime("%Y/%m/%d")

rule SDA_9_sedef_browser:
    input:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sat.count.bed".format(config["prefix"], config["prefix"]),
        fai = config["chr_genome"] + ".fai",
        cen = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.centromere.bed".format(config["prefix"], config["prefix"])
    output:
        fai = temp(config["chr_genome"] + "2.fai"),
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.bed".format(config["prefix"], config["prefix"]),
        lowid = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.lowid.bed".format(config["prefix"], config["prefix"]),
        html = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.html".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"]
    run:
        html = open(f"workflow/report/sedef.html").read()
        open(output["html"], "w+").write(html.format(DATE=DATE, SM={params.prefix}))
        shell("""
            cat {input.fai} | sed "s/_//g" > {output.fai}
            python workflow/scripts/sedef_to_bed.py \
            --fai {output.fai} --cens {input.cen} \
            --sat 0.70 --peri 5000000 --telo 500000 \
            {input.bed} {output.bed} {output.lowid}
        """)

rule SDA_10_sedef_bb:
    input:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.bed".format(config["prefix"], config["prefix"]),
        lowid = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.lowid.bed".format(config["prefix"], config["prefix"]),
        fai = config["chr_genome"] + "2.fai",
    output:
        bb = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.bed.bb".format(config["prefix"], config["prefix"]),
        lowid = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.lowid.bed.bb".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        ucsc = config["software"]["UCSCscripts"]
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        """
        {params.ucsc}
        bedcols=`cat {input.bed} | awk -F"\\t" '{{print NF}}' | sort -V | uniq`
        bedToBigBed -type=bed9+${{bedcols}} -tab {input.bed} {input.fai} {output.bb}
        lowidcols=`cat {input.lowid} | awk -F"\\t" '{{print NF}}' | sort -V | uniq`
        bedToBigBed -type=bed9+${{lowidcols}} -tab {input.lowid} {input.fai} {output.lowid}
        """

rule SDA_11_sum_sedef:
    input:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.bed".format(config["prefix"], config["prefix"]),
        lowid = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.lowid.bed".format(config["prefix"], config["prefix"]),
        cen = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.centromere.bed".format(config["prefix"], config["prefix"]),
        fai = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta.fai".format(config["prefix"], config["prefix"])
    output:
        xlsx = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sedef.summary.xlsx".format(config["prefix"], config["prefix"]),
        lowid = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sedef.lowid.summary.xlsx".format(config["prefix"], config["prefix"])
    threads:1
    shell:
        """
        python workflow/scripts/sedef_summary.py --fai {input.fai} --cen {input.cen} --excel {output.xlsx} {input.bed}
        python workflow/scripts/sedef_summary.py --fai {input.fai} --cen {input.cen} --excel {output.lowid} {input.lowid}
        """

rule SDA_12_enriched:
    input:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.SDs.bed".format(config["prefix"], config["prefix"]),
        fai = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta.fai".format(config["prefix"], config["prefix"])
    output:
        bed = "results/{}/Segmental_Duplication_Annotation/SEDEF/{}.sedef.enriched.bed".format(config["prefix"], config["prefix"])
    threads:1
    shell:
        """
        python workflow/scripts/enriched.py {input.bed} {input.fai} > {output.bed}
        """

######################################################################################################################################################
################################################ biser https://github.com/0xTCG/biser ######################################################
######################################################################################################################################################
rule SDA_13_biser:
    input:
        #hard_genome = "results/{}/RepeatMasker/hardmasked/{}.masked".format(config["prefix"], input_genome_name())
        maksed_fasta_pro = "results/{}/Segmental_Duplication_Annotation/{}_hard_masked.fasta".format(config["prefix"], config["prefix"])
        #maksed_fasta_pro = "results/{}/Segmental_Duplication_Annotation/{}_masked.fasta".format(config["prefix"], config["prefix"])
    output:
        bed = "results/{}/Segmental_Duplication_Annotation/BISER/hard/{}_BISER_SDs.bedpe".format(config["prefix"], config["prefix"]),
        tmdir = temp(directory("results/{}/Segmental_Duplication_Annotation/BISER/hard/TEMP".format(config["prefix"], config["prefix"])))
    threads: 8
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        """
        mkdir {output.tmdir}
        biser --temp {output.tmdir} --hard --gc-heap 10G --max-error=20 --max-edit-error=10 --threads {threads} --output {output.bed} --keep-temp {input.maksed_fasta_pro}
        """

rule SDA_15_biser_SD_size:
    input:
        "results/{}/Segmental_Duplication_Annotation/BISER/{}_BISER_SDs_size_identity_stat.txt".format(config["prefix"], config["prefix"])
    output:
        "results/{}/Segmental_Duplication_Annotation/BISER/{}_BISER_SDs_size_identity_stat.pdf".format(config["prefix"], config["prefix"])
    shell:
        """
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        Rscript workflow/scripts/biser_SD_size.R {input} {output}
        """
