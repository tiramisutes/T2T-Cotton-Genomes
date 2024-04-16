if config["split_genome_fasta"]:
    rule SGS_1_genome_oneline:
        input:
            ref=FASTA
        output:
            temp("temp/{SM}/{SM}.oneline")
        shell:
            '''
            cat {input.ref} | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"%":$0 }}' > {output}
            '''

    rule SGS_2_genome_oneline_split:
        input:
            rules.SGS_1_genome_oneline.output
        output:
            fa = temp("temp/{SM}/{SM}.{chromosome}.fa"),
            fai = temp("temp/{SM}/{SM}.{chromosome}.fa.fai")
        threads: 1
        params:
            chromosome = "{chromosome}"
        shell:
            '''
            if [[ {params.chromosome} == *"Scaffold"* ]]
            then
                cat {input} | awk -F"%" -v chr="{params.chromosome}" 'BEGIN{{OFS="\\n"}} $1 ~ chr {{print $1,$2}}' > {output.fa}
                samtools faidx {output.fa}
            else
                cat {input} | \\grep -w "{params.chromosome}" | sed "s/%/\\n/" > {output.fa}
                samtools faidx {output.fa}
            fi
            '''

    rule SGS_3_make_windows:
        input:
            fai=rules.SGS_2_genome_oneline_split.output.fai
        output:
            bed=temp("temp/{SM}/{SM}.{chromosome}.{W}.bed")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=1
        params:
            slide=f"-s {SLIDE}" if SLIDE > 0 else ""
        shell:
            """
            bedtools makewindows -g {input.fai} -w {wildcards.W} {params.slide} > {output.bed}
            """


    rule SGS_4_split_windows:
        input:
            bed=rules.SGS_3_make_windows.output.bed
        output:
            bed=temp(expand("temp/{SM}/{SM}.{chromosome}.{W}.{ID}.bed", ID=IDS, allow_missing=True))
        threads: 1
        conda:
            "../envs/env.yaml"
        resources:
            mem=1
        shell:
            """
            python workflow/scripts/batch_bed_files.py {input.bed} --outputs {output.bed}
            """


    rule SGS_5_window_fa:
        input:
            ref=rules.SGS_2_genome_oneline_split.output.fa,
            bed=rules.SGS_3_make_windows.output.bed
        output:
            fasta=temp("results/{SM}/{SM}.{chromosome}.{W}.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=4
        shell:
            """
            bedtools getfasta -fi {input.ref} -bed {input.bed} > {output.fasta}
            """


    rule SGS_6_split_fa:
        input:
            ref=rules.SGS_2_genome_oneline_split.output.fa
        output:
            fasta=temp("temp/{SM}/{SM}.{chromosome}.{W}.{REF_ID}.split.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            name=lambda wc: names[int(wc.REF_ID.strip("ref_"))]
        resources:
            mem=4
        shell:
            """
            samtools faidx {input.ref} {params.name} > {output.fasta} 
            """


    rule SGS_7_aln_prep:
        input:
            ref=rules.SGS_6_split_fa.output.fasta if SLIDE > 0 else rules.SGS_5_window_fa.output.fasta
        output:
            split_ref_index=temp("temp/{SM}/{SM}.{chromosome}.{W}.{F}.{REF_ID}.fasta.mmi")
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            S=S,
            MAP_PARAMS=MAP_PARAMS
        shell:
            """
            minimap2 \
                -f {wildcards.F} -s {params.S} \
                {params.MAP_PARAMS} \
                -d {output.split_ref_index} \
                {input.ref}
            """


    rule SGS_8_query_prep:
        input:
            ref=rules.SGS_2_genome_oneline_split.output.fa,
            bed="temp/{SM}/{SM}.{chromosome}.{W}.{ID}.bed"
        output:
            query_fasta=temp("temp/{SM}/{SM}.{chromosome}.{W}.{ID}.query.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=4
        shell:
            """
            bedtools getfasta -fi {input.ref} -bed {input.bed} > {output.query_fasta}
            """


    rule SGS_9_aln:
        input:
            split_ref=rules.SGS_7_aln_prep.output.split_ref_index,
            query=rules.SGS_8_query_prep.output.query_fasta
        output:
            aln=temp("temp/{SM}/{SM}.{chromosome}.{W}.{F}.{ID}.{REF_ID}.bam")
        conda:
            "../envs/env.yaml"
        threads: ALN_T
        params:
            S=S,
            MAP_PARAMS=MAP_PARAMS
        shell:
            """
            minimap2 \
                -t {threads} \
                -f {wildcards.F} -s {params.S} \
                {params.MAP_PARAMS} \
                --dual=yes --eqx \
                {input.split_ref} {input.query} \
                    | samtools sort -o {output.aln}
            """


    rule SGS_10_merge_list:
        input:
            aln=expand(
                "temp/{SM}/{SM}.{chromosome}.{W}.{F}.{ID}.{REF_ID}.bam",
                ID=IDS,
                REF_ID=REF_IDS,
                allow_missing=True,
            ),
        output:
            alns=temp("temp/{SM}/{SM}.{chromosome}.{W}.{F}.list")
        threads: 1
        resources:
            mem=8
        run:
            open(output.alns, "w").write("\n".join(input.aln) + "\n")


    rule SGS_11_merge_aln:
        input:
            alns=rules.SGS_10_merge_list.output.alns,
            aln=expand(
                "temp/{SM}/{SM}.{chromosome}.{W}.{F}.{ID}.{REF_ID}.bam",
                REF_ID=REF_IDS,
                ID=IDS,
                allow_missing=True,
            ),
        output:
            aln=temp("temp/{SM}/{SM}.{chromosome}.{W}.{F}.bam")
        conda:
            "../envs/env.yaml"
        threads: 4
        resources:
            mem=4,
        shell:
            """
            #samtools cat -b {input.alns} -o {output.aln} 
            samtools merge -b {input.alns} {output.aln} 
            """


    rule SGS_12_sort_aln:
        input:
            aln=rules.SGS_11_merge_aln.output.aln
        output:
            aln="results/{SM}/{SM}.{chromosome}.{W}.{F}.sorted.bam"
        conda:
            "../envs/env.yaml"
        threads: 8
        resources:
            mem=SAMTOOLS_MEM
        shell:
            """
            samtools sort -m {resources.mem}G -@ {threads} --write-index \
                -o {output.aln} {input.aln}
            """


    rule SGS_13_identity:
        input:
            aln=rules.SGS_12_sort_aln.output.aln
        output:
            tbl=temp("temp/{SM}/{SM}.{chromosome}.{W}.{F}.tbl.gz")
        conda:
            "../envs/env.yaml"
        threads: 8
        params:
            S=S
        resources:
            mem=8,
        threads: 8
        shell:
            """
            python workflow/scripts/samIdentity.py --threads {threads} \
                --matches  {params.S} --header \
                {input.aln} \
                | pigz -p {threads} > {output.tbl}
            """


    rule SGS_14_pair_end_bed:
        input:
            tbl=rules.SGS_13_identity.output.tbl,
            fai=rules.SGS_2_genome_oneline_split.output.fai
        output:
            bed="results/{SM}/{SM}.{chromosome}.{W}.{F}.bed.gz",
            full="results/{SM}/{SM}.{chromosome}.{W}.{F}.full.tbl.gz"
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            one="--one" if SLIDE > 0 else ""
        shell:
            """
            python workflow/scripts/refmt.py \
                --window {wildcards.W} --fai {input.fai} \
                --full {output.full} \
                {params.one} \
                {input.tbl} {output.bed}
            """

    rule SGS_15_pair_end_bed_merge:
        input:
            bed=expand("results/{SM}/{SM}.{chromosome}.{W}.{F}.bed.gz", chromosome = config["chromosome"], allow_missing=True)
        output:
            bed="results/{SM}/{SM}.{W}.{F}.bed.gz"
        shell:
            """
            zcat {input.bed} | head -n 1 > {output.bed}.header
            zcat {input.bed} | grep -v "^#" > {output.bed}.tmp
            cat {output.bed}.header {output.bed}.tmp | gzip > {output.bed}
            \\rm {output.bed}.header {output.bed}.tmp
            """

else:

    rule SG_3_make_windows:
        input:
            fai=FAI
        output:
            bed=temp("temp/{SM}/{SM}.{W}.bed")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=1
        params:
            slide=f"-s {SLIDE}" if SLIDE > 0 else ""
        shell:
            """
            bedtools makewindows -g {input.fai} -w {wildcards.W} {params.slide} > {output.bed}
            """


    rule SG_4_split_windows:
        input:
            bed=rules.SG_3_make_windows.output.bed
        output:
            bed=temp(expand("temp/{SM}/{SM}.{W}.{ID}.bed", ID=IDS, allow_missing=True))
        threads: 1
        conda:
            "../envs/env.yaml"
        resources:
            mem=1
        shell:
            """
            python workflow/scripts/batch_bed_files.py {input.bed} --outputs {output.bed}
            """


    rule SG_5_window_fa:
        input:
            ref=FASTA,
            bed=rules.SG_3_make_windows.output.bed
        output:
            fasta=temp("results/{SM}/{SM}.{W}.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=4
        shell:
            """
            bedtools getfasta -fi {input.ref} -bed {input.bed} > {output.fasta}
            """


    rule SG_6_split_fa:
        input:
            ref=FASTA
        output:
            fasta=temp("temp/{SM}/{SM}.{W}.{REF_ID}.split.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            name=lambda wc: names[int(wc.REF_ID.strip("ref_"))]
        resources:
            mem=4
        shell:
            """
            samtools faidx {input.ref} {params.name} > {output.fasta} 
            """


    rule SG_7_aln_prep:
        input:
            ref=rules.SG_6_split_fa.output.fasta if SLIDE > 0 else rules.SG_5_window_fa.output.fasta
        output:
            split_ref_index=temp("temp/{SM}/{SM}.{W}.{F}.{REF_ID}.fasta.mmi")
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            S=S,
            MAP_PARAMS=MAP_PARAMS
        shell:
            """
            minimap2 \
                -f {wildcards.F} -s {params.S} \
                {params.MAP_PARAMS} \
                -d {output.split_ref_index} \
                {input.ref}
            """


    rule SG_8_query_prep:
        input:
            ref=FASTA,
            bed="temp/{SM}/{SM}.{W}.{ID}.bed"
        output:
            query_fasta=temp("temp/{SM}/{SM}.{W}.{ID}.query.fasta")
        conda:
            "../envs/env.yaml"
        threads: 1
        resources:
            mem=4
        shell:
            """
            bedtools getfasta -fi {input.ref} -bed {input.bed} > {output.query_fasta}
            """


    rule SG_9_aln:
        input:
            split_ref=rules.SG_7_aln_prep.output.split_ref_index,
            query=rules.SG_8_query_prep.output.query_fasta
        output:
            aln=temp("temp/{SM}/{SM}.{W}.{F}.{ID}.{REF_ID}.bam")
        conda:
            "../envs/env.yaml"
        threads: ALN_T
        params:
            S=S,
            MAP_PARAMS=MAP_PARAMS
        shell:
            """
            minimap2 \
                -t {threads} \
                -f {wildcards.F} -s {params.S} \
                {params.MAP_PARAMS} \
                --dual=yes --eqx \
                {input.split_ref} {input.query} \
                    | samtools sort -o {output.aln}
            """


    rule SG_10_merge_list:
        input:
            aln=expand(
                "temp/{SM}/{SM}.{W}.{F}.{ID}.{REF_ID}.bam",
                ID=IDS,
                REF_ID=REF_IDS,
                allow_missing=True,
            ),
        output:
            alns=temp("temp/{SM}/{SM}.{W}.{F}.list")
        threads: 1
        resources:
            mem=8
        run:
            open(output.alns, "w").write("\n".join(input.aln) + "\n")


    rule SG_11_merge_aln:
        input:
            alns=rules.SG_10_merge_list.output.alns,
            aln=expand(
                "temp/{SM}/{SM}.{W}.{F}.{ID}.{REF_ID}.bam",
                REF_ID=REF_IDS,
                ID=IDS,
                allow_missing=True,
            ),
        output:
            aln=temp("temp/{SM}/{SM}.{W}.{F}.bam")
        conda:
            "../envs/env.yaml"
        threads: 4
        resources:
            mem=4,
        shell:
            """
            #samtools cat -b {input.alns} -o {output.aln} 
            samtools merge -b {input.alns} {output.aln} 
            """


    rule SG_12_sort_aln:
        input:
            aln=rules.SG_11_merge_aln.output.aln
        output:
            aln="results/{SM}/{SM}.{W}.{F}.sorted.bam"
        conda:
            "../envs/env.yaml"
        threads: 8
        resources:
            mem=SAMTOOLS_MEM
        shell:
            """
            samtools sort -m {resources.mem}G -@ {threads} --write-index \
                -o {output.aln} {input.aln}
            """


    rule SG_13_identity:
        input:
            aln=rules.SG_12_sort_aln.output.aln
        output:
            tbl=temp("temp/{SM}/{SM}.{W}.{F}.tbl.gz")
        conda:
            "../envs/env.yaml"
        threads: 8
        params:
            S=S
        resources:
            mem=8,
        threads: 8
        shell:
            """
            python workflow/scripts/samIdentity.py --threads {threads} \
                --matches  {params.S} --header \
                {input.aln} \
                | pigz -p {threads} > {output.tbl}
            """


    rule SG_14_pair_end_bed:
        input:
            tbl=rules.SG_13_identity.output.tbl,
            fai=FAI
        output:
            bed="results/{SM}/{SM}.{W}.{F}.bed.gz",
            full="results/{SM}/{SM}.{W}.{F}.full.tbl.gz"
        conda:
            "../envs/env.yaml"
        threads: 1
        params:
            one="--one" if SLIDE > 0 else ""
        shell:
            """
            python workflow/scripts/refmt.py \
                --window {wildcards.W} --fai {input.fai} \
                --full {output.full} \
                {params.one} \
                {input.tbl} {output.bed}
            """

rule SG_16_make_R_figures:
    input:
        bed=rules.SGS_15_pair_end_bed_merge.output.bed if config["split_genome_fasta"] else rules.SG_14_pair_end_bed.output.bed
    output:
        pdf="results/{SM}/{SM}.{W}.{F}_figures/pdfs/{SM}.{W}.{F}.tri.TRUE__onecolorscale.FALSE__all.pdf",
        png=report(
            expand(
                "results/{SM}/{SM}.{W}.{F}_figures/pngs/{SM}.{W}.{F}.tri.{first}__onecolorscale.{second}__all.png",
                allow_missing=True,
                first=["TRUE", "FALSE"],
                second=["TRUE", "FALSE"],
            ),
            category="figures",
        ),
        facet=report(
            "results/{SM}/{SM}.{W}.{F}_figures/pngs/{SM}.{W}.{F}.facet.all.png",
            category="figures",
        ),
    #conda: "../envs/R.yaml"
    params:
        workdir = config["workdir"],
        outdir = lambda w, output: "/".join(os.path.dirname(output[0]).split("/")[0:2])
    threads: 1
    shell:
        """
        module load R/3.6.0
        Rscript workflow/scripts/aln_plot.R \
            --bed {input.bed} \
            --threads {threads} \
            --prefix {wildcards.SM}.{wildcards.W}.{wildcards.F} \
            --outrootdir {params.workdir}/{params.outdir}
        """

# ##############################################################################################################################
# for large scale visualization in https://higlass.io/
# ##############################################################################################################################
rule SG_17_cooler_strand:
    input:
        bed=rules.SGS_15_pair_end_bed_merge.output.bed if config["split_genome_fasta"] else rules.SG_14_pair_end_bed.output.bed,
        fai=FAI
    output:
        cool="results/{SM}/{SM}.{W}.{F}.strand.cool"
    conda:
        "../envs/env.yaml"
    threads: 1
    resources:
        mem=64,
    shell:
        """
        gunzip -c {input.bed} | tail -n +2 \
            | sed  's/+$/100/g' \
            | sed  's/-$/50/g' \
              | cooler cload pairs \
                -c1 1 -p1 2 -c2 4 -p2 5 \
                        --field count=8:agg=mean,dtype=float \
                --chunksize 50000000000 \
                {input.fai}:{wildcards.W} \
                --zero-based \
                - {output.cool}
        """


rule SG_18_cooler_identity:
    input:
        bed=rules.SGS_15_pair_end_bed_merge.output.bed if config["split_genome_fasta"] else rules.SG_14_pair_end_bed.output.bed,
        fai=FAI
    output:
        cool="results/{SM}/{SM}.{W}.{F}.identity.cool"
    conda:
        "../envs/env.yaml"
    threads: 1
    resources:
        mem=64,
    shell:
        """
        gunzip -c {input.bed} | tail -n +2 \
              | cooler cload pairs \
                -c1 1 -p1 2 -c2 4 -p2 5 \
                        --field count=7:agg=mean,dtype=float \
                --chunksize 50000000000 \
                {input.fai}:{wildcards.W} \
                --zero-based \
                - {output.cool}

        """

rule SG_19_cooler_zoomify_s:
    input:
        s=rules.SG_17_cooler_strand.output.cool
    output:
        s="results/{SM}/{SM}.{W}.{F}.strand.mcool"
    conda:
        "../envs/env.yaml"
    threads: 8
    resources:
        mem=8,
    shell:
        """
        cooler zoomify --field count:agg=mean,dtype=float {input.s} -n {threads} -o {output.s}
        """

rule SG_20_cooler_zoomify_i:
    input:
        i=rules.SG_18_cooler_identity.output.cool
    output:
        i="results/{SM}/{SM}.{W}.{F}.identity.mcool"
    conda:
        "../envs/env.yaml"
    threads: 8
    resources:
        mem=8,
    shell:
        """
        cooler zoomify --field count:agg=mean,dtype=float {input.i} -n {threads} -o {output.i}
        """

# ##############################################################################################
# ##############################################################################################
rule SG_21_bwa_index:
    input:
        fasta=FASTA,
    output:
        done=touch("temp/{SM}/bwa_ref_{SM}")
    conda:
        "../envs/env.yaml"
    threads: 1
    shell:
        """
        bwa index -p temp/{SM}/ref_{wildcards.SM} {input.fasta}
        """


rule SG_22_bwa_aln:
    input:
        bwa_index_done="temp/{SM}/bwa_ref_{SM}",
        reads=expand("results/{SM}/{SM}.{chromosome}.{W}.fasta", chromosome = config["chromosome"], allow_missing=True) if config["split_genome_fasta"] else rules.SG_5_window_fa.output.fasta
    output:
        sam="results/{SM}/output_{SM}_{W}.sam.gz"
    conda:
        "../envs/env.yaml"
    params:
        num_dups=NUM_DUPS,
    threads: ALN_T
    shell:
        """
        bwa aln -t {threads} temp/{SM}/ref_{wildcards.SM} {input.reads} \
            | bwa samse -n {params.num_dups} temp/{SM}/ref_{wildcards.SM} - {input.reads} \
            | gzip > {output.sam}
        """


rule SG_23_create_contacts:
    input:
        sam=rules.SG_22_bwa_aln.output.sam
    output:
        contacts="results/{SM}/contacts_{SM}_{W}.gz"
    conda:
        "../envs/env.yaml"
    threads: 1
    shell:
        """
        gunzip -c {input.sam} | python workflow/scripts/ntn_bam_to_contacts.py - | gzip > {output.contacts}
        """


rule SG_24_cooler_density_d:
    input:
        contacts=rules.SG_23_create_contacts.output.contacts,
        fai=FAI,
    output:
        cool="results/{SM}/{SM}.{W}.{CW}.density.cool"
    conda:
        "../envs/env.yaml"
    threads: 1
    resources:
        mem=64,
    shell:
        """
        cooler cload pairs \
                -c1 1 -p1 2 -c2 4 -p2 5 \
                --chunksize 50000000000 \
                {input.fai}:{wildcards.CW} \
                --zero-based \
                {input.contacts} {output.cool}
        """


rule SG_25_cooler_zoomify_d:
    input:
        d=rules.SG_24_cooler_density_d.output.cool
    output:
        d="results/{SM}/{SM}.{W}.{CW}.density.mcool"
    conda:
        "../envs/env.yaml"
    threads: 8
    resources:
        mem=8
    shell:
        """
        cooler zoomify {input.d} -n {threads} -o {output.d}
        """
