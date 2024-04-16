localrules:
    Centromere_Analysis_2_genome_oneline,
    Centromere_Analysis_3_genome_oneline_split,
    Centromere_Analysis_5_RepeatExplorer2_install,
    Centromere_Analysis_6_RepeatExplorer2_install,
    Centromere_Analysis_11_EDTA_Format,
    Centromere_Analysis_14_intact_LTR_internal_domain_Seq,
    Centromere_Analysis_15_centromere_chromosome,
    Centromere_Analysis_17_T2T_Centromere_Kmer_genome_centromere_specific_db,
    Centromere_Analysis_18_T2T_Centromere_Kmer_genome_centromere_specific_db_count,
    Centromere_Analysis_60_list

##################################################################################################################################
#################################################### 1.  ######################################################
# Typically, eukaryotic centromere regions contain repeat units ranging in length from ∼150 to ∼210 bases, approximately the length required to form a single nucleosome. 
##################################################################################################################################
rule Centromere_Analysis_2_whole_pyTanFinder:
    input:
        "results/{}/Centromere_CENH3/{}.centromere.fasta".format(config["prefix"], config["prefix"])
        #config["genome"]
    output:
        trfdat = "results/{}/Centromere_Analysis/Tandem_Repeats_Finder/{}_{}_{}ML/{}.centromere.fasta.2.7.7.80.10.20.{}.dat".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["centromere_satellite_arrays_minMonLength"], config["Centromere_Analysis"]["centromere_satellite_arrays_maxMonLength"], config["prefix"], config["Centromere_Analysis"]["centromere_satellite_arrays_maxMonLength"]),
        repeatAbundance = "results/{}/Centromere_Analysis/Tandem_Repeats_Finder/{}_{}_{}ML/out_pyTanFinder/out_merger/merger_output_{}/RepeatAbundance_{}.txt".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["centromere_satellite_arrays_minMonLength"], config["Centromere_Analysis"]["centromere_satellite_arrays_maxMonLength"], config["prefix"], config["prefix"])
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        centromere_satellite_arrays_minMonLength = config["Centromere_Analysis"]["centromere_satellite_arrays_minMonLength"],
        centromere_satellite_arrays_maxMonLength = config["Centromere_Analysis"]["centromere_satellite_arrays_maxMonLength"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/Centromere_Analysis_pyTanFinder.yaml"
    shell:
        '''
        module load TRF/4.0.9
        trfpath=`which trf`

        cd {params.outdir}
        python {params.workdir}/workflow/bin/pyTanFinder-master/pyTanFinder.py --minMonLength {params.centromere_satellite_arrays_minMonLength} \
          --maxMonLength {params.centromere_satellite_arrays_maxMonLength} -px {params.prefix} -tp ${{trfpath}} {params.workdir}/{input}
        '''

##################################################################################################################################
################################################## 2.  ##################################################
##################################################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TandemTools for ONT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule Centromere_Analysis_7_meryl_count_genome:
    input:
        genome = config["chr_genome"],
    output:
        merylDB = temp(directory("results/{}/Centromere_Analysis/Read_Coverage/{}_{}kmer.meryl".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["winnowmap_kmer"]))),
        repetitive = temp("results/{}/Centromere_Analysis/Read_Coverage/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["winnowmap_kmer"]))
    threads: 1
    params:
        winnowmap_kmer = config["Centromere_Analysis"]["winnowmap_kmer"]
    conda: "../envs/Centromere_Analysis.yaml"
    shell:
        '''
        meryl count k={params.winnowmap_kmer} output {output.merylDB} {input.genome}
        meryl print greater-than distinct=0.9998 {output.merylDB} > {output.repetitive}
        '''

rule Centromere_Analysis_8_winnowmap_mapping_ont:
    input:
        repetitive = "results/{}/Centromere_Analysis/Read_Coverage/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["winnowmap_kmer"]),
        genome = config["chr_genome"],
        ont = config["Centromere_Analysis"]["ONT_reads"]
    output:
        ontsam = temp("results/{}/Centromere_Analysis/Read_Coverage/{}_ont2genome.sam".format(config["prefix"], config["prefix"]))
    params:
        winnowmap_kmer = config["Centromere_Analysis"]["winnowmap_kmer"]
    threads: 10
    conda: "../envs/Centromere_Analysis.yaml"
    shell:
        '''
        winnowmap -k {params.winnowmap_kmer} -W {input.repetitive} --MD -ax map-ont {input.genome} {input.ont} > {output.ontsam}
        '''

rule Centromere_Analysis_9_winnowmap_mapping_ont_extracted:
    input:
        sam = "results/{}/Centromere_Analysis/Read_Coverage/{}_ont2genome.sam".format(config["prefix"], config["prefix"]),
        bed = "results/{}/Centromere_CENH3/{}.centromere.bed".format(config["prefix"], config["prefix"]),
        ont = config["Centromere_Analysis"]["ONT_reads"]
    output:
        readsId = "results/{}/Centromere_Analysis/Read_Coverage/{}_centromeric_ont_reads.id".format(config["prefix"], config["prefix"]),
        readsFa = "results/{}/Centromere_Analysis/Read_Coverage/{}_centromeric_ont_reads.fasta".format(config["prefix"], config["prefix"])
    threads: 1
    #conda: "../envs/Centromere_Analysis.yaml"
    shell:
        '''
        samtools view -L {input.bed} {input.sam} | cut -f 1 | sort -V | uniq > {output.readsId}
        xargs samtools faidx {input.ont} < {output.readsId} > {output.readsFa}
        '''

rule Centromere_Analysis_10_genome_Chromosome:
    input:
        "results/{}/Centromere_CENH3/{}.centromere_2Mb.fasta".format(config["prefix"], config["prefix"])
    output:
        expand("results/{prefix}/Centromere_Analysis/Read_Coverage/Chromosome_Centromere/{chromosome}.fasta", prefix = config["prefix"], chromosome = config["Centromere_Analysis"]["chromosome"])
    threads: 1
    params:
        workdir = config["workdir"],
        UCSCscripts = config["software"]["UCSCscripts"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        {params.UCSCscripts}
        faSplit byname {params.workdir}/{input} {params.workdir}/{params.outdir}/
        cd {params.workdir}/{params.outdir}/
        for i in `ls | xargs`
        do
        mv ${{i}} ${{i}}sta
        done
        '''

rule Centromere_Analysis_11_winnowmap_mapping_tandemquast_ont_Chromosome:
    input:
        ont_reads = "results/{prefix}/Centromere_Analysis/Read_Coverage/{prefix}_centromeric_ont_reads.fasta",
        genome = "results/{prefix}/Centromere_Analysis/Read_Coverage/Chromosome_Centromere/{chromosome}.fasta"
    output:
        directory("results/{prefix}/Centromere_Analysis/Read_Coverage/TandemTools/TandemQUAST_ONT_{chromosome}")
    threads: 10
    conda: "../envs/Centromere_Analysis_TandemTools.yaml"
    shell:
        '''
        module load Jellyfish/2.3.0
        export MPLBACKEND='Agg'
        python workflow/bin/TandemTools/tandemquast.py -t {threads} --nano {input.ont_reads} -o {output} {input.genome}
        '''

rule Centromere_Analysis_12_winnowmap_mapping_tandemmapper_ont_list:
    input:
        expand("results/{prefix}/Centromere_Analysis/Read_Coverage/TandemTools/TandemQUAST_ONT_{chromosome}", prefix = config["prefix"], chromosome = config["Centromere_Analysis"]["chromosome"])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NucFreq for HiFi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule Centromere_Analysis_13_Read_Coverage_hifi:
    input:
        repetitive = "results/{}/Centromere_Analysis/Read_Coverage/{}_repetitive_k{}.txt".format(config["prefix"], config["prefix"], config["Centromere_Analysis"]["winnowmap_kmer"]),
        genome = config["chr_genome"],
        hifi = config["Centromere_Analysis"]["hifi_reads"]
    output:
        hifisam = "results/{}/Centromere_Analysis/Read_Coverage/NucFreq/{}.hifi.winnowmap.bam".format(config["prefix"], config["prefix"])
    threads: 10
    params:
        winnowmap_kmer = config["Centromere_Analysis"]["winnowmap_kmer"]
    conda: "../envs/Centromere_Analysis.yaml"
    shell:
        '''
        winnowmap -k {params.winnowmap_kmer} -W {input.repetitive} --MD -Ha -ax map-pb {input.genome} {input.hifi} > {output.hifisam}
        '''

rule Centromere_Analysis_14_winnowmap_mapping_filter_hifi:
    input:
        hifi_sam = "results/{}/Centromere_Analysis/Read_Coverage/NucFreq/{}.hifi.winnowmap.bam".format(config["prefix"], config["prefix"])
    output:
        hifi_bam = "results/{}/Centromere_Analysis/Read_Coverage/NucFreq/{}.hifi.winnowmap_F2308.bam".format(config["prefix"], config["prefix"])
    threads: 10
    shell:
        '''
        samtools view -@ {threads} -F 2308 -bS {input.hifi_sam} | samtools sort -@ {threads} -o {output.hifi_bam} -
        samtools index -@ {threads} {output.hifi_bam} 
        '''

rule Centromere_Analysis_15_centromere_chromosome:
    input:
        region_bed = "results/{}/Centromere_CENH3/{}.centromere.bed".format(config["prefix"], config["prefix"])
    output:
        temp("results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{chromosome}.centromere.bed")
    params:
        chromosome = "{chromosome}"
    shell:
        '''
        cat {input.region_bed} | \\grep -w {params.chromosome} > {output}
        '''

rule Centromere_Analysis_16_winnowmap_mapping_filter_hifi_NucFreq:
    input:
        hifi_bam = "results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{prefix}.hifi.winnowmap_F2308.bam",
        region_bed = "results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{chromosome}.centromere.bed"
    output:
        bed = "results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{chromosome}.hifi.centromere_points.bed",
        plotpng = "results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{chromosome}.hifi.centromere.png"
    threads: 2
    params:
        chromosome = "{chromosome}"
    shell:
        '''
        python workflow/scripts/NucPlot.py --height 6 --width 16 --ylim 100 --threads {threads} --obed {output.bed} --bed {input.region_bed} --minobed 2 {input.hifi_bam} {output.plotpng}
        '''

rule Centromere_Analysis_16_winnowmap_mapping_hifi_NucFreq_list:
    input:
        expand("results/{prefix}/Centromere_Analysis/Read_Coverage/NucFreq/{chromosome}.hifi.centromere.png", prefix = config["prefix"], chromosome = config["Centromere_Analysis"]["chromosome"])

##################################################################################################################################
########################################################## 3. TE #############################################################
##################################################################################################################################
rule Centromere_Analysis_18_centromere_TEs:
    input:
        intact_gff3 = "results/{{prefix}}/EDTA/{}.mod.EDTA.intact.gff3".format(input_genome_name()),
        cent_bed = "results/{prefix}/Centromere_CENH3/{prefix}.{centype}.bed"
    output:
        cent_tes = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_{centype}_intact_TEs_LTR.bed"
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        '''
        cat {input.intact_gff3} | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $1,$4,$5,$3,$6,$7}}' | \
          grep -v -E "Low_complexity|Simple_repeat|Satellite|target_site_duplication|long_terminal_repeat|repeat_region" | sort -k1,1V -k2,2n | uniq | \
          \\grep "LTR" | sort -V | uniq | \
          bedtools intersect -a {input.cent_bed} -b - -wb | awk -F"\\t" 'BEGIN{{OFS="\\t"}} {{print $4,$5,$6,$7}}' > {output.cent_tes}
        '''

rule Centromere_Analysis_19_centromere_TEs_LTR_Lineages_fasta:
    input:
        intact_LTR_bed = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_{centype}_intact_TEs_LTR.bed",
        genome = config["genome"]
    output:
        cent_ltr_fa = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_{centype}_intact_TEs_LTR.fasta"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    conda: 
        "../envs/Segmental_Duplication_Annotation.yaml"
    shell:
        '''
        bedtools getfasta -fi {input.genome} -bed {input.intact_LTR_bed} -fo {output.cent_ltr_fa}
        sed -i "s|:|_|g" {output.cent_ltr_fa}
        '''

rule Centromere_Analysis_20_centromere_TEs_LTR_DANTE_Finder:
    input:
        cent_ltr_fa = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_{centype}_intact_TEs_LTR.fasta"
    output:
        cent_ltr_gff3 = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains.gff3"
    threads: 1
    params:
        workdir = config["workdir"],
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        module load last/1411
        cd {params.outdir}
        python {params.workdir}/workflow/bin/dante/dante.py --domain_gff {params.workdir}/{output.cent_ltr_gff3} -q {params.workdir}/{input.cent_ltr_fa} \
          -pdb {params.workdir}/workflow/bin/dante/tool-data/protein_domains/Viridiplantae_v3.0_pdb \
          -cs {params.workdir}/workflow/bin/dante/tool-data/protein_domains/Viridiplantae_v3.0_class
        '''

rule Centromere_Analysis_21_centromere_TEs_LTR_DANTE_Filter:
    input:
        cent_ltr_gff3 = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains.gff3"
    output:
        cent_filtered_domains_gff = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}.gff3",
        cent_filtered_domains_fa = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}.fa"
    threads: 1
    params:
        workdir = config["workdir"],
        domain = "{domain}",
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        python workflow/bin/dante/dante_gff_output_filtering.py \
          --dom_gff {input.cent_ltr_gff3} \
          --domains_filtered {output.cent_filtered_domains_gff} \
          --domains_prot_seq {output.cent_filtered_domains_fa} \
          --selected_dom {params.domain}
        '''

rule Centromere_Analysis_22_centromere_TEs_LTR_DANTE_Extract_Domains_Nucleotide_Sequences:
    input:
        cent_ltr_fa = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_{centype}_intact_TEs_LTR.fasta",
        cent_filtered_domains_gff = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}.gff3",
        cent_filtered_domains_fa = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}.fa"
    output:
        cent_ltr_domains = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}_Seq/domains_counts.txt"
    threads: 1
    params:
        workdir = config["workdir"],
        cent_outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        python workflow/bin/dante/dante_gff_to_dna.py --domains_gff {input.cent_filtered_domains_gff} --input_dna {input.cent_ltr_fa} \
          --classification workflow/bin/dante/tool-data/protein_domains/Viridiplantae_v3.0_class --out_dir {params.cent_outdir} --extended True
        '''

rule Centromere_Analysis_23_centromere_TEs_LTR_Domains_Tree:
    input:
        "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}_Seq/domains_counts.txt"
    output:
        Copia_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][0]),
        Copia_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][1]),
        Gypsy_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][0]),
        Gypsy_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][1])
    threads: 5
    params:
        workdir = config["workdir"],
        species = "{prefix}",
        copia_lineages_outgroup = config["Centromere_Analysis"]["copia_lineages_outgroup"],
        gypsy_lineages_outgroup = config["Centromere_Analysis"]["gypsy_lineages_outgroup"],
        indir = lambda w, input: os.path.dirname(input[0]),
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        set +eu
        module load mafft/7.453
        module load iqtree/1.6.12
        module load last/1133
        module load FastTree/2.1.11
        module load R/3.6.0
        copia=`cat {input} | \\grep -i "LTR" | awk -F"|" 'BEGIN{{OFS="\\t"}} {{print $3,$NF}}' | cut -f 2 -d / | cut -f 1 -d : | \\grep -i "copia" | cut -f 2 | sort -V | xargs`
        gypsy=`cat {input} | \\grep -i "LTR" | awk -F"|" 'BEGIN{{OFS="\\t"}} {{print $3,$NF}}' | cut -f 2 -d / | cut -f 1 -d : | \\grep -i "gypsy" | cut -f 2 | sort -V | xargs`
        ########################################################################################
        ## Ty1/Copia
        ########################################################################################
        if [ -z "${{copia}}" ]
        then
            echo "##################################################################"
            echo "##################### Ty1/Copia is empty #########################"
            echo "##################################################################"
            touch {output.Copia_R10_og1} {output.Copia_R10_og2}
        else
            echo "##################################################################"
            echo "############### Ty1/Copia is NOT empty and Run It ################"
            echo "##################################################################"
            if grep -i -q "{params.copia_lineages_outgroup[0]}" <<< "${{copia}}"
            then
                # SIRE Ty1/Copia sequences were included as an outgroup for all other Ty1/Copia trees
                echo "SIRE Ty1/Copia sequences were included as an outgroup for all other Ty1/Copia trees"
                cat {params.indir}/{params.copia_lineages_outgroup[0]}.fasta | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"%":$0 }}' | sed -n '1,10p' | sed "s/%/\\n/g" | sed -e '$a\\' > {output.Copia_R10_og1}
                # For rooting of SIRE Ty1/Copia tree, Angela Ty1/Copia elements were used as outgroup.
                if grep -i -q "{params.copia_lineages_outgroup[1]}" <<< "${{copia}}"
                then
                    cat {params.indir}/{params.copia_lineages_outgroup[1]}.fasta | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"%":$0 }}' | sed -n '1,10p' | sed "s/%/\\n/g" | sed -e '$a\\' > {output.Copia_R10_og2}
                else
                    touch {output.Copia_R10_og2}
                fi
                cat {output.Copia_R10_og2} {params.indir}/{params.copia_lineages_outgroup[0]}.fasta | sed "s/:/_/g" > {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}.fasta
                ##  Aligment and Tree
                seqnum=`grep -c ">" {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}.fasta`
                if [ "${{seqnum}}" -gt 2 ]
                then
                    mafft --globalpair --maxiterate 100 --thread {threads} {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}.fasta > {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta
                    iqtree -s {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta -m TN93 -nt {threads}
                    FastTree -nt -gtr -gamma {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta > {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft.tree
                    Rscript workflow/scripts/ApeTree.R {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_TN93_NJ.pdf {params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_TN93_NJ.tree
                fi

                for i in `echo ${{copia}} | sed "s/ /\\n/g" | grep -v -i "{params.copia_lineages_outgroup[0]}" | sort -V | xargs`
                do
                    echo "##################################################################"
                    echo "##################### Copia: ${{i}} #############################"
                    echo "##################################################################"
                    cat {output.Copia_R10_og1} {params.indir}/${{i}}.fasta | sed "s/:/_/g" > {params.outdir}/Copia_${{i}}.fasta
                    seqnum=`grep -c ">" {params.outdir}/Copia_${{i}}.fasta`
                    if [ "${{seqnum}}" -gt 2 ]
                    then
                        mafft --globalpair --maxiterate 100 --thread {threads} {params.outdir}/Copia_${{i}}.fasta > {params.outdir}/Copia_${{i}}_mafft_aln.fasta
                        iqtree -s {params.outdir}/Copia_${{i}}_mafft_aln.fasta -m TN93 -nt {threads}
                        FastTree -nt -gtr -gamma {params.outdir}/Copia_${{i}}_mafft_aln.fasta > {params.outdir}/Copia_${{i}}_mafft.tree
                        Rscript workflow/scripts/ApeTree.R {params.outdir}/Copia_${{i}}_mafft_aln.fasta {params.outdir}/Copia_${{i}}_mafft_TN93_NJ.pdf {params.outdir}/Copia_${{i}}_mafft_TN93_NJ.tree
                    fi
                done
            else
                echo "The Copia: {params.copia_lineages_outgroup[0]} Ty1/Copia outgroup is miss!!!"
                touch {output.Copia_R10_og1} {output.Copia_R10_og2}
                for i in ${{copia}}
                    do
                        echo "##################################################################"
                        echo "##################### Copia: ${{i}} #############################"
                        echo "##################################################################"
                        seqnum=`grep -c ">" {params.indir}/${{i}}.fasta`
                        if [ "${{seqnum}}" -gt 2 ]
                        then
                            mafft --globalpair --maxiterate 100 --thread {threads} {params.indir}/${{i}}.fasta > {params.outdir}/Copia_${{i}}_mafft_aln.fasta
                            iqtree -s {params.outdir}/Copia_${{i}}_mafft_aln.fasta -m TN93 -nt {threads}
                            FastTree -nt -gtr -gamma {params.outdir}/Copia_${{i}}_mafft_aln.fasta > {params.outdir}/Copia_${{i}}_mafft.tree
                            Rscript workflow/scripts/ApeTree.R {params.outdir}/Copia_${{i}}_mafft_aln.fasta {params.outdir}/Copia_${{i}}_mafft_TN93_NJ.pdf {params.outdir}/Copia_${{i}}_mafft_TN93_NJ.tree
                        fi
                    done
            fi
        fi
        ########################################################################################
        ## Ty3/Gypsy
        ########################################################################################
        if [ -z "${{gypsy}}" ]
        then
            echo "##################################################################"
            echo "##################### Ty3/Gypsy is empty #########################"
            echo "##################################################################"
            touch {output.Gypsy_R10_og1} {output.Gypsy_R10_og2}
        else
            echo "##################################################################"
            echo "############### Ty3/Gypsy is NOT empty and Run It ################"
            echo "##################################################################"
            if grep -i -q "{params.gypsy_lineages_outgroup[0]}" <<< "${{gypsy}}"
            then
                cat {params.indir}/{params.gypsy_lineages_outgroup[0]}.fasta | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"%":$0 }}' | sed -n '1,10p' | sed "s/%/\\n/g" | sed -e '$a\\' > {output.Gypsy_R10_og1}
                if grep -i -q "{params.gypsy_lineages_outgroup[1]}" <<< "${{gypsy}}"
                then
                    cat {params.indir}/{params.gypsy_lineages_outgroup[1]}.fasta | awk '/^>/&&NR>1{{print "";}}{{ printf "%s",/^>/ ? $0"%":$0 }}' | sed -n '1,10p' | sed "s/%/\\n/g" | sed -e '$a\\' > {output.Gypsy_R10_og2}
                else
                    touch {output.Gypsy_R10_og2}
                fi
                cat {output.Gypsy_R10_og2} {params.indir}/{params.gypsy_lineages_outgroup[0]}.fasta | sed "s/:/_/g" > {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}.fasta
                ##  Aligment and Tree
                seqnum=`grep -c ">" {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}.fasta`
                if [ "${{seqnum}}" -gt 2 ]
                then
                    mafft --globalpair --maxiterate 100 --thread {threads} {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}.fasta > {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta
                    iqtree -s {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta -m TN93 -nt {threads}
                    FastTree -nt -gtr -gamma {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta > {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft.tree
                    Rscript workflow/scripts/ApeTree.R {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_TN93_NJ.pdf {params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_TN93_NJ.tree
                fi

                for i in `echo ${{gypsy}} | sed "s/ /\\n/g" | grep -v -i "{params.gypsy_lineages_outgroup[0]}" | sort -V | xargs`
                do
                    echo "##################################################################"
                    echo "##################### Gypsy: ${{i}} #############################"
                    echo "##################################################################"
                    cat {output.Gypsy_R10_og1} {params.indir}/${{i}}.fasta > {params.outdir}/Gypsy_${{i}}.fasta
                    seqnum=`grep -c ">" {params.outdir}/Gypsy_${{i}}.fasta`
                    if [ "${{seqnum}}" -gt 2 ]
                    then
                        mafft --globalpair --maxiterate 100 --thread {threads} {params.outdir}/Gypsy_${{i}}.fasta | sed "s/:/_/g" > {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta
                        iqtree -s {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta -m TN93 -nt {threads}
                        FastTree -nt -gtr -gamma {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta > {params.outdir}/Gypsy_${{i}}_mafft.tree
                        Rscript workflow/scripts/ApeTree.R {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta {params.outdir}/Gypsy_${{i}}_mafft_TN93_NJ.pdf {params.outdir}/Gypsy_${{i}}_mafft_TN93_NJ.tree
                    fi
                done
            else
                echo "The Gypsy: {params.gypsy_lineages_outgroup[0]} Ty3/Gypsy outgroup is miss!!!"
                touch {output.Gypsy_R10_og1} {output.Gypsy_R10_og2}
                for i in ${{gypsy}}
                    do
                        echo "##################################################################"
                        echo "##################### Gypsy: ${{i}} #############################"
                        echo "##################################################################"
                        seqnum=`grep -c ">" {params.indir}/${{i}}.fasta`
                        if [ "${{seqnum}}" -gt 2 ]
                        then
                            mafft --globalpair --maxiterate 100 --thread {threads} {params.indir}/${{i}}.fasta | sed "s/:/_/g" > {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta
                            iqtree -s {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta -m TN93 -nt {threads}
                            FastTree -nt -gtr -gamma {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta > {params.outdir}/Gypsy_${{i}}_mafft.tree
                            Rscript workflow/scripts/ApeTree.R {params.outdir}/Gypsy_${{i}}_mafft_aln.fasta {params.outdir}/Gypsy_${{i}}_mafft_TN93_NJ.pdf {params.outdir}/Gypsy_${{i}}_mafft_TN93_NJ.tree
                        fi
                    done
            fi
        fi
        '''

rule Centromere_Analysis_24_centromere_TEs_LTR_Domains_Time_PATHd8:
    input:        
        Copia_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][0]),
        Copia_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][1]),
        Gypsy_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][0]),
        Gypsy_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][1])
    output:
        "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_Protein_Domains_filtered_{domain}_Seq_Tree_Time/{prefix}_{centype}_TE_{domain}_Lineages_PATHd8.Age.txt"
    params:
        workdir = config["workdir"],
        species = "{prefix}",
        domain = "{domain}",
        indir = lambda w, input: os.path.dirname(input[0]),
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        set +e
        for i in `ls {params.indir} | \\grep "treefile" | sed "s/\\.treefile//" | xargs`
        do
        class=`echo ${{i}} | cut -f 1 -d _`
        lineages=`echo ${{i}} | cut -f 2 -d _`
        aln_seq_length=`bioawk -c fastx '{{ print $name, length($seq) }}' < {params.indir}/${{i}} | cut -f 2 | sed -n '1p'`
        echo "Sequence length = ${{aln_seq_length}};" >> {params.outdir}/${{i}}.PATHd8
        cat {params.indir}/${{i}}.treefile >> {params.outdir}/${{i}}.PATHd8
        workflow/bin/PATHd8/PATHd8 -n {params.outdir}/${{i}}.PATHd8 -r {params.outdir}/${{i}}_PATHd8.txt -pa -pn
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            echo "PATHd8 Error for file: ${{i}}"
        else
            cat {params.outdir}/${{i}}_PATHd8.txt | \\grep "d8 tree" | sed "s/^.*: //g" > {params.outdir}/${{i}}_PATHd8_ultrametric.tre
            cat {params.outdir}/${{i}}_PATHd8.txt | sed '/./!d' | sed 's/[ ][ ]*/ /g' | \
            sed -n '/^d8 tree/,/ *) Rate = MPL/p' | sed '$d' | sed "1,2d" | \
            awk -v species={params.species} -v class=${{class}} -v lineages="${{lineages}}" -v domain="{params.domain}" 'BEGIN{{OFS="\\t"}} $1 ~ lineages {{print species,domain,class,lineages,$4}}' >> {output}
        fi
        done
        sed -i '1 i\\Species\\tType\\tClass\\tLineages\\tAge' {output}
        '''

rule Centromere_Analysis_25_centromere_TEs_LTR_Domains_Ultrametric_Tree_Plot:
    input:
        PATHd8 = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_Protein_Domains_filtered_{domain}_Seq_Tree_Time/{prefix}_{centype}_TE_{domain}_Lineages_PATHd8.Age.txt",
        domains = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_intact_LTR_Protein_Domains_filtered_{domain}_Seq/domains_counts.txt",
        Copia_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][0]),
        Copia_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Copia_R10_{}.fasta".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][1]),
        Gypsy_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][0]),
        Gypsy_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Tree/Gypsy_R10_{}.fasta".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][1])
    output:
        Copia_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Ultrametric_Tree/Copia_{}_outgroup.txt".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][0]),
        Copia_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Ultrametric_Tree/Copia_{}_outgroup.txt".format(config["Centromere_Analysis"]["copia_lineages_outgroup"][1]),
        Gypsy_R10_og1 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Ultrametric_Tree/Gypsy_{}_outgroup.txt".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][0]),
        Gypsy_R10_og2 = "results/{{prefix}}/Centromere_Analysis/TEs_Distribution/{{prefix}}_{{centype}}/DANTE/{{prefix}}_{{centype}}_Protein_Domains_filtered_{{domain}}_Seq_Ultrametric_Tree/Gypsy_{}_outgroup.txt".format(config["Centromere_Analysis"]["gypsy_lineages_outgroup"][1])
    params:
        workdir = config["workdir"],
        species = "{prefix}",
        copia_lineages_outgroup = config["Centromere_Analysis"]["copia_lineages_outgroup"],
        gypsy_lineages_outgroup = config["Centromere_Analysis"]["gypsy_lineages_outgroup"],
        indir = lambda w, input: os.path.dirname(input["PATHd8"]),
        outdir = lambda w, output: os.path.dirname(output[0])
    shell:
        '''
        set +e
        module load R/3.6.0
        copia=`cat {input.domains} | \\grep -i "LTR" | awk -F"|" 'BEGIN{{OFS="\\t"}} {{print $3,$NF}}' | cut -f 2 -d / | cut -f 1 -d : | \\grep -i "copia" | cut -f 2 | sort -V | xargs`
        gypsy=`cat {input.domains} | \\grep -i "LTR" | awk -F"|" 'BEGIN{{OFS="\\t"}} {{print $3,$NF}}' | cut -f 2 -d / | cut -f 1 -d : | \\grep -i "gypsy" | cut -f 2 | sort -V | xargs`
        cd {params.outdir}
        cat {params.workdir}/{input.Copia_R10_og1} | \\grep ">" | sed "s/>//g" | sed "s/:/_/g" | sed "s/\\[/_/g" | sed "s/\\]/_/g" | sed "s/|/_/g" | sed "s|/|_|g" | sed "1i id" > {params.workdir}/{output.Copia_R10_og1}
        cat {params.workdir}/{input.Copia_R10_og2} | \\grep ">" | sed "s/>//g" | sed "s/:/_/g" | sed "s/\\[/_/g" | sed "s/\\]/_/g" | sed "s/|/_/g" | sed "s|/|_|g" |  sed "1i id" > {params.workdir}/{output.Copia_R10_og2}
        cat {params.workdir}/{input.Gypsy_R10_og1} | \\grep ">" | sed "s/>//g" | sed "s/:/_/g" | sed "s/\\[/_/g" | sed "s/\\]/_/g" | sed "s/|/_/g" | sed "s|/|_|g" |  sed "1i id" > {params.workdir}/{output.Gypsy_R10_og1}
        cat {params.workdir}/{input.Gypsy_R10_og2} | \\grep ">" | sed "s/>//g" | sed "s/:/_/g" | sed "s/\\[/_/g" | sed "s/\\]/_/g" | sed "s/|/_/g" | sed "s|/|_|g" |  sed "1i id" > {params.workdir}/{output.Gypsy_R10_og2}

        ##################################################################
        ################################ Ty1/copia #######################
        ##################################################################
        if [ -z "${{copia}}" ]
        then
            echo "##################################################################"
            echo "##################### Ty1/Copia is empty #########################"
            echo "##################################################################"
        else
            echo "##################################################################"
            echo "############### Ty1/Copia is NOT empty and Run It ################"
            echo "##################################################################"
            cat {params.workdir}/{params.indir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.txt | \\grep "d8 tree" | sed "s/d8 tree    : //g" > {params.workdir}/{params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.tre
            Rscript {params.workdir}/workflow/scripts/TreePlot.R {params.workdir}/{params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.tre {params.workdir}/{output.Copia_R10_og2} {params.workdir}/{params.outdir}/Copia_{params.copia_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.pdf

            for i in `echo ${{copia}} | sed "s/ /\\n/g" | grep -v -i "{params.copia_lineages_outgroup[0]}" | sort -V | xargs`
            do
                echo "##################################################################"
                echo "##################### Copia: ${{i}} #############################"
                echo "##################################################################"
                cat {params.workdir}/{params.indir}/Copia_${{i}}_mafft_aln.fasta_PATHd8.txt | \\grep "d8 tree" | sed "s/d8 tree    : //g" > Copia_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.tre
                Rscript {params.workdir}/workflow/scripts/TreePlot.R Copia_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.tre {params.workdir}/{output.Copia_R10_og1} Copia_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.pdf
            done
        fi

        ##################################################################
        ################################ Ty3/gypsy #######################
        ##################################################################
        if [ -z "${{gypsy}}" ]
        then
            echo "##################################################################"
            echo "##################### Ty3/Gypsy is empty #########################"
            echo "##################################################################"
        else
            echo "##################################################################"
            echo "############### Ty3/Gypsy is NOT empty and Run It ################"
            echo "##################################################################"
            if [ -s {params.workdir}/{output.Gypsy_R10_og1} ]; then
                # The file is not-empty.
                cat {params.workdir}/{params.indir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.txt | \
                    \\grep "d8 tree" | sed "s/d8 tree    : //g" > {params.workdir}/{params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.tre
                if [ -s {params.workdir}/{output.Gypsy_R10_og2} ]; then
                    # The file is not-empty.
                    Rscript {params.workdir}/workflow/scripts/TreePlot.R \
                        {params.workdir}/{params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.tre TRUE \
                        {params.workdir}/{output.Gypsy_R10_og2} {params.workdir}/{params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.pdf
                else
                    Rscript {params.workdir}/workflow/scripts/TreePlot.R \
                        {params.workdir}/{params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.tre FALSE \
                        {params.workdir}/{params.outdir}/Gypsy_{params.gypsy_lineages_outgroup[0]}_mafft_aln.fasta_PATHd8.ultrametric.pdf
                fi
            else
                echo "The {params.workdir}/{output.Gypsy_R10_og1} is empty."
            fi

            for i in `echo ${{gypsy}} | sed "s/ /\\n/g" | grep -v -i "{params.gypsy_lineages_outgroup[0]}" | sort -V | xargs`
            do
                echo "##################################################################"
                echo "##################### Gypsy: ${{i}} #############################"
                echo "##################################################################"
                cat {params.workdir}/{params.indir}/Gypsy_${{i}}_mafft_aln.fasta_PATHd8.txt | \
                  \\grep "d8 tree" | sed "s/d8 tree    : //g" > Gypsy_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.tre
                if [ -s {params.workdir}/{output.Gypsy_R10_og1} ]; then
                    # The file is not-empty.
                    Rscript {params.workdir}/workflow/scripts/TreePlot.R Gypsy_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.tre TRUE \
                      {params.workdir}/{output.Gypsy_R10_og1} Gypsy_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.pdf
                else
                    # The file is empty.
                    Rscript {params.workdir}/workflow/scripts/TreePlot.R Gypsy_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.tre FALSE \
                      Gypsy_${{i}}_mafft_aln.fasta_PATHd8.ultrametric.pdf
                fi
            done
        fi
        '''

rule Centromere_Analysis_26_centromere_TEs_LTR_Domains_Age_Plot:
    input:
        "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_Protein_Domains_filtered_{domain}_Seq_Tree_Time/{prefix}_{centype}_TE_{domain}_Lineages_PATHd8.Age.txt"
    output:
        "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_Protein_Domains_filtered_{domain}_Seq_Tree_Time/{prefix}_{centype}_TE_{domain}_Lineages_PATHd8.Age.pdf"
    shell:
        '''
        module load R/3.6.0
        Rscript workflow/scripts/TE_Lineages_Age_Plot.R {input} {output}
        '''

rule Centromere_Analysis_27_centromere_TEs_LTR_ORF_Age_Plot:
    input:
        expand("results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/DANTE/{prefix}_{centype}_Protein_Domains_filtered_{domain}_Seq_Tree_Time/{prefix}_{centype}_TE_{domain}_Lineages_PATHd8.Age.txt", prefix = config["prefix"], domain = config["Centromere_Analysis"]["LTR_classification_domain"], centype = config["Centromere_Analysis"]["LTR_classification_domain_tree_centype"])
    output:
        RT_INT = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_Protein_Domains_filtered_ORF_Seq_Tree_Time_PATHd8.Age.txt",
        RT_INT_Plot = "results/{prefix}/Centromere_Analysis/TEs_Distribution/{prefix}_{centype}/{prefix}_Protein_Domains_filtered_ORF_Seq_Tree_Time_PATHd8.Age.pdf"
    shell:
        '''
        module load R/3.6.0
        cat {input[0]} | sed "1d" | cat {input[1]} - > {output.RT_INT}
        Rscript workflow/scripts/TE_Lineages_Age_Plot.R {output.RT_INT} {output.RT_INT_Plot}
        '''
