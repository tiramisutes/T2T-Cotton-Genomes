def select_BioNano_input_genome():
    BioNano_input_genome = "results/{}/BioNano/BioNano_genome.fasta".format(config["prefix"])
    print("\033[31m[BioNano]:\033[0m \033[0;36mThe input genome is " + BioNano_input_genome + "\033[0m")
    return BioNano_input_genome

def input_genome_name():
    input_file = select_BioNano_input_genome()
    name = os.path.basename(input_file)
    return name

rule BioNano_1_bnx_stats:
    input:
        bnx = "results/{prefix}/Resources/raw_bionano_list.csv"
    output:
        report("results/{prefix}/BioNano/bnx_stats/{prefix}_MapStatsHistograms.pdf", caption="../report/BNX_stat.rst", category="5. BioNano stats")
    threads: 1
    conda:
        "../envs/BioNano.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/{prefix}/BioNano/bnx_stats.log"
    shell:
        '''
        mkdir -p {params.outdir}
        cd {params.outdir}
        bnx=`cat {params.workdir}/{input.bnx} | xargs`
        /bin/perl {params.workdir}/workflow/bin/bnx_stats.pl --min_length_kb 150 ${{bnx}} > {params.workdir}/{log} 2>&1
        mv {params.workdir}/{params.outdir}/MapStatsHistograms.pdf {params.workdir}/{output}
        '''

rule BioNano_2_runBNG_QC:
    input:
        bnx = "results/{prefix}/Resources/raw_bionano_list.csv"
    output:
        report("results/{prefix}/BioNano/runBNG/{prefix}_BNX_stat.list", caption="../report/BNX_stat.rst", category="5. BioNano stats")
    threads: 1
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/{prefix}/BioNano/runBNG.log"
    shell:
        '''
        export PATH="/bin:$PATH"
        mkdir -p {params.outdir}
        for i in `cat {input.bnx} | xargs`
        do
        bnx_name=`basename ${{i}}`
        bash {params.workdir}/workflow/bin/runBNG/runBNG bnxstats -b ${{i}} -p ${{bnx_name}} -o {params.outdir} >> {log} 2>&1
        done
        ls {params.outdir} > {output}
        '''

rule BioNano_3_bnx_filter_molecules_labels:
    input:
        "results/{prefix}/Resources/raw_bionano_list.csv"
    output:
        filterbnx_ml = "results/{prefix}/BioNano/BNX_filter/filter_molecules_labels.bnx",
        filterbnx_SNR = "results/{prefix}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx"
    threads: 1
    params:
        workdir = config["workdir"],
        threshold = config["BioNano"]["bnx_filter_threshold"],
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/{prefix}/BioNano/bnx_filter_molecules_labels.log"
    shell:
        '''
        mkdir -p {params.workdir}/{params.outdir}
        filterbnx_ml_name=`basename {output.filterbnx_ml}`
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel
        ./RefAligner -if {params.workdir}/{input} -merge -bnx -o {params.workdir}/{params.outdir}/${{filterbnx_ml_name%.*}} {params.threshold} >> {params.workdir}/{log} 2>&1
        
        export PATH="$PATH:{params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019
        perl filter_SNR_dynamic.pl -i {params.workdir}/{output.filterbnx_ml} -o {params.workdir}/{output.filterbnx_SNR} >> {params.workdir}/{log} 2>&1
        '''

rule BioNano_4_bnx_filter_mapping_genome:
    input:
        select_BioNano_input_genome(),
        filterbnx_SNR = "results/{}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx".format(config["prefix"])
    output:
        rcmap = "results/{}/BioNano/BNX_filter/{}_{}_{}kb_{}labels.cmap".format(config["prefix"], os.path.basename(select_BioNano_input_genome()).split(".")[0], config["BioNano"]["enzyme"], config["BioNano"]["fa2cmap"]["min_molecule_length"], config["BioNano"]["fa2cmap"]["min_number_enzymes_molecule"]),
        rcmap_summary = report("results/{}/BioNano/BNX_filter/{}_{}_{}kb_{}labels_summary.txt".format(config["prefix"], os.path.basename(select_BioNano_input_genome()).split(".")[0], config["BioNano"]["enzyme"], config["BioNano"]["fa2cmap"]["min_molecule_length"], config["BioNano"]["fa2cmap"]["min_number_enzymes_molecule"]), caption="../report/BNX_stat.rst", category="5. BioNano stats")
    threads: 20
    conda:
        "../envs/BioNano.yaml"
    params:
        outdir = lambda w, output: os.path.dirname(output['rcmap']),
        workdir = config["workdir"],
        prefix = config["prefix"],
        enzyme = config["BioNano"]["enzyme"],
        min_molecule_length = config["BioNano"]["fa2cmap"]["min_molecule_length"],
        min_number_enzymes_molecule = config["BioNano"]["fa2cmap"]["min_number_enzymes_molecule"]
    log:
        "logs/{}/BioNano/bnx_filter_mapping_genome.log".format(config["prefix"])
    shell:
        '''
        echo "#########################################################################################"
        echo "1. fa to cmap"
        echo "#########################################################################################"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/HybridScaffold/12162019/scripts
        perl fa2cmap_multi_color.pl -i {params.workdir}/{input[0]} -e {params.enzyme} 1 -o {params.workdir}/{params.outdir} -m {params.min_number_enzymes_molecule} -M {params.min_molecule_length} >> {params.workdir}/{log} 2>&1

        echo "##################### Stats of the generated cmap file ############################" >> {params.workdir}/{output.rcmap_summary}
        perl calc_cmap_stats.pl {params.workdir}/{output.rcmap} >> {params.workdir}/{output.rcmap_summary}
        
        density=`grep -v "^#" {params.workdir}/{output.rcmap} | cut -f1-3 | uniq | awk '{{sum+=$2;SUM+=$3}}END{{print SUM/sum*100000}}'`
        echo "#################################################" >> {params.workdir}/{output.rcmap_summary}
        echo "Label density (/100Kb) = `printf "%.3f" $density`" >> {params.workdir}/{output.rcmap_summary}
        '''

rule BioNano_5_filter_bnx_DeNovo_Without_reference:
    input:
        filterbnx_mg = "results/{}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx".format(config["prefix"])
    output:
        DeNovoCmap = "results/{}/BioNano/DeNovo/Without_reference/contigs/{}_No_Rcmap_refineFinal1/{}_NO_RCMAP_REFINEFINAL1.cmap".format(config["prefix"], config["prefix"], config["prefix"].upper()),
        errbin = "results/{}/BioNano/DeNovo/Without_reference/contigs/auto_noise/autoNoise1.errbin".format(config["prefix"])
    threads: 30
    conda:
        "../envs/BioNano.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        deNovo_optArguments_XML = config["BioNano"]["deNovo_optArguments_XML_DLE"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] == "CTTAAG" 
                             else config["BioNano"]["deNovo_optArguments_XML"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] != "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype_DLE"] if config["BioNano"]["assembly_type"] == "haplotype" and config["BioNano"]["enzyme"] == "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype"],
        wooutdir = lambda w, output: "/".join(os.path.dirname(output[0]).split("/")[0:5])
    log:
        "logs/{}/BioNano/DeNovo_Without_reference.log".format(config["prefix"])
    shell:
        '''
        export PATH="/bin:$PATH"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019
        
        echo "Without reference"
        mkdir -p {params.workdir}/{params.wooutdir}
        python pipelineCL.py -T {threads} -N 4 -i 5 -b {params.workdir}/{input.filterbnx_mg} \
         -l {params.workdir}/{params.wooutdir} \
         -t {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel \
         -e {params.prefix}_No_Rcmap -c 0 \
         -a {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/{params.deNovo_optArguments_XML} \
         -j 10 -je 10 -jp 10 -J 10 -TJ 10 -Te 10 \
         -R > {params.workdir}/{log} 2>&1
        '''

rule BioNano_6_Hybrid_Scaffold_xml:
    output:
        "results/{}/BioNano/{}_hybridScaffold_config.xml".format(config["prefix"], config["prefix"])
    threads: 1
    params:
        workdir = config["workdir"],
        hybrid_Scaffold_xml = config["BioNano"]["hybrid_Scaffold_xml_DLE"] if config["BioNano"]["enzyme"] == "CTTAAG" else config["BioNano"]["hybrid_Scaffold_xml"],
        enzyme = config["BioNano"]["enzyme"]
    shell:
        '''
        sed "s/attr=\\"enzyme\\" val0=\\"CTTAAG\\"/attr=\\"enzyme\\" val0=\\"{params.enzyme}\\"/g" {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/HybridScaffold/12162019/{params.hybrid_Scaffold_xml} > {output}
        '''

rule BioNano_7_Hybrid_Scaffold:
    input:
        select_BioNano_input_genome(),
        filterbnx_mg = "results/{}/BioNano/BNX_filter/filter_molecules_labels_SNR.bnx".format(config["prefix"]),
        DeNovoCmap = "results/{}/BioNano/DeNovo/Without_reference/contigs/{}_No_Rcmap_refineFinal1/{}_NO_RCMAP_REFINEFINAL1.cmap".format(config["prefix"], config["prefix"], config["prefix"].upper()),
        errbin = "results/{}/BioNano/DeNovo/Without_reference/contigs/auto_noise/autoNoise1.errbin".format(config["prefix"]),
        hybridScaffold_config = "results/{}/BioNano/{}_hybridScaffold_config.xml".format(config["prefix"], config["prefix"])
    output:
        bg = report("results/{}/BioNano/Hybrid_Scaffold/hybrid_scaffolds/{}_NO_RCMAP_REFINEFINAL1_bppAdjust_cmap_{}_fa_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta".format(config["prefix"], config["prefix"].upper(), os.path.basename(select_BioNano_input_genome()).split(".")[0]), caption="../report/BioNano_stat.rst", category="3. BioNano"),
        bgn = report("results/{}/BioNano/Hybrid_Scaffold/hybrid_scaffolds/{}_NO_RCMAP_REFINEFINAL1_bppAdjust_cmap_{}_fa_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta".format(config["prefix"], config["prefix"].upper(), os.path.basename(select_BioNano_input_genome()).split(".")[0]), caption="../report/BioNano_stat.rst", category="3. BioNano"),
        bg_report = report("results/{}/BioNano/Hybrid_Scaffold/hybrid_scaffolds/hybrid_scaffold_informatics_report.txt".format(config["prefix"]), caption="../report/BioNano_stat.rst", category="3. BioNano")
    threads: 30
    conda:
        "../envs/BioNano.yaml"
    params:
        workdir = config["workdir"],
        prefix = config["prefix"],
        conflict_cut = config["BioNano"]["conflict_cut"],
        deNovo_optArguments_XML = config["BioNano"]["deNovo_optArguments_XML_DLE"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] == "CTTAAG" 
                             else config["BioNano"]["deNovo_optArguments_XML"] if config["BioNano"]["assembly_type"] == "nonhaplotype" and config["BioNano"]["enzyme"] != "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype_DLE"] if config["BioNano"]["assembly_type"] == "haplotype" and config["BioNano"]["enzyme"] == "CTTAAG"
                             else config["BioNano"]["deNovo_optArguments_XML_haplotype"],
        outdir = lambda w, output: "/".join(os.path.dirname(output[0]).split("/")[0:4])
    log:
        "logs/{}/BioNano/Hybrid_Scaffold.log".format(config["prefix"])
    shell:
        '''
        export PATH="/bin:$PATH"
        echo "Single-Enzyme Hybrid Scaffold Pipeline"
        cd {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/HybridScaffold/12162019
        
        mkdir -p {params.workdir}/{params.outdir}
        perl hybridScaffold.pl -n {params.workdir}/{input[0]} -b {params.workdir}/{input.DeNovoCmap} \
         -c {params.workdir}/{input.hybridScaffold_config} \
         -r {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/RefAligner \
         -o {params.workdir}/{params.outdir} {params.conflict_cut} -f -x -y -m {params.workdir}/{input.filterbnx_mg} \
         -p {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/Pipeline/12162019 \
         -q {params.workdir}/workflow/bin/Bionano_Tools_v1.5.3/pipeline/Solve3.5.1_01142020/RefAligner/10330.10436rel/{params.deNovo_optArguments_XML} \
         -e {params.workdir}/{input.errbin} > {params.workdir}/{log} 2>&1
        '''

rule BioNano_8_assembly_genome:
    input:
        bg = "results/{}/BioNano/Hybrid_Scaffold/hybrid_scaffolds/{}_NO_RCMAP_REFINEFINAL1_bppAdjust_cmap_{}_fa_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta".format(config["prefix"], config["prefix"].upper(), os.path.basename(select_BioNano_input_genome()).split(".")[0]),
        bgn = "results/{}/BioNano/Hybrid_Scaffold/hybrid_scaffolds/{}_NO_RCMAP_REFINEFINAL1_bppAdjust_cmap_{}_fa_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta".format(config["prefix"], config["prefix"].upper(), os.path.basename(select_BioNano_input_genome()).split(".")[0])
    output:
        bg1 = "results/{}/BioNano/BioNano_genome_org.fasta".format(config["prefix"]),
        bg2 = "results/{}/BioNano/BioNano_genome.fasta".format(config["prefix"])
    shell:
        '''
        cp {input.bg} {output.bg1}
        cat {input.bg} {input.bgn} > {output.bg2}
        '''

rule BioNano_9_assembly_stat:
    input:
        bg1 = "results/{}/BioNano/BioNano_genome_org.fasta".format(config["prefix"]),
        bg2 = "results/{}/BioNano/BioNano_genome.fasta".format(config["prefix"])
    output:
        bg1 = report("results/{}/BioNano/BioNano_genome_org.fasta.stat".format(config["prefix"]), caption="../report/BioNano_stat.rst", category="3. BioNano"),
        bg2 = report("results/{}/BioNano/BioNano_genome.fasta.stat".format(config["prefix"]), caption="../report/BioNano_stat.rst", category="3. BioNano")
    params:
        workdir = config["workdir"]
    shell:
        '''
        cd {params.workdir}/workflow/bin
        ./assemblathon_stats.pl {params.workdir}/{input.bg1} > {params.workdir}/{output.bg1}
        ./assemblathon_stats.pl {params.workdir}/{input.bg2} > {params.workdir}/{output.bg2}
        '''

rule BioNano_10_assembly_BUSCO2:
    input:
        "results/{}/BioNano/BioNano_genome.fasta".format(config["prefix"])
    output:
        summary = report("results/{}/BioNano/BUSCO/BioNano_genome.fasta/short_summary.specific.{}.BioNano_genome.fasta.txt".format(config["prefix"], config["BUSCO_databsse"]), caption="../report/BUSCO.rst", category="9. Assessment_Finally_Assembly", subcategory="BUSCO"),
        figure = report("results/{}/BioNano/BUSCO/BioNano_genome.fasta/busco_figure.png".format(config["prefix"]), caption="../report/BUSCO.rst", category="9. Assessment_Finally_Assembly", subcategory="BUSCO")
    threads: 10
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    log:
        "logs/{}/BioNano/BioNano_8_assembly_BUSCO.log".format(config["prefix"])
    params:
        workdir = config["workdir"],
        prefix = lambda w, input: os.path.basename(input[0]),
        database = config["BUSCO_databsse"],
        outdir = lambda w, output: "/".join(os.path.dirname(output['summary']).split("/")[0:4])
    shell:
        '''
        cd {params.outdir}
        busco -i {params.workdir}/{input} --out_path {params.workdir}/{params.outdir} -o {params.prefix} -l {params.database} -f -m genome --augustus -c {threads} --long > {params.workdir}/{log} 2>&1
        generate_plot.py -wd {params.workdir}/{params.outdir}/{params.prefix}
        '''

rule BioNano_11_assembly_BUSCO:
    input:
        "results/{}/BioNano/BioNano_genome_org.fasta".format(config["prefix"])
    output:
        summary = report("results/{}/BioNano/BUSCO/BioNano_genome_org.fasta/short_summary.specific.{}.BioNano_genome_org.fasta.txt".format(config["prefix"], config["BUSCO_databsse"]), caption="../report/BUSCO.rst", category="9. Assessment_Finally_Assembly", subcategory="BUSCO"),
        figure = report("results/{}/BioNano/BUSCO/BioNano_genome_org.fasta/busco_figure.png".format(config["prefix"]), caption="../report/BUSCO.rst", category="9. Assessment_Finally_Assembly", subcategory="BUSCO")
    threads: 10
    conda:
        "../envs/Assessment_Finally_Assembly.yaml"
    log:
        "logs/{}/BioNano/BioNano_8_assembly_BUSCO2.log".format(config["prefix"])
    params:
        workdir = config["workdir"],
        prefix = lambda w, input: os.path.basename(input[0]),
        database = config["BUSCO_databsse"],
        outdir = lambda w, output: "/".join(os.path.dirname(output['summary']).split("/")[0:4])
    shell:
        '''
        cd {params.outdir}
        busco -i {params.workdir}/{input} --out_path {params.workdir}/{params.outdir} -o {params.prefix} -l {params.database} -f -m genome --augustus -c {threads} --long > {params.workdir}/{log} 2>&1
        generate_plot.py -wd {params.workdir}/{params.outdir}/{params.prefix}
        '''

rule BioNano_12_genome_finally_mapping:
    input:
        select_BioNano_input_genome(),
        "results/{}/BioNano/BioNano_genome.fasta".format(config["prefix"])
    output:
        "results/{}/BioNano/NGSgenome2Bionano_genome.paf".format(config["prefix"])
    threads: 20
    shell:
        '''
        minimap2 -x asm5 -t {threads} {input[1]} {input[0]} > {output}
        '''

rule BioNano_13_genome_finally_mapping_plot:
    input:
        "results/{}/BioNano/NGSgenome2Bionano_genome.paf".format(config["prefix"])
    output:
        directory("results/{}/BioNano/Check_Synteny".format(config["prefix"]))
    threads: 1
    params:
        chromosome = config["params"]["chromosome"]
    shell:
        '''
        module load R/3.6.0
        for i in {params.chromosome}
        do
        Rscript workflow/scripts/paf_plot_chr.R {input} ${{i}} ${{i}} {output}
        done
        '''

rule BioNano_14_list:
    input:
        "results/{}/BioNano/bnx_stats/{}_MapStatsHistograms.pdf".format(config["prefix"], config["prefix"]),
        "results/{}/BioNano/runBNG/{}_BNX_stat.list".format(config["prefix"], config["prefix"]),
        "results/{}/BioNano/BNX_filter/{}_{}_{}kb_{}labels_summary.txt".format(config["prefix"], os.path.basename(select_BioNano_input_genome()).split(".")[0], config["BioNano"]["enzyme"], config["BioNano"]["fa2cmap"]["min_molecule_length"], config["BioNano"]["fa2cmap"]["min_number_enzymes_molecule"]),
        "results/{}/BioNano/BioNano_genome_org.fasta.stat".format(config["prefix"]),
        "results/{}/BioNano/BioNano_genome.fasta.stat".format(config["prefix"]),
        "results/{}/BioNano/BUSCO/BioNano_genome_org.fasta/short_summary.specific.{}.BioNano_genome_org.fasta.txt".format(config["prefix"], config["BUSCO_databsse"]),
        "results/{}/BioNano/BUSCO/BioNano_genome.fasta/short_summary.specific.{}.BioNano_genome.fasta.txt".format(config["prefix"], config["BUSCO_databsse"]),
        "results/{}/BioNano/Check_Synteny".format(config["prefix"])
