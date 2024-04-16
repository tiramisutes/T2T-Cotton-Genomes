############################################################################################################################################################
############################################################################################################################################################
rule GMM_1_genmap_mappability_index:
    input:
        genome = "results/{}/Genome_Mappability/Genomes/{}".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        directory("results/{}/Genome_Mappability/GenMap/genmap_{}_idx".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]))
    conda:
        "../envs/Genome_Mappability.yaml"
    shell:
        """
        genmap index -F {input} -I {output}
        """

rule GMM_2_genmap_mappability_map:
    input:
        "results/{}/Genome_Mappability/GenMap/genmap_{}_idx".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}.bedgraph".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    conda:
        "../envs/Genome_Mappability.yaml"
    params:
        kmer = "{kmer}",
        mismatches = "{mismatches}",
        outpre = lambda w, output: ".".join(output[0].split(".")[0:-1])
    shell:
        """
        genmap map -I {input} -O {params.outpre} -K {params.kmer} -E {params.mismatches} -w -bg
        """

rule GMM_3_genmap_mappability_plot_prepare_size:
    input:
        genome = "results/{}/Genome_Mappability/Genomes/{}".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        genome_size = "results/{}/Genome_Mappability/GenMap/{}_genome_size.txt".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    shell:
        """
        bioawk -c fastx '{{ print $name, "1", length($seq) }}' < {input.genome} | sort -k1,1V | sed '1i chr\\tstart\\tend' > {output.genome_size}
        """

rule GMM_4_genmap_mappability_plot_prepare_windows:
    input:
        genome = "results/{}/Genome_Mappability/Genomes/{}".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        windows_bed = "results/{}/Genome_Mappability/GenMap/{}_genome_w{{windows}}.bed".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    conda:
        "../envs/Genome_Mappability.yaml"
    params:
        windows = "{windows}"
    shell:
        """
        module load BEDTools/2.27
        samtools faidx {input.genome}
        bedtools makewindows -g {input.genome}.fai -w {params.windows} > {output.windows_bed}
        """

rule GMM_5_genmap_mappability_plot:
    input:
        bedgraph = "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}.bedgraph".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        genome_size = "results/{}/Genome_Mappability/GenMap/{}_genome_size.txt".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        windows_bed = "results/{}/Genome_Mappability/GenMap/{}_genome_w{{windows}}.bed".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        windows_bed_mean = "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}_w{{windows}}.bedgraph.mean".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        plot = "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}_w{{windows}}.pdf".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        plot2_pdf = "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}_w{{windows}}_plot.pdf".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        #plot2_png = "results/{}/Genome_Mappability/GenMap/genmap_{}_k{{kmer}}_e{{mismatches}}_w{{windows}}_plot.png".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    conda:
        "../envs/Genome_Mappability.yaml"
    params:
        windows = "{windows}"
    shell:
        """
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        module load BEDTools/2.27
        bedtools map -a {input.windows_bed} -b {input.bedgraph} -c 4 -o mean > {output.windows_bed_mean}
        Rscript workflow/scripts/genmap_mappability_plot.R {input.genome_size} {output.windows_bed_mean} {output.plot}
        Rscript workflow/scripts/genmap_mappability_plot2.R {output.windows_bed_mean} {output.plot2_pdf}
        """

rule GMM_6_genmap_mappability_plot_onepage:
    input:
        windows_bed_mean = expand("results/{prefix}/Genome_Mappability/GenMap/genmap_{genome}_k{kmer}_e{mismatches}_w{windows}.bedgraph.mean", prefix = config["prefix"], genome = config["Genome_Mappability"]["T2T_MUK_GENOMES"], kmer = config["Genome_Mappability"]["genmap_mappability_kmer"], mismatches = config["Genome_Mappability"]["genmap_mappability_mismatches"], allow_missing=True),
        genome_size = "results/{}/Genome_Mappability/GenMap/{}_genome_size.txt".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    output:
        windows_bed_mean_all = "results/{}/Genome_Mappability/GenMap/genmap_{}_allKmer_w{{windows}}.mean.bedgraph".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        plot = "results/{}/Genome_Mappability/GenMap/genmap_{}_allKmer_w{{windows}}.pdf".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        plot2_pdf = "results/{}/Genome_Mappability/GenMap/genmap_{}_allKmer_w{{windows}}_plot.pdf".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"]),
        #plot2_png = "results/{}/Genome_Mappability/GenMap/genmap_{}_allKmer_w{{windows}}_plot.png".format(config["prefix"], config["Genome_Mappability"]["T2T_MUK_GENOMES"])
    conda:
        "../envs/Genome_Mappability.yaml"
    shell:
        """
        module load R/3.6.0
        module load GCC/7.2.0-2.29
        module load BEDTools/2.27
        for i in {input.windows_bed_mean}
        do
        cat ${{i}} | awk -F"\\t" -v file=${{i##*/}} 'BEGIN{{OFS="\\t"}} {{print $0,file}}' >> {output.windows_bed_mean_all}
        done
        Rscript workflow/scripts/genmap_mappability_plot_onepage.R {input.genome_size} {output.plot} {input.windows_bed_mean}
        Rscript workflow/scripts/genmap_mappability_plot2_onepage.R {output.windows_bed_mean_all} {output.plot2_pdf}
        """

rule GMM_7_genmap_mappability_plot_list:
    input:
        expand("results/{prefix}/Genome_Mappability/GenMap/genmap_{genome}_k{kmer}_e{mismatches}.bedgraph", prefix = config["prefix"], genome = config["Genome_Mappability"]["T2T_MUK_GENOMES"], kmer = config["Genome_Mappability"]["genmap_mappability_kmer"], mismatches = config["Genome_Mappability"]["genmap_mappability_mismatches"]),
        expand("results/{prefix}/Genome_Mappability/GenMap/genmap_{genome}_k{kmer}_e{mismatches}_w{windows}.pdf", prefix = config["prefix"], genome = config["Genome_Mappability"]["T2T_MUK_GENOMES"], kmer = config["Genome_Mappability"]["genmap_mappability_kmer"], mismatches = config["Genome_Mappability"]["genmap_mappability_mismatches"], windows = config["Genome_Mappability"]["genmap_mappability_plot_windows"]),
        expand("results/{prefix}/Genome_Mappability/GenMap/genmap_{genome}_allKmer_w{windows}.pdf", prefix = config["prefix"], genome = config["Genome_Mappability"]["T2T_MUK_GENOMES"], windows = config["Genome_Mappability"]["genmap_mappability_plot_windows"])
