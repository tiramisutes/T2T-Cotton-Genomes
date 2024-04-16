setwd("F:/HZAU/CottonGenome/")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure2 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure2a #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library("karyoploteR")
library("rtracklayer")

custom.genome <- toGRanges("Jin668_genome_size.txt")
ont <- read.table("ont.depth", header =T)
hifi <- read.table("hifi.depth", header =T)
ngs <- read.table("ngs.depth", header =T)
ont.bg <- toGRanges(ont)
hifi.bg <- toGRanges(hifi)
ngs.bg <- toGRanges(ngs)

## Area
pdf("Figure2a.pdf", width = 15, height = 20)
kp <- plotKaryotype(genome = custom.genome, plot.type = 2)
# kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")

kp <- kpArea(kp, data=ont.bg, r0=0, r1=0.27, y=ont.bg$Avg_depth, ymin=0, ymax=173, border="#0099e5", col="#0099e5", data.panel = 1)
kpAxis(kp, ymax=173, r0=0, r1=0.27, cex=0.8, data.panel = 1)

kp <- kpArea(kp, data=hifi.bg, r0=0.34, r1=0.61, y=hifi.bg$Avg_depth, ymin=0, ymax=25, border="#ff4c4c", col="#ff4c4c", data.panel = 1)
kpAxis(kp, ymax=25, r0=0.34, r1=0.61, cex=0.8, data.panel = 1)

kp <- kpArea(kp, data=ngs.bg, r0=0.68, r1=1, y=ngs.bg$Avg_depth, ymin=0, ymax=100, border="#34bf49", col="#34bf49", data.panel=1)
kpAxis(kp, ymax=100, r0=0.68, r1=1, cex=0.8, data.panel=1)

legend(x = "bottomright", fill = c("#0099e5", "#ff4c4c", "#34bf49"), legend = c(paste0("ONT_Windows_", ont[1,3]), paste0("HiFi_Windows_", hifi[1,3]), paste0("Illumina_Windows_", ngs[1,3])))
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure3 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure3a #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library("karyoploteR")
library(rtracklayer)

custom.genome <- toGRanges(argv[2])
cenh.bg <- import(argv[3], format="bedGraph")
data.bg50 <- import(argv[4], format="bedGraph")
data.bg100 <- import(argv[5], format="bedGraph")
data.bg150 <- import(argv[6], format="bedGraph")
data.bg200 <- import(argv[7], format="bedGraph")
data.bg300 <- import(argv[8], format="bedGraph")

pdf("Figure3a.pdf", width = 15, height = 8)
kp <- plotKaryotype(genome = custom.genome, plot.type = 2)


kpArea(kp, data=cenh.bg, data.panel=1, y=cenh.bg$score, ymin=0, ymax=max(cenh.bg$score), col="red", r0=0, r1=1)
kpArea(kp, data=data.bg50, data.panel=2, y=data.bg50$score, ymin=0, ymax=max(data.bg50$score), col="#ca5572", r0=0, r1=0.2)
kpArea(kp, data=data.bg100, data.panel=2, y=data.bg100$score, ymin=0, ymax=max(data.bg100$score), col="#72a553", r0=0.2, r1=0.4)
kpArea(kp, data=data.bg150, data.panel=2, y=data.bg150$score, ymin=0, ymax=max(data.bg150$score), col="#a265c2", r0=0.4, r1=0.6)
kpArea(kp, data=data.bg200, data.panel=2, y=data.bg200$score, ymin=0, ymax=max(data.bg200$score), col="#c57c3d", r0=0.6, r1=0.8)
kpArea(kp, data=data.bg300, data.panel=2, y=data.bg300$score, ymin=0, ymax=max(data.bg300$score), col="#6097ce", r0=0.8, r1=1)
## add legend
legend(x = "bottomright", fill = c("#ca5572", "#72a553", "#a265c2", "#c57c3d", "#6097ce"), legend = c("50 nucleotides", "100 nucleotides", "150 nucleotides", "200 nucleotides", "300 nucleotides"))
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure3e #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(ggplot2)
library(dplyr)
library(ggbreak) 
age_file <- "Jin668_centromere_TE_Protein_Domains_Age.txt"
age_plot_count <- "Jin668_centromere_TE_Protein_Domains_Age_count.pdf"
age_plot <- "Figure3e.pdf"

df <- read.table(age_file, header = T)
#df <- filter(df,Age<0.25)
df_count <- df %>% dplyr::count(Species,Type,Lineages)
df_count$Type <- factor(df_count$Type, levels = c("RT","INT"))
df_count_table <- dcast(df_count, Lineages ~ Species+Type)

pdf(age_plot_count, width = 4, heigh = 10)
ggplot(df_count, aes(x=Type, y=n, fill=Type)) +
  geom_bar(stat="identity", alpha=0.4) +
  facet_grid(rows = vars(Lineages)) + 
  scale_fill_brewer(palette="Set1") +
  geom_text(data=df_count, aes(x = Type, y = n, label = n), size = 4) +
  theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +
  theme(panel.background=element_rect(fill='transparent', color='black')) +
  theme_bw()
dev.off()

pdf(age_plot, width = 12, heigh = 9)
df$Type <- factor(df$Type, levels = c("RT","INT"))
ggplot(df, aes(x = Age, fill = Type)) +
  #geom_area(stat = "bin") +
  geom_density(alpha=0.4) +
  #facet_grid(rows = vars(Lineages), scales="free") +
  facet_grid(rows = vars(Lineages)) +
  scale_x_break(c(0.25, 0.95)) +
  scale_fill_brewer(palette="Set1") +
  theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +
  theme(panel.background=element_rect(fill='transparent', color='black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure3f #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(ggpubr)
library(dplyr)
library(plyr)
library(ggstatsplot)
library("ggthemes")

df <- read.table("Jin668_centromere_OutCentromere_Identity.txt", header = F)
ddply(df, c("V4","V5"), summarise, grp.mean=mean(V3))
ddply(df, c("V4","V5"), summarise, grp.median=median(V3))
df %>% dplyr::group_by(V4, V5) %>% dplyr::summarise(count = n())

pdf("Figure3f.pdf", width = 6, heigh = 6)
grouped_ggbetweenstats(
  data = df,
  x = V4,
  y = V3,
  grouping.var = V5,
  bf.message = FALSE
) + labs(x="", y = "LTR % Identity")

dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure4 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure4b #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library("karyoploteR")
library(rtracklayer)
library(tidyverse)

custom.genome <- toGRanges(argv[1])
custom.cytobands <- toGRanges(argv[2])
lengthf1 <- as.numeric(argv[3])
lengthf2 <- as.numeric(argv[4])
sddf <- read.table(argv[5], header =F) %>% mutate(length1 = V3 - V2, length2 = V6 - V5)

sddf$id <- row.names(sddf)
sddfg <- sddf %>% filter(length1 > lengthf1 & length2 > lengthf1)
sddfm <- sddf %>% filter(length1 > lengthf2 & length1 < lengthf1 | length2 > lengthf2 & length2 < lengthf1)
tmdf <- rbind(sddfg, sddfm)
sddfl <- sddf %>% filter(!(id %in% tmdf$id))
sddata1 <- toGRanges(sddfg[,c(1:3)])
sddata2 <- toGRanges(sddfg[,c(4:6)])
sddata3 <- toGRanges(sddfm[,c(1:3)])
sddata4 <- toGRanges(sddfm[,c(4:6)])
sddata5 <- toGRanges(sddfl[,c(1:3)])
sddata6 <- toGRanges(sddfl[,c(4:6)])
plot <- argv[6]
plotl <- argv[7]

pdf(plot, width = 10, height = 10)
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, plot.type = 2)
kpPlotLinks(kp, data=sddata1, data2=sddata2, data.panel = 1, col="#fac7ffaa", r0=0.1)
kpPlotLinks(kp, data=sddata3, data2=sddata4, data.panel = 2, col="#8d9aff", r0=0.1)
legend(x = "bottomright", fill = c("#fac7ffaa", "#8d9aff"), legend = c(paste0("> ", lengthf1), paste0("< ", lengthf1, " & > ", lengthf2)))
dev.off()

pdf(plotl, width = 10, height = 10)
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, plot.type = 1)
kpPlotLinks(kp, data=sddata5, data2=sddata6, data.panel = 1, col="#fac7ffaa", r0=0.1)
legend(x = "bottomright", fill = c("#fac7ffaa"), legend = paste0("< ", lengthf2))
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure5 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure5b #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(dplyr)
library(stringr)
library(ggpubr)
library(ggh4x)

df <- read.table("genome2genome_SV.txt", header = F)
colnames(df) <- c("chr", "start", "end", "species", "svtype")

genome_szie <- read.table("genome.size", header = F)
colnames(genome_szie) <- c("chr", "chr_length")
genome_szie_chr <- genome_szie %>% dplyr::filter(!grepl("Scaffol", chr))
genome_szie_chr_sum <- sum(genome_szie_chr$chr_length)

svtypelabs <- c("SNPs" = "SNPs", "InDels" = "InDels", "Inversions" = "Invs", "Translocations" = "Trans", "PAVs" = "PAVs")
title_color <- c("#b5d3b4", "#cab6d7")
title_color_list <- rep(title_color, each=13)
strip <- strip_themed(background_x = elem_list_rect(fill = title_color_list), 
                      background_y = element_rect(fill = "white"))


sv_chr_stat <- df %>% dplyr::group_by(chr, species, svtype) %>% dplyr::summarise(sv_chr_count = n()) %>% 
  mutate(sv_chr_count_scale = case_when(grepl("SNPs", svtype) ~ sv_chr_count/1000,
                                        grepl("InDels", svtype) ~ sv_chr_count/1000,
                                        !grepl("SNPs", svtype) ~ sv_chr_count,
                                        !grepl("InDels", svtype) ~ sv_chr_count))
sv_chr_stat$svtype <- factor(sv_chr_stat$svtype, levels=c("SNPs", "InDels", "Inversions", "Translocations", "PAVs"))
sv_chr_stat_plot <- ggplot(sv_chr_stat %>% dplyr::filter(!grepl("Scaffol", chr)), aes(x = chr, y = sv_chr_count_scale, fill = svtype)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(rows = vars(species), scales = "free", space = "free") +
  labs(x="", y = "Number of SVs") +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="top",
        strip.background=element_rect(fill="white"))

sv_species_stat <- df %>% dplyr::group_by(species, svtype) %>% dplyr::summarise(sv_species_count = n()) %>% 
  mutate(sv_species_count_scale = case_when(grepl("SNPs", svtype) ~ sv_species_count/1000,
                                            grepl("InDels", svtype) ~ sv_species_count/1000,
                                            !grepl("SNPs", svtype) ~ sv_species_count,
                                            !grepl("InDels", svtype) ~ sv_species_count))
sv_species_stat$svtype <- factor(sv_species_stat$svtype, levels=c("SNPs", "InDels", "Inversions", "Translocations", "PAVs"))
sv_species_stat$species <- factor(sv_species_stat$species, levels=c("ZM24", "YZ1", "Jin668"))
sv_species_stat_plot <- ggplot(sv_species_stat, aes(x = species, y = sv_species_count_scale, fill = svtype)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="", y = "Number of SVs") +
  scale_fill_npg() +
  theme_classic() + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="top",
        strip.background=element_rect(fill="white"))


sv_chr_length <- df %>% dplyr::mutate(sv_chr_length = case_when(grepl("SNPs", svtype) ~ as.integer(1),
                                                                !grepl("SNPs", svtype) ~ end - start)) %>% 
  ddply(c("chr", "species","svtype"), summarise, sv_chr_length_sum=sum(sv_chr_length)) %>% 
  dplyr::left_join(genome_szie, "chr") %>% dplyr::mutate(sv_chr_length_sum_pre = (sv_chr_length_sum/chr_length)*100)
sv_chr_length$svtype <- factor(sv_chr_length$svtype, levels=c("SNPs", "InDels", "Inversions", "Translocations", "PAVs"))

sv_chr_length_plot <- ggplot(sv_chr_length %>% dplyr::filter(!grepl("Scaffol", chr)), aes(x = chr, y = sv_chr_length_sum_pre, fill = svtype)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_grid(rows = vars(species), scales = "free", space = "free") +
  labs(x="", y = "Percentage of SVs in chromosome size (%)") +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="top",
        strip.background=element_rect(fill="white"))


sv_species_length <- sv_chr_length %>% dplyr::filter(!grepl("Scaffol", chr)) %>% 
  ddply(c("species","svtype"), summarise, sv_species_length_sum=sum(sv_chr_length_sum))

sv_species_length$genome_szie_sum <- genome_szie_chr_sum
sv_species_length <- sv_species_length %>% mutate(sv_species_length_sum_pre = (sv_species_length_sum/genome_szie_sum)*100)

sv_species_length$svtype <- factor(sv_species_length$svtype, levels=c("SNPs", "InDels", "Inversions", "Translocations", "PAVs"))
sv_species_length$species <- factor(sv_species_length$species, levels=c("ZM24", "YZ1", "Jin668"))
sv_species_length_plot <- ggplot(sv_species_length, aes(x = species, y = sv_species_length_sum_pre, fill = svtype)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="", y = "Percentage of SVs in genome size (%)") +
  scale_fill_npg() +
  theme_classic() + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="top",
        strip.background=element_rect(fill="white"))


sv_chr_length_distribution <- df %>% dplyr::filter(!grepl("Scaffol", chr)) %>% 
  mutate(sv_chr_length = case_when(grepl("SNPs", svtype) ~ as.integer(1),
                                   !grepl("SNPs", svtype) ~ end - start)) %>%
  dplyr::filter(svtype != "SNPs", svtype != "InDels")

sv_chr_length_distribution$svtype <- factor(sv_chr_length_distribution$svtype, levels=c("Inversions", "Translocations", "PAVs"))
sv_chr_length_distribution$species <- factor(sv_chr_length_distribution$species, levels=c("Jin668", "YZ1", "ZM24"))

sv_chr_length_distribution_plot <- ggplot(sv_chr_length_distribution, aes(x = sv_chr_length, color=svtype, fill=svtype)) + 
  geom_histogram(aes(y=..density..), bins = 30, position="identity")+
  #geom_density(alpha=.2) +
  #facet_grid(rows = vars(species+svtype) ,cols = vars(chr), scales = "free", space = "free")
  #facet_grid(species+svtype ~ chr, labeller = labeller(svtype = svtypelabs)) +
  ggh4x::facet_grid2(species+svtype ~ chr, labeller = labeller(svtype = svtypelabs), strip = strip) +
  labs(x="", y = "Density") + 
  #scale_x_continuous(breaks=seq(0,40000000,3)) + scale_y_continuous(breaks=seq(0,1.5e-06,3)) +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="top",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #strip.background=element_rect(fill="white")
  )



sv_chr_length_distribution_stable <- desc_statby(sv_chr_length_distribution, measure.var = "sv_chr_length", 
                                                 grps = c("species", "svtype")) %>%
  select(c("species", "svtype", "min", "max", "median", "mean", "sd")) %>% 
  dplyr::arrange(species, svtype) %>% 
  dplyr::mutate(across('svtype', str_replace, 'ersion', '')) %>% 
  dplyr::mutate(across('svtype', str_replace, 'locations', ''))

sv_chr_length_distribution_stable_plot <- ggtexttable(sv_chr_length_distribution_stable, rows = NULL, 
                                                      height = rep(unit(12, "mm") ,nrow(sv_chr_length_distribution_stable)),
                                                      theme = ttheme(
                                                        colnames.style = colnames_style(color = "white", fill = "#8cc257"),
                                                        tbody.style = tbody_style(
                                                          color = "black", fill = c("#e8f3de", "#d3e8bb"),
                                                          hjust = as.vector(matrix(c(0, 1, 1, 1, 1, 1, 1), ncol = 7, nrow = nrow(sv_chr_length_distribution_stable), byrow = TRUE)),
                                                          x = as.vector(matrix(c(.1, .9, .9,.9, .9, .9, .9), ncol = 7, nrow = nrow(sv_chr_length_distribution_stable), byrow = TRUE))
                                                        )
                                                      ))


pdf("Figure5b.pdf", width = 18, height = 10)
ggpubr::ggarrange(sv_chr_stat_plot, sv_species_stat_plot, 
                  sv_chr_length_plot, sv_species_length_plot,
                  sv_chr_length_distribution_plot, sv_chr_length_distribution_stable_plot,
                  ncol = 2, nrow = 3, 
                  widths = c(2, 1), heights = c(1, 1, 2),
                  common.legend = TRUE, legend = "top")
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure5d #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

set.seed(2023)

df <- read.table("Jin668_TPM_mean.tsv", header = T) %>% 
  dplyr::select(-c("Jin668")) %>% 
  tibble::rownames_to_column("index") %>% 
  tidyr::unite(genes, index, TM1, sep = "", remove = T) %>% 
  tibble::column_to_rownames(var="genes")

dfm <- df %>% dplyr::select(-c("type","class")) %>% dplyr::filter_all(any_vars(. != 0))  
annot <- df %>% dplyr::select(c("type","class"))
annots <- annot %>% dplyr::filter(row.names(annot) %in% row.names(dfm))

scaled_dfm <- t(scale(t(dfm), center = T))
scaled_z_score_row_df <- t(apply(dfm, 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
z_score_col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))


pdf("Figure5d.pdf", width = 10, height = 10)
Heatmap(scaled_z_score_row_df,
        col = z_score_col_fun,
        heatmap_legend_param = list(title = "Z-score"),
        cluster_rows = F, 
        cluster_columns = FALSE, 
        show_row_dend = FALSE, show_column_dend = FALSE,
        right_annotation = HeatmapAnnotation(type = annots$type, class = annots$class, which = "row"),
        row_split = factor(annots$class),
        column_split = rep(c("Jin668", "TM-1"), each = 5),
        show_row_names = F, 
        border = TRUE
)
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure6 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure6a #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set.seed(2023)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A','#E42A2A', '#A31D1D'))(100)

dfjatac <- read.delim("Jin668.txt", header = T, row.names = 1, sep='\t', stringsAsFactors = F, check.names = F)
head(dfjatac)
dim(dfjatac)
dfjatac <- dfjatac %>% dplyr::filter_all(any_vars(. != 0))  # remove row that all is zero
dim(dfjatac)
column_order <- c("J0", "J1", "J2", 'J3', "J4")
dfjatac <- dfjatac[column_order]

scaled_z_score_row_dfjatac <- t(apply(dfjatac, 1, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
z_score_col_fun <- colorRamp2(seq(-2, 2, by = 4/99), pal_atac)


dfjatach <- Heatmap(
  scaled_z_score_row_dfjatac, 
  col = z_score_col_fun,
  row_km = length(colnames(dfjatac)),      # Split by k-means clustering
  cluster_rows = FALSE, cluster_columns = FALSE, 
  show_row_dend = FALSE, show_column_dend = FALSE,
  column_order = column_order,
  border = TRUE,
  show_row_names = F,
  column_names_rot = 360, 
  heatmap_legend_param = list(
    title = "ATAC-seq RPKM-Jin668\nZ-score",
    legend_width  = unit(3, "cm"),
    legend_direction = "horizontal",
    labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 12)
  )
)


pdf("Figure6a.pdf", width = 6, height = 8)
dfjatach <- draw(object = dfjatach, newpage = TRUE, merge_legends = TRUE, column_title_gp = gpar(fontsize = 12),heatmap_legend_side = "bottom")
dev.off()


rcl.list <- row_order(dfjatach)  #Extract clusters (output is a list)
dfjatac$PeakID <- row.names(dfjatac)
dfjataccluster <- lapply(names(rcl.list), function(i){
  out <- data.frame(PeakID = rownames(scaled_z_score_row_dfjatac[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i), stringsAsFactors = FALSE)
  return(out)
}) %>% do.call(rbind, .) %>% left_join(dfjatac, "PeakID")

head(dfjataccluster)
dim(dfjatac)
dim(dfjataccluster)

write.table(dfjataccluster, file= "Figure6a_clusters.csv", sep=",", quote=F, row.names=FALSE)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure6b #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(monaLisa)
library(JASPAR2022)
library(TFBSTools)
library(Biostrings)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)

## loading genome and homer
bedpeaks <- "Jin668_ATAC_peaks.bed"
genome <- "Jin668_genome.fasta"
myhomerfile <- "/public/home/bin/findMotifsGenome.pl"
myoutdir <- "/public/home/Jin668_ATAC_peaks_cluster_homer_Motifs"
motiffile <- "/public/home/ATAC_motif/Comprehensive_JASPAR2022_plants_known_homer.motifs"
prefix <- "Jin668"


gr <- rtracklayer::import(con = bedpeaks, format = "bed")
bins <- factor(gr$name)

se <- calcBinnedMotifEnrHomer(gr = gr, b = bins, genomedir = genome,
                              outdir = myoutdir, 
                              motifFile = motiffile, homerfile = myhomerfile,
                              regionsize = 50, Ncpu = 2, verbose = TRUE, verbose.Homer = TRUE)

my_sig_value <- as.numeric(1.3)
for (i in 1:length(unique(gr$name))){
  print(paste0("Significantly enriched in bin/cluster ", unique(gr$name)[i]))
  selbin <- assay(se, "negLog10Padj")[, i] > my_sig_value
  print(sum(selbin, na.rm = TRUE))
  if(sum(selbin, na.rm = TRUE) < 2){
    message(paste0("Total enriched motif less than 2 in ", unique(gr$name)[i], " when set significantly to ", my_sig_value, " !!!"))
  }else{
    message(paste0("Total enriched ", sum(selbin, na.rm = TRUE), " motif in ", unique(gr$name)[i], " when set significantly to ", my_sig_value, " !!!"))
    selbinSel <- se[selbin, ]
    SimMatSel <- motifSimilarity(rowData(selbinSel)$motif.pfm, BPPARAM = BiocParallel::MulticoreParam(2))
    range(SimMatSel)
    pdf(paste0(myoutdir, "/Figure6b_", prefix, "_ATAC_peaks_", unique(gr$name)[i], "_Sig", my_sig_value, "Enriched_Motif_heatmap.pdf"), width = 10, heigh = 10)
    hcl <- hclust(as.dist(1 - SimMatSel), method = "average")
    print(
      plotMotifHeatmaps(x = selbinSel, which.plots = c("log2enr", "negLog10Padj"), 
                        width = 2, cluster = hcl, maxEnr = 2, maxSig = 10,
                        show_dendrogram = TRUE, show_seqlogo = TRUE, show_motif_GC = TRUE,
                        width.seqlogo = 1.2)
    )
    dev.off()
  }
}
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure7b ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(ggpubr)
df <- read.delim("data.csv", header = T, row.names = 1, sep=',', stringsAsFactors = F, check.names = F)

## jin668
dfj <- df %>% filter(species == "Jin668")
jin668 <- ggbarplot(dfj, x = "group", y = "CPR", add = "mean_se",
                    fill = "group", palette = "jco", 
                    position = position_dodge(0.8))+
  stat_compare_means(aes(group = group), label = "p.signif", ref.group = "P7N", method = "t.test") +
  labs(x = "Jin668")+
  ylab("Callus proliferation rate (CPR %)") +
  scale_fill_npg() +
  ylim(c(0,1.25))+
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          legend.position="top",
                          strip.background=element_rect(fill="white"))

## TM1
dft <- df %>% filter(species == "TM1")
tm1 <- ggbarplot(dft, x = "group", y = "CPR", add = "mean_se",
                 fill = "group", palette = "jco", 
                 position = position_dodge(0.8))+
  stat_compare_means(aes(group = group), label = "p.signif", ref.group = "P7N", method = "t.test") +
  ylim(c(0,1.25))+
  labs(x = "TM1")+
  ylab("Callus proliferation rate (CPR %)") +
  scale_fill_npg() +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          legend.position="top",
                          strip.background=element_rect(fill="white"))

## YZ1
dfy <- df %>% filter(species == "YZ1")
yz <- ggbarplot(dfy, x = "group", y = "CPR", add = "mean_se",
                fill = "group", palette = "jco", 
                position = position_dodge(0.8))+
  stat_compare_means(aes(group = group), label = "p.signif", ref.group = "P7N", method = "t.test") +
  ylim(c(0,1.25))+
  labs(x = "YZ1")+
  ylab("Callus proliferation rate (CPR %)") +
  scale_fill_npg() +
  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          legend.position="top",
                          strip.background=element_rect(fill="white"))

pdf("Figure7b.pdf", width = 10, height = 8)
ggpubr::ggarrange(jin668, yz, tm1, ncol = 3, nrow = 1, 
                  common.legend = TRUE, legend = "top")
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Figure8 ####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure8a #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
title_color1 <- c("#e5aecd")
title_color_list1 <- rep(title_color1, times=c(2))
strip1 <- strip_themed(background_x = elem_list_rect(fill = title_color_list1), 
                       background_y = element_rect(fill = "white"))

title_color2 <- c("#8fd9c8")
title_color_list2 <- rep(title_color2, times=c(2))
strip2 <- strip_themed(background_x = elem_list_rect(fill = title_color_list2), 
                       background_y = element_rect(fill = "white"))

title_color3 <- c("#99c2ec")
title_color_list3 <- rep(title_color3, times=c(2))
strip3 <- strip_themed(background_x = elem_list_rect(fill = title_color_list3), 
                       background_y = element_rect(fill = "white"))

title_color4 <- c("#d6cb9c")
title_color_list4 <- rep(title_color4, times=c(12))
strip4 <- strip_themed(background_x = elem_list_rect(fill = title_color_list4), 
                       background_y = element_rect(fill = "white"))


dfatcg <- read.table("genome_and_CDS_ATCG_percentage.txt", header = F)
dfatcg$V2 <- factor(dfatcg$V2, levels = c("genome","CDS"))
dfatcg$V3 <- factor(dfatcg$V3, levels = c("A","T","C","G"))

dfatcg_plot <- ggplot(dfatcg, aes(x=V3, y=V1, group=genome)) +
  geom_point(aes(color=genome, size=V4), shape=19) +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(V2), scales = "free", space = "free", strip = strip2) +
  labs(x = "", y = "") +
  scale_color_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank()
  )

## PAM
dfpam <- read.table("genome_and_CDS_NGG_TTTN_PAM.stat", header = F)
dfpam$V3 <- factor(dfpam$V3, levels = c("NGG","TTTN"))
dfpam$V2 <- factor(dfpam$V2, levels = c("genome","CDS"))

dfpam_plot <- ggplot(dfpam, aes(x=V3, y=V1, group=genome)) +
  geom_point(aes(color=genome, size=V4), shape=15) +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(V2), scales = "free", space = "free", strip = strip3) +
  labs(x = "", y = "") +
  scale_color_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank()
  )

## sgRNA
dfsgRNAcount <- read.table("sgRNA_and_Offtarget_Count.txt", header = F)
colnames(dfsgRNAcount) <- c("strand", "type", "GrafEtAlStatus", "PAM", "species", "PAMsum")
dfsgRNAcount <- dfsgRNAcount %>%
  ddply(c("type","GrafEtAlStatus","PAM","species","genome"), summarise, NewPAMsum=sum(PAMsum))

dfsgRNAcount$type <- factor(dfsgRNAcount$type, levels = c("sgRNA","Offtargets"))
dfsgRNAcount$GrafEtAlStatus <- factor(dfsgRNAcount$GrafEtAlStatus, levels = c("GrafOK","ggc","tt"))
dfsgRNAcount$PAM <- factor(dfsgRNAcount$PAM, levels = c("Cas9","Cas12"))

dfsgRNAcount_plot <- ggplot(dfsgRNAcount, aes(x=GrafEtAlStatus, y=species, group=genome)) +
  geom_point(aes(color=genome, size=NewPAMsum), shape=19) +
  ggh4x::facet_grid2(genome ~ PAM + GrafEtAlStatus, scales = "free", space = "free", strip = strip4) +
  labs(x = "", y = "") +
  scale_color_npg() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

mylabeller <- labeller(
  PAM = ~ paste("PAM: ", .),
  GrafEtAlStatus = ~ paste("GrafEtAlStatus: ", .),
  .multi_line = FALSE
)
dfsgRNAcount_plot2 <- ggplot(data=dfsgRNAcount, aes(x=species, y=NewPAMsum, fill=genome)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  ggh4x::facet_grid2(genome ~ PAM + GrafEtAlStatus, scales = "free", space = "free_y", labeller = mylabeller, strip = strip4) +
  labs(x = "", y = "") +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        #strip.text.y = element_blank()
  )

dfsgRNAcount_plot3 <- ggplot(data=dfsgRNAcount, aes(x=species, y=NewPAMsum)) +
  ggplot2::geom_segment(aes(x = species, xend = species, y = 0, yend = NewPAMsum), color="skyblue") +
  geom_point( color="blue", size=1, alpha=0.6) +
  coord_flip() +
  ggh4x::facet_grid2(genome ~ PAM + GrafEtAlStatus, scales = "free", space = "free_y", labeller = mylabeller, strip = strip4) +
  labs(x = "", y = "") +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
  )


## genome size
dfsizecount <- read.table("genome_and_CDS_size_count.txt", header = F)
colnames(dfsizecount) <- c("species", "type", "size_count")
dfgenomesize_plot <- ggplot(data=dfsizecount %>% filter(type=="genome_size"), aes(x=species, y=size_count, fill=genome)) +
  geom_bar(stat="identity") + coord_flip() +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(type), scales = "free", space = "free", strip = strip1) +
  labs(x = "", y = "") +
  scale_fill_npg() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text.y = element_blank()
  )

dfgenecount_plot <- ggplot(data=dfsizecount %>% filter(type=="gene_count"), aes(x=species, y=size_count, fill=genome)) +
  geom_bar(stat="identity") + coord_flip() +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(type), scales = "free", space = "free", strip = strip1) +
  labs(x = "", y = "") +
  scale_fill_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank()
  )



pdf("Figure8a.pdf", width = 20, height = 8)
ggpubr::ggarrange(dfgenomesize_plot, dfgenecount_plot, 
                  dfatcg_plot, dfpam_plot,
                  dfsgRNAcount_plot2,
                  ncol = 5, nrow = 1, 
                  widths = c(0.65, 0.25, 0.5, 0.5, 2), 
                  heights = c(1, 1, 1, 1, 1),
                  common.legend = TRUE, legend = "top")
ggpubr::ggarrange(dfgenomesize_plot, dfgenecount_plot, 
                  dfatcg_plot, dfpam_plot,
                  dfsgRNAcount_plot3,
                  ncol = 5, nrow = 1, 
                  widths = c(0.65, 0.25, 0.5, 0.5, 2), 
                  heights = c(1, 1, 1, 1, 1),
                  common.legend = TRUE, legend = "top")
dev.off()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Figure8d #####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(Rtsne)
library(ggfortify)
title_color_list_y <- c("#b5d3b4", "#cab6d7", "#d6cb9c")
## ~~~~~~~~~~~~~~~~~~~~~~ Cas9
dftmcas <- read.csv("Cas9_targets.txt", sep = "\t", header = T)

## tsne
dftmcas_tsne_results <- Rtsne(dftmcas_data, perplexity=30, check_duplicates = FALSE)

## PCA
cas9_pca_res <- prcomp(dftmcas_data, scale. = TRUE)
cas9_pca_plot <- autoplot(cas9_pca_res, data = dftmcas, colour = 'genome', 
                          loadings = TRUE, loadings.colour = 'blue',
                          loadings.label = TRUE, loadings.label.size = 2, 
                          frame = TRUE, frame.type = 'norm') + 
  scale_fill_manual(values = title_color_list_y) +
  scale_color_manual(values = title_color_list_y) +
  theme_few()


## ~~~~~~~~~~~~~~~~~~~~~~ Cas12
dftmcpf <- read.csv("Cas12_targets.txt", sep = "\t", header = T)

## tsne
dftmcpf_tsne_results <- Rtsne(dftmcpf_data, perplexity=30, check_duplicates = FALSE)

## PCA
cpf_pca_res <- prcomp(dftmcpf_data, scale. = TRUE)
cpf_pca_plot <- autoplot(cpf_pca_res, data = dftmcpf, colour = 'genome', 
                         loadings = TRUE, loadings.colour = 'blue',
                         loadings.label = TRUE, loadings.label.size = 2, 
                         frame = TRUE, frame.type = 'norm') + 
  scale_fill_npg() +
  theme_few()


pdf("Figure8d_tsne.pdf", width = 4, height = 4)
plot(dftmcas_tsne_results$Y, col = "black", bg= factor(dftmcas$genome), pch = 21, cex = 1)
plot(dftmcpf_tsne_results$Y, col = "blue", bg= factor(dftmcpf$genome), pch = 21, cex = 1)
dev.off()

pdf("Figure8d_pca.pdf", width = 8, height = 4)
ggpubr::ggarrange(cas9_pca_plot, cpf_pca_plot,
                  ncol = 2, nrow = 1, 
                  common.legend = TRUE, legend = "top")
dev.off()

library(stringr)
title_color_list_X_cas <- c("#FEF8C6")
title_color_list_X_cpf <- c("#CCE5F8")
title_color_list_y <- c("#b5d3b4", "#cab6d7", "#d6cb9c")
strip_cas <- strip_themed(background_x = elem_list_rect(fill = title_color_list_X_cas),
                          background_y = elem_list_rect(fill = title_color_list_y))

strip_cpf <- strip_themed(background_x = elem_list_rect(fill = title_color_list_X_cpf),
                          background_y = elem_list_rect(fill = title_color_list_y))
# ~~~~~~~~~~~~~ Cas9
dftmcas <- read.csv("Cas9_targets.txt", sep = "\t", header = T)
dftmcas <- dftmcas %>% 
  mutate(strand = case_when(grepl("forw", guideId) ~ "forw",
                            grepl("rev", guideId) ~"rev")) %>% 
  dplyr::mutate(across('guideId', str_replace, 'forw', '')) %>% 
  dplyr::mutate(across('guideId', str_replace, 'rev', '')) %>% 
  dplyr::mutate(guideId = as.integer(guideId)) %>% 
  dplyr::mutate(sgRNAlength = as.integer(20)) %>% 
  dplyr::mutate(end = guideId + sgRNAlength) %>%
  dplyr::mutate(Cas = "Cas9") %>%
  dplyr::arrange(guideId)

dftmcas$myindex <- row.names(dftmcas)

cas_plot <- ggplot(dftmcas) +
  geom_point(aes(x=guideId, y=myindex, group=strand, color=strand), size = 0.1) +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(Cas), scales = "free", space = "free_y", strip = strip_cas) +
  labs(x = "", y = "") +
  scale_color_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
  )

# ~~~~~~~~~~~~~ Cas12
dftmcpf <- read.csv("Cas12_targets.txt", sep = "\t", header = T)
dftmcpf <- dftmcpf %>% 
  mutate(strand = case_when(grepl("forw", guideId) ~ "forw",
                            grepl("rev", guideId) ~"rev")) %>% 
  dplyr::mutate(across('guideId', str_replace, 'forw', '')) %>% 
  dplyr::mutate(across('guideId', str_replace, 'rev', '')) %>% 
  dplyr::mutate(guideId = as.integer(guideId)) %>% 
  dplyr::mutate(sgRNAlength = as.integer(23)) %>% 
  dplyr::mutate(end = guideId + sgRNAlength) %>%
  dplyr::mutate(Cas = "Cas12") %>%
  dplyr::arrange(guideId)

dftmcpf$myindex <- row.names(dftmcpf)

cpf_plot <- ggplot(dftmcpf) +
  geom_point(aes(x=guideId, y=myindex, group=strand, color=strand), size = 0.1) +
  ggh4x::facet_grid2(rows = vars(genome), cols = vars(Cas), scales = "free", space = "free_y", strip = strip_cpf) +
  labs(x = "", y = "") +
  scale_color_npg() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank()
  )

pdf("Figure8d_targets_located_site.pdf", width = 4, height = 6)
ggpubr::ggarrange(cpf_plot,cas_plot,
                  ncol = 2, nrow = 1, 
                  common.legend = TRUE, legend = "top")
dev.off()