This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
After adapter removal with `Trimmomatic-0.39 <http://www.usadellab.org/cms/?page=trimmomatic>`_, reads were mapped and gene counts were generated with `HISAT2 <https://daehwankimlab.github.io/hisat2/>`_, and `StringTie <https://ccb.jhu.edu/software/stringtie/>`_.
Gene counts of replicated were summed up.
Integrated normalization and differential expression analysis was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ following standard procedure as outlined in the manual.
