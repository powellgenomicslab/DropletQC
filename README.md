
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dropletQC

<!-- badges: start -->
<!-- badges: end -->

This is a simple R package to calculate a QC metric, the nuclear
fraction score, for single cell RNA-seq (scRNA-seq) datasets generated
using the 10x Genomics Chromium Single Cell Gene Expression platform.
This statistic is simply:

    nuclear fraction = intronic reads / (intronic  reads  +  exonic  reads)

In words: for each cell barcode provided, the proportion of reads that
originated from intronic regions is calculated. These RNA fragments
likely originate from unspliced (nuclear) pre-mRNA, hence the name
“nuclear fraction”. This metric can be useful to help identify two
populations that you may want to exclude from your dataset prior to
downstream analysis:

1.  “Empty” droplets containing ambient RNA, characterised by a low
    nuclear fraction score

2.  Droplets containing damaged cells, characterised by a high nuclear
    fraction score

The biological principle behind this is intuitive. Sheared cell
membranes of damaged cells in the input cell suspension release
cytoplasmic RNA into solution while the nuclear envelope will often
remain intact. As a result, RNA released from stressed or damaged cells
will consist of mostly mature cytoplasmic mRNA and will be relatively
depleted of unspliced nuclear precursor mRNA.

## Installation

You can install dropletQC from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("WalterMuskovic/dropletQC")
```

## Examples

The simplest way to calculate the nuclear fraction score is to simply
point to the ‘outs’ directory produced by the Cell Ranger software:

``` r
library(dropletQC)
nf1 <- nuclear_fraction(
    outs = system.file("extdata", "outs", package = "dropletQC"),
     tiles = 10, cores = 1, verbose = FALSE)
head(nf1)
#>                    nuclear_fraction
#> AAAAGTCACTTACTTG-1       0.05025126
#> AAAAGTGGATCTCTAA-1       0.03804348
#> AAACACGTTCTCATCG-1       0.02985075
#> AAACAGGCAGCGACTG-1       0.06451613
#> AAAGCAGTTACGAAGA-1       0.04624277
#> AAAGCGGATGCATGGT-1       0.03821656
```

This assumes the following three files are present in the specified
directory:

``` r
list.files(system.file("extdata", "outs", package = "dropletQC"), recursive = TRUE)
#> [1] "filtered_feature_bc_matrix/barcodes.tsv.gz"
#> [2] "possorted_genome_bam.bam"                  
#> [3] "possorted_genome_bam.bam.bai"
```

Alternatively, if you don’t have this directory structure or your files
have been renamed e.g. they were given to you by a collaborator, you can
specify the paths to the required files directly:

``` r
nf2 <- nuclear_fraction(
   bam = system.file("extdata", "outs","possorted_genome_bam.bam", package =
   "dropletQC"),
   barcodes = c("AAAAGTCACTTACTTG-1",
                "AAAAGTGGATCTCTAA-1",
                "AAACACGTTCTCATCG-1"),
   tiles = 10, cores = 1,
   verbose = FALSE)
nf2
#>                    nuclear_fraction
#> AAAAGTCACTTACTTG-1       0.05025126
#> AAAAGTGGATCTCTAA-1       0.03804348
#> AAACACGTTCTCATCG-1       0.02985075
```

Note that here we have provided a vector of requested barcode IDs to the
`barcodes` argument rather than the path to a file on disk
`barcodes.tsv.gz`. Either is fine, just make sure the format of your
barcodes matches the BAM file.

## More information

For more details see our paper published in **Journal Name**:

[paper title here](https://www.google.com)

and the associated GitHub repo:

[public GitHub repo](https://www.google.com)

For more information about how the package works and some tips on how to
interpret the nuclear fraction metric using real-world examples, see
`vignette("doubletQC-vignette")`.
