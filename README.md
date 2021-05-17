
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DropletQC

<!-- badges: start -->
<!-- badges: end -->

This is a simple R package to calculate, for every requested cell
barcode in a provided scRNA-seq BAM file, the nuclear fraction score:

    nuclear fraction = intronic reads / (intronic  reads  +  exonic  reads)

The score captures the proportion of reads from intronic regions. These
RNA fragments originate from unspliced (nuclear) pre-mRNA, hence the
name “nuclear fraction”. This score can be used to help identify:

1.  “Empty” droplets containing ambient RNA: low nuclear fraction score
    and low UMI count

2.  Droplets containing damaged cells: high nuclear fraction score and
    low UMI count

## Installation

You can install DropletQC with:

``` r
# install.packages("devtools")
devtools::install_github("WalterMuskovic/DropletQC", build_vignettes = TRUE)
```

## Calculating the nuclear fraction

There are two functions which can be used to calculate the nuclear
fraction; `nuclear_fraction_tags` and `nuclear_fraction_annotation`.

If your BAM file contains region tags which identify aligned reads as
intronic or exonic, such as those produced by 10x Genomics’ Cell Ranger
software, then the simplest and fastest way to calculate the nuclear
fraction is to point `nuclear_fraction_tags` to the directory:

``` r
library(DropletQC)
nf1 <- nuclear_fraction_tags(
    outs = system.file("extdata", "outs", package = "DropletQC"),
     tiles = 1, cores = 1, verbose = FALSE)
head(nf1)
#>                    nuclear_fraction
#> AAAAGTCACTTACTTG-1        0.9032698
#> AAAAGTGGATCTCTAA-1        0.4032761
#> AAAGCAGTTACGAAGA-1        0.3957704
#> AACGACTTCAATATGT-1        0.4004525
#> AACGGCGTCATCTGGA-1        0.8845109
#> AAGCAGGGGTCGCGAA-1        0.3929376
```

Alternatively, you can point `nuclear_fraction_annotation` to a gene
annotation, BAM and barcode files:

``` r
nf2 <- nuclear_fraction_annotation(
 annotation_path = system.file("extdata/outs/chr1.gff3",package = "DropletQC"),
 bam = system.file("extdata/outs/possorted_genome_bam.bam",package = "DropletQC"),
 barcodes = system.file("extdata/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",package = "DropletQC"),
 tiles = 1, cores = 1, verbose = FALSE)
head(nf2)
#>                    nuclear_fraction
#> AAAAGTCACTTACTTG-1        0.9032698
#> AAAAGTGGATCTCTAA-1        0.4032761
#> AAAGCAGTTACGAAGA-1        0.3957704
#> AACGACTTCAATATGT-1        0.4004525
#> AACGGCGTCATCTGGA-1        0.8845109
#> AAGCAGGGGTCGCGAA-1        0.3929376
```

This methods is more flexible, as it makes no assumptions about how your
BAM file was produced - but it will take longer. Take care that the
provided barcodes match the barcode structure in the BAM file.

## Identifying empty droplets and damaged cells

Once the nuclear fraction score has been calculated, the
`identify_empty_drops` and `identify_damaged_cells` functions can be
used to assist in identifying each these populations. Empty or damaged
cells are flagged, not removed.

## More information

For a detailed discussion see our paper published in **Journal Name**:

[paper title here](https://www.google.com)

For more information about the functions included in the package,
including tips on how to assess the nuclear fraction score using
real-world examples, see the [package
vignette](https://powellgenomicslab.github.io/DropletQC/articles/DropletQC.html).
