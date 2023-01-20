#!/usr/bin/env Rscript

library(optparse)

parser = OptionParser(description='Run SNP pileup and phasing with 1000G')
parser = add_option(parser, '--label', default = 'subject', type = "character", help = "Individual label. One per run.")
parser = add_option(parser, '--samples', default = 'sample', type = "character", help = "Sample name(s); comma delimited if multiple. All samples must belong to the same individual.")
parser = add_option(parser, '--gmap', default = NULL, type = "character", help = "Path to genetic map provided by Eagle2 (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)")
parser = add_option(parser, '--outdir', default = './pileup_and_phase', type = "character", help = "Output directory")
args <- parse_args(parser)

suppressPackageStartupMessages({
    library(glue)
    library(stringr)
    library(data.table)
    library(dplyr)
    library(vcfR)
    library(Matrix)
    library(numbat)
})

# required args
if (any(is.null(c(args$gmap)))) {
    stop('Missing one or more required arguments: gmap')
}

label = args$label
samples = str_split(args$samples, ',')[[1]]
outdir = args$outdir
n_samples = length(samples)
gmap = args$gmap
genome = ifelse(str_detect(args$gmap, 'hg19'), 'hg19', 'hg38')
message(paste0('Using genome version: ', genome))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
dir.create(glue('{outdir}/phasing'), showWarnings = FALSE)
for (sample in samples) {
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

cat('Reading genetic map\n')
genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

for (sample in samples) {
    
    # read in phased VCF
    cat('Reading topmed VCF files\n')
    vcf_phased = lapply(1:22, function(chr) {
            vcf_file = glue('{outdir}/topmed_phasing/chr{chr}.dose.vcf.gz')
            if (file.exists(vcf_file)) {
                fread(vcf_file, skip = '#CHROM') %>%
                    rename(CHROM = `#CHROM`) %>%   
                    mutate(CHROM = str_remove(CHROM, 'chr'))
            } else {
                stop('Phased VCF not found')
            }
        }) %>%
        Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    cat('Reading pileup vcf file\n')
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
        rename(CHROM = `#CHROM`) %>%
        mutate(CHROM = str_remove(CHROM, 'chr'))

    # count matrices
    cat('Reading count matrices\n')
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

    cat('Reading cell barcodes file\n')
    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

    cat('Generating actual allele count dataframes\n')
    df = numbat:::preprocess_allele(
        sample = label,
        vcf_pu = vcf_pu,
        vcf_phased = vcf_phased,
        AD = AD,
        DP = DP,
        barcodes = cell_barcodes,
        gtf = gtf,
        gmap = genetic_map
    ) %>%
    filter(GT %in% c('1|0', '0|1'))
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
}
