library(biomaRt)  # biomaRt_2.30.0
library(optparse)

### handle command line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="rsID file name (1 rsID per line)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.bed",
              help="output file name [default= %default]", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("at least one argument must be supplied (rsID file)", call.=FALSE)
}
###

snp_ids <- list(read.table(file = opt$file, header = FALSE, stringsAsFactors = FALSE)$V1)
snp_attributes <- c("chr_name", "chrom_start", "chrom_end", "refsnp_id")  # ask for chr number, start coordinate, and rsID

# connect to the hg19 ensembl database
snp_mart <- useMart("ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")

# retrieve requested attributes from database
snp_locations <- getBM(attributes=snp_attributes, filters="snp_filter", values=snp_ids, mart=snp_mart)

# format results returned from function
snp_locations$chr_name <- paste("chr", snp_locations$chr_name, sep = "")
snp_locations$chrom_start <- snp_locations$chrom_start - 1  # because https://www.biostars.org/p/84686/

# save results to BED file
write.table(snp_locations, file = opt$out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
