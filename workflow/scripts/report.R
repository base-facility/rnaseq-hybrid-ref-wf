suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description="Generate RNAseq report")
parser$add_argument("-w", "--id",  dest="sample_id", help="Sample ID")
parser$add_argument("-q", "--quant",  dest="quant", help="Salmon quant output (quant.sf)")
parser$add_argument("-d", "--output_dir",  dest="output_dir", help="Report output directory path")
parser$add_argument("-f", "--output_file",  dest="output_file", help="Report output file path")
parser$add_argument("-k", "--knit_root",  dest="knitroot", help="Knit root directory path")
parser$add_argument("-r", "--synth",  dest="synth_name", help="Synthetic construct name", nargs="+")
args <- parser$parse_args()

# debug
message("Sample ID:", args$sample_id)
message("Quant file:", args$quant)
message("Output directory:", args$output_dir)
message("Output file:", args$output_file)
message("Knit root directory:", args$knitroot)
message("Synthetic construct name:", args$synth_name)

rmarkdown::render("workflow/scripts/report.Rmd",
                  params=list(sample_id=args$sample_id,
                              quant=args$quant,
                              synth_name=args$synth_name),
                  output_format="html_document",
                  output_dir=args$output_dir,
                  output_file=args$output_file,
                  knit_root_dir=args$knitroot)