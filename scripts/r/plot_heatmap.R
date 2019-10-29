

if(!require(gplots)){
  install.packages("gplots", repos='http://cran.us.r-project.org')
  library(gplots)
}

if(!require(readr)){
  install.packages("readr", repos='http://cran.us.r-project.org')
  library(readr)
}

if(!require(optparse)){
  install.packages("optparse", repos='http://cran.us.r-project.org')
  library(optparse)
}


parseArgs <- function() {
  
  option_list <- list(
    make_option(c("-i", "--input"),
                action = "store",
                dest = "input",
                metavar = "FILE",
                type = "character",
                help = "Matrix file"),
    
    make_option(c("-m", "--max"),
                action = "store",
                dest = "max_val",
                type = "integer",
                default = 255,
                help = "Max value for heatmap")
    
  )
  
  args <- parse_args(OptionParser(option_list = option_list))
  
  if (!file.exists(args$input)) {
    stop(paste0("Cannot open file ", args$input, "!"))
  }
  
  args
}


plot_heatmap <- function(input, output_png, output_pdf, max_val){
  
  value_table <- read_delim(input, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  ramp <- colorRamp(c("white","pink","red","black"))
  
  
  png(output_png, width = 8, height = 8, units = 'in', res = 300)
  a <- dev.cur()
  pdf(output_pdf, width = 8, height = 8)
  dev.control("enable")
  a_heatmap <- heatmap.2(cexRow=0.5,cexCol=0.5,as.matrix(value_table), Rowv=FALSE, Colv=FALSE, 
                         col = rgb( ramp(seq(0, 1, length = 200)), max = max_val),scale="none",density.info="none",trace="none",dendrogram="none")
  dev.copy(which = a)
  dev.off()
  dev.off()

  
}


main <- function(){
  args <- parseArgs()
  f_name <- unlist(strsplit(args$input, "\\."))
  f_png <- paste(c(f_name[-(length(f_name))], "png"), collapse = ".")
  f_pdf <- paste(c(f_name[-(length(f_name))], "pdf"), collapse = ".")
  
  plot_heatmap(args$input, f_png, f_pdf, args$max_val)
  
}

main()