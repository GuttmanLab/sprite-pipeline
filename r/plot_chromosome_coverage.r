#!/usr/bin/env Rscript

# This program plots SPRITE clusters as a barplot, with each bar corresponding
# to the read-coverage over a chromosome. If passed only a single clusters
# file, this program will plot absolute coverage   

library(ggplot2)
library(optparse)

refs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
          "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
          "chr16", "chr17", "chr18", "chr19", "chrX")

parseArgs <- function() {
  
    option_list <- list(
        make_option(c("-i", "--input"),
                    action = "store",
                    dest = "input",
                    metavar = "FILE",
                    type = "character",
                    help = "blah"),

        make_option(c("-n", "--normalize"),
                    action = "store",
                    dest = "normalize",
                    metavar = "FILE",
                    type = "character",
                    default = NULL,
                    help = "blah"),
  
        make_option(c("-o", "--output"),
                    action = "store",
                    dest = "output",
                    metavar = "FILE",
                    type = "character",
                    help = "blah")
      )

      args <- parse_args(OptionParser(option_list = option_list))

      if (!file.exists(args$input)) {
          stop(paste0("Cannot open file ", args$input, "!"))
      }

      if (!is.null(args$normalize) && !file.exists(args$normalize)) {
          stop(paste0("Cannot open file ", args$normalize, "!"))
      }
  
      args
}

processFile <- function(input) {
    hash <- new.env()
    sapply(refs, function(x) hash[[x]] <- 0)
  
    tryCatch({
        con <- file(input, open = "r")
        while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
            fields <- unlist(strsplit(oneLine, "\t"))[-1]
            for (i in 1:length(fields)) {
                chrom <- getChrom(fields[i])
                count <- hash[[chrom]]
                if (!is.null(count)) {
                    hash[[chrom]] <- count + 1
                }
            }
        }
    }, finally = {
    close(con)
    })
    
    hash
}

plotHash <- function(hash, output) {
    chroms <- names(hash)
    df <- data.frame(matrix(ncol = 2, nrow = length(chroms)))
    colnames(df) <- c("chromosome", "count")
    for (i in 1:length(chroms)) {
        df[i, ] <- c(chroms[i], hash[[chroms[i]]])
    }
    
    df$chromosome <- factor(df$chromosome, levels = refs)
    
    p <- ggplot(df, aes(x = chromosome, y = as.numeric(count)))
    p <- p + geom_bar(stat = "identity")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p <- p + ylab("coverage (absolute)")

    ggsave(paste0(output, ".pdf"), plot = p)
    ggsave(paste0(output, ".eps"), plot = p)
}

plotHashes <- function(numHash, denomHash, output) {
    chroms <- names(numHash)
    df <- data.frame(matrix(ncol = 3, nrow = length(chroms)))
    colnames(df) <- c("chromosome", "num", "denom")
    for (i in 1:length(chroms)) {
        df[i, ] <- c(chroms[i], numHash[[chroms[i]]], denomHash[[chroms[i]]])
    }

    df$num <- as.numeric(df$num)
    df$denom <- as.numeric(df$denom)
    df$num.frac <- df$num / sum(df$num)
    df$denom.frac <- df$denom / sum(df$denom)
    df$oe <- df$num.frac / df$denom.frac

    df$chromosome <- factor(df$chromosome, levels = refs)

    p <- ggplot(df, aes(x = chromosome, y = oe))
    p <- p + geom_bar(stat = "identity")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p <- p + ylab("O/E")

    ggsave(paste0(output, ".pdf"), plot = p)
    ggsave(paste0(output, ".eps"), plot = p)
}

getChrom <- function(s) unlist(strsplit(s, ":"))[1]

main <- function() {
    args <- parseArgs()
    hash <- processFile(args$input)
    if (is.null(args$normalize)) {
        plotHash(hash, args$output)
    } else {
        hash2 <- processFile(args$normalize)
        plotHashes(hash, hash2, args$output)
    }
}

main()
