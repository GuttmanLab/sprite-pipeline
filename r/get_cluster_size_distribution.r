#!/usr/bin/env Rscript

library(ggplot2)

SINGLETON        <- 1
FROM_2_TO_10     <- 2
FROM_11_TO_100   <- 3
FROM_101_TO_1000 <- 4
OVER_1000        <- 5

getCount <- function(x, n) x[[n]]

incrementCount <- function(x, cat, num) {
    x[[cat]] <-  x[[cat]] + num
    x
}

normalizeCount  <- function(x) {
    total <- sum(unlist(x))
    lapply(x, function(l) l / total)
}

processFile <- function(f) {
    sizes <- list(0, 0, 0, 0, 0)

    tryCatch({
        con <- file(f, open = 'r')
        while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
            numReads <- length(unlist(strsplit(oneLine, "\t"))) - 1  # First column is barcode
            if      (numReads == 1)    sizes <- incrementCount(sizes, SINGLETON, numReads)
            else if (numReads <= 10)   sizes <- incrementCount(sizes, FROM_2_TO_10, numReads)
            else if (numReads <= 100)  sizes <- incrementCount(sizes, FROM_11_TO_100, numReads)
            else if (numReads <= 1000) sizes <- incrementCount(sizes, FROM_101_TO_1000, numReads)
            else                       sizes <- incrementCount(sizes, OVER_1000, numReads)
        }
    }, finally = {
        close(con)
    })

    sizes <- normalizeCount(sizes)
    df <- data.frame(unlist(sizes))
    colnames(df) <- c("count")
    fs <- c("SINGLETON", "FROM_2_TO_10", "FROM_11_TO_100", "FROM_101_TO_1000", "OVER_1000")
    df$category <- factor(fs, levels = c(fs[1], fs[2], fs[3], fs[4], fs[5]))
    df$filename <- basename(f)
    df
}

parseArgs <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) stop("get_cluster_sizes <directory> <pattern>")
    args
}

plotCounts <- function(df) {
    p <- ggplot(df, aes(x = filename, y = count, fill = category))
    p <- p + geom_bar(stat = "identity")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    png("cluster_sizes.png")
    print(p)
    graphics.off()

    pdf("cluster_sizes.pdf")
    print(p)
    graphics.off()
}

main <- function() {
    args <- parseArgs()
    files <- list.files(path = args[1], pattern = args[2], full.names = TRUE, recursive = FALSE)
    dfs <- lapply(files, processFile)
    df <- do.call("rbind", dfs)
    plotCounts(df)
}

main()
