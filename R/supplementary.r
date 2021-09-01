make.hit.ratio.figure <- function(outfile) {
  HASH_PATH <- "./results/hash_comparison/"

  species <- c("hg38", "pantro6", "mm10", "danre11", "galgal6", "tair10")
  species.names <- c("H. sapiens", "P. troglodytes",
                     "M. musculus", "D. rerio",
                     "G. gallus", "A. thaliana")

  tables <- list()

  pdf(outfile, width = 9, height = 5)
  par(mfrow = c(2, 3), family = "Times")
  for (i in 1:length(species)) {
    tbl <- read.table(paste0(HASH_PATH, species[i], "_hash.tsv"),
                              header = 1, row.names=NULL)

    ks <- tbl[,1]
    real <- tbl[,7]
    theoretical <- tbl[,8]
    upper.bound <- tbl[,12]

    # real collision by counting kmers and summing squares
    plot(ks, real, type = 'o', xlab = "number of bits", ylab = "Expected hit rate",
         col = "black", lty = 1, cex = .5, pch = 19, ylim = c(0,1))

    # theoretical collision by getting weak and strong probabilities and
    # using iid theory
    points(ks, theoretical, type = 'o', xlab = "number of bits",
           ylab = "ratio of expected hit rates", col = "#666666", lty = 2, cex = .5, pch = 19)

    # upper bound devised in the equation
    points(ks, upper.bound, type = 'o',
           lty = 1, cex = .5, pch = 19, col = "#FF0000")
    title(main = species.names[i])
    legend("topright", lty = c(1,2,1), col = c("black", "#666666","#FF0000"),
           legend = c("empirical", "theoretical", "0.93^b"))
  }
  dev.off()
}
