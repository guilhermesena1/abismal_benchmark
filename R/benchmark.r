source("R/fig_label.r")

BENCHMARK_PATH =
"/panfs/qcb-panasas/desenabr/ab/results/benchmark-files"

# whether to redo plots
plot <- T

datasets <- read.table("metadata/tests.txt",
                       header = 1,row.names = NULL,stringsAsFactors = F)
mappers <- c("bismark", "bsmap", "walt")

##########################
# Aux functions
##########################

# turns absolute number of mapped reads into percentage
get.pct <- function(ans) {
  ans$pct.unique <- ans$mapped.unique / ans$total.reads
  ans$pct.total <- ans$mapped.total / ans$total.reads
  ans
}

#########################
# Parse map statistics
#########################

# get useful information from bismark_PE_report.txt
parse.bismark <- function(path) {
  lines <- readLines(path)
  if(length(lines) != 42 & length(lines) != 40)
    stop(paste0("bad bismark file: ", path))

  ans <- list()

  # mapping statistics
  ans$total.reads <- as.integer(gsub("^.*\t","",lines[7]))
  ans$mapped.unique <- as.integer(gsub("^.*\t","",lines[8]))

  ans$unmapped <- as.integer(gsub("^.*\t","",lines[10]))
  ans$mapped.total <- ans$total.reads - ans$unmapped

  ans <- get.pct(ans)
}

# get useful information from walt/abismal outputs
parse.walt <- function(path) {
  if(!file.exists(path))
    stop(paste0(c("bad file: " , path)))
  lines <- readLines(path)
  ans <- list()

  tot <- strsplit(lines[1], " ")[[1]]
  if (length(tot) == 6) { # PE
    tot <- tot[6]
    tot <- gsub("]", "", tot)
    ans$total.reads <- as.integer(tot)
    ans$mapped.unique <- as.integer(strsplit(lines[2], " ")[[1]][5])
    ambig <- as.integer(strsplit(lines[3], " ")[[1]][5])
    ans$mapped.total <- ans$mapped.unique + ambig 
  }
  else if (length(tot) == 5) { # SE
    tot <- tot[5]
    tot <- gsub("]", "", tot)
    ans$total.reads <- as.integer(tot)
    ans$mapped.unique <- as.integer(strsplit(lines[2], " ")[[1]][4])
    ambig <- as.integer(strsplit(lines[3], " ")[[1]][4])
    ans$mapped.total <- ans$mapped.unique + ambig 
  }
  else {
    stop("bad walt file")
  }
  ans <- get.pct(ans)
  ans
}

# get useful information from walt/abismal outputs
parse.abismal <- function(path) {
  if(!file.exists(path))
    stop(paste0(c("bad file: " , path)))
  data <- yaml::yaml.load_file(path)
  ans <- list()

  # paired ended
  if (!is.null(data$pairs))  {
    data <- data$pairs
    ans$total.reads <- data$total_read_pairs
  }

  # single end
  else {
    ans$total.reads <- data$total_reads
  }

  ans$mapped.unique <- data$mapped$num_unique
  ans$unmapped <- data$unmapped
  ans$mapped.total <- ans$mapped.unique + data$mapped$ambiguous

  ans <- get.pct(ans)
  ans
}


# get useful information from the bsmap stderr when -V flag
# is set to 1
parse.bsmap <- function(path) {
  if(!file.exists(path))
    stop(paste0("bad bsmap file: ", path))
  ans <- list()
  data <- readLines(path)
  if (grepl("total read pairs", data[length(data)-3])) {
    total.data <- strsplit(data[length(data) - 3], " +")[[1]]
    mapped.data <- strsplit(data[length(data) - 2], " +")[[1]]
    ans$total.reads <- as.integer(total.data[10])
    ans$mapped.unique <- as.integer(mapped.data[7])
    ans$mapped.total <- as.integer(mapped.data[3]) + as.integer(mapped.data[12])
  } else {
    total.data <- strsplit(data[length(data) - 1], " +")[[1]]
    mapped.data <- strsplit(data[length(data)], " +")[[1]]
    ans$total.reads <- as.integer(total.data[9])
    ans$mapped.unique <- as.integer(mapped.data[7])
    ans$mapped.total <- as.integer(mapped.data[3]) + as.integer(mapped.data[12])
  }
  ans <- get.pct(ans)
  ans
}

#################################################
# parse grount truth comparisons
#################################################

# get useful information from walt/abismal outputs
parse.truth <- function(path) {
  if(!file.exists(path))
    stop(paste0(c("bad file: " , path)))
  data <- yaml::yaml.load_file(path)

  tot <- data$total$input_only + data$total$both + data$single$missed_uniq

  ans <- list()
  ans$accuracy <- data$total$accuracy
  ans$truth.superoptimal <- data$single$superoptimal_uniq / tot
  ans$truth.optimal <- data$single$optimal_uniq / tot
  ans$truth.unknown <- data$single$unknown_uniq / tot
  ans$truth.suboptimal <- data$single$suboptimal_uniq / tot
  ans$truth.false <- data$single$false_uniq / tot
  ans$truth.missed <- data$single$missed_uniq / tot
  ans$truth.dangling <- data$single$dangling_uniq / tot
  ans
}

#################################################
# parse methpipe outputs
#################################################

# get average error rate across all bases in bsrate file
parse.bsrate <- function(bsrate.path) {
  if(!file.exists(bsrate.path))
    stop(paste0("bad bsrate file: ", bsrate.path))

  bs.data <- readLines(bsrate.path)
  return(as.numeric(strsplit(bs.data[1], " ")[[1]][5]))
}

# get useful information from levels
parse.levels <- function(levels.path) {
  if(!file.exists(levels.path))
    stop(paste0("bad levels file: ", levels.path))

  suppressWarnings({
    lev <- yaml::yaml.load_file(levels.path)
  })

  ans <- list()
  ans$meth.cytosine <- lev$cytosines$mean_meth
  ans$meth.cpg <- lev$cpg$mean_meth

  ans$coverage.cytosine <- lev$cytosines$sites_covered_fraction
  ans$coverage.cpg <- lev$cpg$sites_covered_fraction

  ans$cpg.depth <- lev$cpg$mean_depth
  ans$cpg.depth.cov <- lev$cpg$mean_depth_covered
  ans
}

# get methylation levels from symetric CpGs
parse.sym.cpg <- function(sym.cpg.path) {
  if (!file.exists(sym.cpg.path))
    stop(paste0("bad sym cpg file: ", sym.cpg.path))
  a <- read.table(sym.cpg.path)
  return(data.frame(level = a$V5, row.names = paste0(a$V1, "_", a$V2)))
}

########################################
# main function: flatten all statistics
########################################

get.row.for.mapper <- function(srr, mapper, mapstats.suffix, f,
                               species, protocol) {
  ans <- data.frame(
    total.reads = NA,
    map.unique = NA,
    map.total = NA,

    map.accuracy = NA,
    truth.superoptimal = NA,
    truth.optimal = NA,
    truth.unknown = NA,
    truth.suboptimal = NA,
    truth.false = NA,
    truth.missed = NA,
    truth.dangling = NA,

    bs.err = NA,
    coverage.cytosine = NA,
    coverage.cpg = NA,
    cpg.depth.total = NA,
    cpg.depth.cov = NA,
    #tpg.mutation.rate = NA,
    row.names = paste0(srr,"_",species)
  )
  # path to outputs directory
  base.path <- paste0(BENCHMARK_PATH, "/", mapper, "/")
  truth.path <- paste0(BENCHMARK_PATH, "/ground_truth/", mapper, "/")

  # fetch mapping data
  path <- paste0(base.path, srr, "_", species, mapstats.suffix)
  tryCatch({
    mapdata <- f(path)
    ans$map.unique = mapdata$pct.unique
    ans$map.total = mapdata$pct.total
    ans$total.reads = mapdata$total.reads
  }, error = function(e){
    write(paste0("failure opening file: ", path), stderr());
  })

  # bsrate statistics
  path <- paste0(base.path, srr, "_", species, ".bsrate")
  tryCatch({
    bsdata <- parse.bsrate(path)
    ans$bs.err = bsdata
  }, error = function(e) {
    write(paste0("failure opening file: ", path), stderr());
  })

  # accuracy statistics
  path <- paste0(base.path, srr, "_", species,
                 ifelse( #protocol != "wgbs_single" & mapper == "bismark",
                          protocol == "wgbs_paired" & mapper != "walt",
                        ".paired_stats", ".single_stats"))
  tryCatch({
    acc <- parse.truth(path)
    ans$map.accuracy = acc$accuracy
    ans$truth.superoptimal <- acc$truth.superoptimal
    ans$truth.optimal <- acc$truth.optimal
    ans$truth.unknown <- acc$truth.unknown
    ans$truth.suboptimal <- acc$truth.suboptimal
    ans$truth.false <- acc$truth.false
    ans$truth.missed <- acc$truth.missed
    ans$truth.dangling <- acc$truth.dangling
  }, error = function(e) {
    write(paste0("failure opening file: ", path), stderr());
  })

  # levels statistics (coverage for all Cs, covered Cs, genome)
  path <- paste0(base.path, srr, "_", species, ".levels")
  tryCatch({
    levelsdata <- parse.levels(path)
    ans$meth.cytosine <- levelsdata$meth.cytosine
    ans$meth.cpg <- levelsdata$meth.cpg
    ans$coverage.cytosine = levelsdata$coverage.cpg
    ans$cpg.depth.total = levelsdata$cpg.depth
    ans$cpg.depth.cov = levelsdata$cpg.depth.cov
    ans$coverage.cpg = levelsdata$coverage.cpg
    #ans$tpg.mutation.rate = levelsdata$tpg.mutation

  }, error = function(e) {
    write(paste0("failure opening file: ", path), stderr());
  })

  colnames(ans) <- paste0(colnames(ans), ".", mapper)
  ans
}

# get all data for a particular sample
get.row <- function(srr, species, protocol) {
  ans <- list()

  ################# ABISMAL #######################
  a.data <- get.row.for.mapper(srr, "abismal", ".sam.mapstats",
                               parse.abismal, species, protocol)
  tot.reads <- a.data$total.reads
  a.data$total.reads.abismal <- NULL
  ################# WALT #######################
  if ("walt" %in% mappers) {
    w.data <- get.row.for.mapper(srr, "walt", ".mr.mapstats",
                                 parse.walt, species, protocol)

    if (is.na(tot.reads))
      tot.reads <- w.data$total.reads

    if (!is.na(tot.reads) & !is.na(w.data$total.reads))
      if (w.data$total.reads != tot.reads)
        stop(paste0("inconsistent # of reads between walt and abismal: ",
                  srr))
    w.data$total.reads.walt <- NULL
  }
  ################# BISMARK #######################
  if ("bismark" %in% mappers) {
    k.data <- get.row.for.mapper(srr, "bismark",
                                 ".mapstats",
                                 parse.bismark, species, protocol)

    if (is.na(tot.reads))
      tot.reads <- k.data$total.reads

    if (!is.na(tot.reads) & !is.na(k.data$total.reads))
      if (k.data$total.reads != tot.reads)
        stop(paste0("inconsistent # of reads between bismark and abismal: ",
                  srr))
    k.data$total.reads.bismark <- NULL
  }
  ################# BSMAP #######################
  if ("bsmap" %in% mappers) {
    p.data <- get.row.for.mapper(srr, "bsmap", ".out",
                               parse.bsmap, species, protocol)

    if (is.na(tot.reads))
      tot.reads <- p.data$total.reads

    if (!is.na(tot.reads) & !is.na(p.data$total.reads))
      if(p.data$total.reads != tot.reads)
        stop(paste0("inconsistent # of reads between bsmap and abismal: ",
                  srr))
    p.data$total.reads.bsmap <- NULL
  }
  # flatten information from all mappers
  row <- a.data
  if ("bismark" %in% mappers) row <- cbind(row, k.data)
  if ("bsmap" %in% mappers) row <- cbind(row, p.data)
  if ("walt" %in% mappers) row <- cbind(row, w.data)
  row <- row[order(names(row))]
  row
}

# rbinds rows from all samples
make.table <- function(datasets) {
  stats <- NULL
  all.names <- c()
  for(i in 1:nrow(datasets)) {
    srr <- datasets[i,1]
    species <- datasets[i,2]
    protocol <- datasets[i,3]
    id <- datasets[i,4]
    readlen <- datasets[i,5]
    total.reads <- datasets[i,6]
    row <- list()
    row$srr <- srr
    row$species <- species
    row$protocol <- protocol
    row$id <- id
    row$read.length <- readlen
    row$total.reads <- total.reads
    row <- c(row, get.row(srr, species, protocol))
    stats <- rbind(stats, row)
    all.names <- c(all.names, paste0(id," (",species, ")"))
  }
  rownames(stats) <- all.names
  data.frame(stats)
}

# get the elapsed time value from /usr/bin/time -v
# given the file name as input
usr.bin.time <- function(filename) {
  # GS: use readlines here because of mappers
  # that say "1 day, ..." on it
  line <- readLines(filename)[2]
  time <- strsplit(line, "\t")[[1]][1]
  as.numeric(time)
}

# in gigabytes
usr.bin.mem <- function(filename) {
  line <- readLines(filename)[2]
  mem <- strsplit(line, "\t")[[1]][3]
  as.numeric(mem)
}

get.times.for.srr <- function(srr, species) {
  mappers <- c("abismal", mappers)
  times <- c()
  for (i in mappers) {
    file <- paste0(BENCHMARK_PATH, "/", i, "/snakemake_time_",
                   srr, "_", species, ".txt")
    if (!file.exists(file)) {
      write(paste0("failure opening file: ", file), stderr())
      times <- c(times, NA)
    }
    else {
      times <- c(times, usr.bin.time(file))
    }
  }
  names(times) <- mappers
  times
}

get.mems.for.srr <- function(srr, species) {
  mappers <- c("abismal", mappers)
  mems <- c()
  for (i in mappers) {
    file <- paste0(BENCHMARK_PATH, "/", i, "/snakemake_time_",
                   srr, "_", species, ".txt")

    if (!file.exists(file))  mems <- c(mems, NA)
    else mems <- c(mems, usr.bin.mem(file))

  }
  names(mems) <- mappers
  mems
}

add.abismal.index.time <- function(times) {
  cur.abismal.times <- unlist(times[,"abismal"])
  species <- gsub("^.*_","",names(cur.abismal.times))
  index.times <- list()

  # get indexing times for all species
  ref.path <- "~/panasas/ref_genomes"
  for (i in unique(species)) {
    line <- readLines(paste0(ref.path, "/", i, "_usr_bin_time.txt"))[5]
    line <- gsub("^.* ", "", line)
    minutes <- as.numeric(gsub(":.*$","",line))
    seconds <- as.numeric(gsub("^.*:","",line))
    index.times[[i]] <- 60*minutes + seconds
  }
  wt(unlist(index.times), "results/tables/supp_table_4_index_times.tsv")

  # add read indexing times
  for (i in 1:length(cur.abismal.times))
    cur.abismal.times[i] <- cur.abismal.times[i] + index.times[[species[i]]]
  
  times[,"abismal"] <- cur.abismal.times
  times
}

make.times <- function(tbl) {
  ans <- NULL
  srrs <- tbl$srr
  ids <- tbl$id
  species <- tbl$species

  for (i in 1:length(srrs))
    ans <- rbind(ans, get.times.for.srr(srrs[i], species[i]))

  # add the "num reads" column from tbl which will be the x axis
  ans <- cbind(ans, tbl$total.reads)
  ans <- cbind(ans, tbl$protocol)
  rownames(ans) <- paste0(ids, "_", species)
  colnames(ans) <- c("abismal", mappers, "total.reads", "protocol")

  ans <- add.abismal.index.time(ans)
  ans
}

make.mems <- function(tbl) {
  ans <- NULL
  srrs <- tbl$srr
  ids <- tbl$id
  species <- tbl$species
  if(length(srrs) != length(species))
    stop("unequal number of srrs and species")
  for (i in 1:length(srrs))
    ans <- rbind(ans, get.mems.for.srr(srrs[i], species[i]))

  rownames(ans) <- paste0(ids, "_", species)
  ans
}


# rounds all numeric elements of a table to 3 sig digitis
round.df <- function(x, digits = 3) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

sv <- function(tbl, file = "dataset_comparison.tsv") {
  write.table(tbl, file = file, quote=F, sep="\t", row.names=T, col.names=T)
}

#############################################################
#### PLOTS
#############################################################

# plots mapping %
plot.percentage.line <- function(tbl,
                     mapper = mappers,
                     what,
                     axis.label,
                     start.at = 0,
                     fig.legend = "",
                     width = 8, height = 8,
                     which.keep = rep(TRUE, nrow(tbl)),
                     zero.to.hundred = F) {

  which.keep <- ifelse(is.na(which.keep), FALSE, which.keep)
  tbl <- tbl[which.keep,]
  x <- tbl[, paste0(what, ".", mapper)]
  y <- tbl[, paste0(what, ".abismal")]

  species <- factor(unlist(tbl$species[which.keep]))

  # axis limits
  lim.min <- ifelse(zero.to.hundred, start.at,
                    min(c(min(unlist(x), na.rm = T),
                          min(unlist(y), na.rm = T))))
  lim.max <- ifelse(zero.to.hundred, 100,
                    max(c(max(unlist(x), na.rm = T),
                          max(unlist(y), na.rm = T))))
  plot(NA,
       xlim = c(lim.min, lim.max),
       ylim = c(lim.min, lim.max),
       ylab = paste0("abismal ", axis.label),
       xlab = paste0(mapper," ", axis.label),
       main = mapper,
       cex.axis = 1.5,
       cex.lab = 1.5,
       cex.main=1.5)
  points(x, y,
         pch = 19,
         col = as.integer(species));
  legend("bottomright",
         legend = levels(species),
         col = 1:length(unique(species)),
         pch = 19,
         cex = 1)
  #legend("bottomright",
  #       legend = c("pos", "neg"),
  #       pch = c(19, 1), col = "black")
  abline(coef =c(0,1), col = "red", lty = 2)
  fig_label(fig.legend, cex=2)
}

plot.collision <- function(species) {
  df <- read.table(paste0("hit_statistics/",species,".tsv"),
                   header=T, row.names=NULL)
  plot(df$k,
       df$p2.p3,
       pch = 19,
       cex = .5,
       ylim = c(0,1.1),
       xlab = "number of bits",
       ylab = "two- / three-letter hit ratio",
       main = species)
  points(df$k, df$p2_theo.p3_theo, col = "red", pch = 19, cex = 0.5)
  abline(h = 1, col = "red", lty = 2)
  legend("topright", pch = c(19,19), col=c("black","red"),
         legend = c("real","theoretical"),
         bg = "white")
}

readable.by.barplot <- function(x) {
  rn <- rownames(x)
  cn <- colnames(x)
  x <- as.numeric(unlist(x))
  x <- matrix(x, nrow = length(rn))
  rownames(x) <- rn
  colnames(x) <- cn
  x
}

plot.map <- function(tbl, protocols, main, what, ylab, fig.label) {
  tbl <- tbl[unlist(tbl$protocol) %in% protocols,]
  tbl <- get.data(tbl, what)
  tbl[is.na(tbl)] <- 0

  tbl <- t(readable.by.barplot(tbl))
  # remove species name
  colnames(tbl) <- gsub(" .*$", "", colnames(tbl))
  # put a space instead of underscore for style purposes
  colnames(tbl) <- gsub("_", " ", colnames(tbl))

  palette(plot.good.colors)
  num.mappers <- 1 + length(mappers)
  cols <- seq(1, num.mappers)
  barplot(tbl, beside = T, horiz = T, xlab = ylab,
          col = cols, main = main, xlim = c(0,1),
          cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.5, cex.main = 1.5)

  fig_label(fig.label, cex=2)
}

plot.mems <- function(mems) {
  # only plot what was successfully run in all mappers
  has.na <- apply(mems, 1, function(x) any(is.na(x)))
  mems <- mems[!has.na,]

  b <- reshape::melt(mems / 1024)
  b$value <- as.numeric(as.character(b$value))
  boxplot(b$value ~ b$X2, na.rm = T,
          xlab = "mapper",
          ylab = "memory used (GB)",
          main = "maximum memory usage benchmark",
          col = c(6,7,8,9),
          cex.axis = 1.5,
          cex.lab = 1.5,
          cex.main=1.5)
  abline(h = 16, col = "red", lty = 2)
}


get.data <- function(tbl, what, erase.prefix = F) {
	a <- tbl[, grepl(what, colnames(tbl))]

  if (erase.prefix)
    colnames(a) <- gsub("^.*\\.", "", colnames(a))
  a
}

##################################
# MAKE FIGURES
##################################
make.accuracy.figure <- function(tbl) {
  pdf("results/figures/accuracy.pdf", width = 12, height = 5)
  layout(mat = matrix(c(1,2,3,4,5,6,6,6,6,6), nrow = 2, byrow = T),
         heights = c(.9, .1))

  par(mar = c(1.5, 1.5, 1.5, 1.5), family = "Times")
  plot.map(tbl, "wgbs_single", "Traditional (SE)",
           "accuracy", "accuracy", "A")
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "accuracy", "accuracy", "B")
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "accuracy", "accuracy", "C")
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "map.total", "% concordant pairs mapped", "D")
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "map.total", "% concordant pairs mapped", "E")

  # common legend to all plots
  num.mappers <- 1 + length(mappers)
  cols <- seq(1, num.mappers)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", legend = c("abismal", mappers), cex = 1.5, pt.cex = 3, pch = 22,
         pt.bg = cols, col = rep("black", num.mappers), horiz = T, bty = "n")
  dev.off()
}

# aux function to turn times to reads per sec
time.to.reads.per.sec <- function(times) {
  tot.reads <- unlist(times[, "total.reads"])
  for(i in c("abismal", mappers)) {
    col <- unlist(times[,i])
    times[,i] <- (tot.reads / col) * (3600 / 1000000)
  }
  times
}

get.top.times <- function(times, top.n) {
  # paired datasets have twice as many reads
  protocols <- unlist(times[,"protocol"])
  paired <- which(protocols %in% c("wgbs_paired", "wgbs_rpbat"))
  for (i in paired)
    times[i, "total.reads"][[1]] <- 2*times[i, "total.reads"][[1]]

  times[order(-unlist(times[,"total.reads"]))[1:top.n],]
}

make.resources.figure <- function(times, mems) {
  times <- time.to.reads.per.sec(times)
  
  # get the top 10 largest datasets
  # times <- get.top.times(times, 15)

  # we'll use this for mem distribution
  paired.datasets <- which(unlist(times[,"protocol"]=="wgbs_paired"))

  # convert to readable form for barplot
  times <- t(readable.by.barplot(times[,1:(length(mappers) + 1)]))

  # convert mems to GB
  mems <- mems / (1024)

  # for plotting purposes, cap memory values to 32GB
  max.mem.to.report <- 16
  for (i in 1:nrow(mems))
    for (j in 1:ncol(mems))
      mems[i,j] = min(mems[i,j], max.mem.to.report)

  cn <- colnames(times)
  cn <- gsub("_[^_]+$","",cn)
  cn <- gsub("_", " ", cn)

  palette(plot.good.colors)
  cols <- seq(1, 1 + length(mappers))

  pdf("results/figures/resources.pdf", width = 12, height = 4)
  layout(matrix(c(1,2,3,4,5,16,16,
                  6,7,8,9,10,16,16,
                  11,12,13,14,15,16,16), nrow = 3, byrow = T))
  par(family = "Times")
  par(mar = c(1.5, 3, 1.5, 3))
  for (i in 1:ncol(times)) {
    barplot(times[,i], beside = T, horiz = F, col = cols, names.arg = NA,
            main = cn[i], ylab = "reads per hour (M)", cex.lab = 1.5,
            font.main = 1, cex.axis = 1.5)

    if (i == 1) fig_label("A", cex = 2)
  }

  boxplot(mems[paired.datasets,], col = cols, horizontal = F, las = 2,
          ylab = "Memory (GB)", cex.axis = 1.5, cex.lab = 1.5)

  fig_label("B", cex = 2)
  dev.off()
}

wt <- function(x, file) {
  write.table(as.matrix(x), file, quote = F, sep="\t", row.names = T, col.names = T)
}

################# MAIN ################
# write benchmark table
tbl <- make.table(datasets)

# resources
times <- make.times(tbl)
mems <- make.mems(tbl)
if (plot) {
  # Supp table 1 is accuracy
  accuracy <- get.data(tbl, "accuracy")
  conc.pairs <- get.data(tbl, "map.total")
  truth.table <- get.data(tbl, "truth", F)
  wt(cbind(accuracy, conc.pairs, truth.table),
     "results/tables/supp_table_1_accuracy.tsv")

  # Supp table 3 is bsrate and coverage
  bsrate <- get.data(tbl, "bs.err")
  coverage <- get.data(tbl, "coverage.cpg")
  depth <- get.data(tbl, "cpg.depth.total")
  meth.cytosine <- get.data(tbl, "meth.cytosine")
  meth.cpg <- get.data(tbl, "meth.cpg")
  wt(cbind(bsrate, coverage, depth, meth.cytosine, meth.cpg),
     "results/tables/supp_table_2_bsrate.tsv")

  # Supp table 4 are the times
  mm <- c("abismal", mappers)
  a <- times[,mm]
  b <- mems[,mm]
  colnames(a) <- paste0("time_", colnames(a))
  colnames(b) <- paste0("mem_", colnames(b))
  wt(cbind(a,b), "results/tables/supp_table_3_resources.tsv")

  ########### FIGURES ###########
  primary <- unlist(datasets$primary) == "yes"
  tbl <- tbl[primary,]
  times <- times[primary,]
  mems <- mems[primary,]

  # for the main figure we filter only the primary tests
  datasets <- datasets[as.character(datasets$primary_test) == "yes",]

  # accuracy
  make.accuracy.figure(tbl)

  # resources
  make.resources.figure(times, mems)
}

