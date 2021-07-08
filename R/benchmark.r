source("R/fig_label.r")

BENCHMARK_PATH =
"results/benchmark-files"

datasets <- read.table("metadata/tests.txt",
                       header = 1,row.names = NULL,stringsAsFactors = F)
mappers <- c("bismark", "bsmap", "walt")

##########################
# Aux functions
##########################

# Write a nice log
write.log <- function(...) {
  write(paste0("[", Sys.time(), "] ", ...), stderr())
}

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
  if(!(length(lines) %in% c(40, 42, 43)))
    stop(paste0("bad bismark file: ", path, " number of lines: ", length(lines)))

  tot.line <- 7
  uniq.line <- 8
  unmapped.line <- 10

  if (length(lines) == 43) {
    tot.line <- 10
    uniq.line <- 11
    unmapped.line <- 13
  }
  
  ans <- list()

  # mapping statistics
  ans$total.reads <- as.integer(gsub("^.*\t","",lines[tot.line]))
  ans$mapped.unique <- as.integer(gsub("^.*\t","",lines[uniq.line]))

  ans$unmapped <- as.integer(gsub("^.*\t","",lines[unmapped.line]))
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
# parse downstream analysis outputs
#################################################

parse.samstats <- function(samstats.path) {
  dat <- readLines(samstats.path)
  dat <- dat[grepl("error rate:", dat)]
  dat <- strsplit(dat, "[ \t]")[[1]][4]

  as.numeric(dat)
}

parse.samstats.minimap2 <- function(samstats.path) {
  if (!file.exists(samstats.path))
    return (list(sam.err.rate = NA, map.total = NA, map.unique = NA))
  dat <- readLines(samstats.path)

  ans <- NULL
  tot.reads <- NULL
  mapped.reads <- NULL
  unmapped.reads <- NULL
  supp.reads <- NULL
  for (i in 1:length(dat)) {
    if (grepl("^SN", dat[i])) {
      explode <- strsplit(dat[i], "[ \t]")[[1]]
      if (grepl("error rate:", dat[i]))
        ans$sam.err.rate <- as.numeric(explode[4])
      else if (grepl("reads unmapped:", dat[i]))
        unmapped.reads <- as.numeric(explode[4])
      else if (grepl("reads mapped:", dat[i]))
        mapped.reads <- as.numeric(explode[4])
      else if (grepl("supplementary alignments:", dat[i]))
        supp.reads <- as.numeric(explode[4])
      else if (grepl("raw total sequences:", dat[i]))
        tot.reads <- as.numeric(explode[5])
    }
  }
  ans$map.total <- (mapped.reads)/tot.reads
  ans$map.unique <- (mapped.reads - supp.reads)/tot.reads
  return(ans)
}


# get average error rate across all bases in bsrate file
parse.bsrate <- function(bsrate.path) {
  if(!file.exists(bsrate.path))
    stop(paste0("bad bsrate file: ", bsrate.path))

  bs.data <- readLines(bsrate.path)
  ans <- list()
  ans$bs.conv.tot <- as.numeric(strsplit(bs.data[1], " ")[[1]][5])

  # 34 = 30th base
  the.line <- strsplit(bs.data[34], "[ \t]")[[1]]
  ans$bs.conv.rand <- (as.numeric(the.line[3]) + as.numeric(the.line[6]))/
                      (as.numeric(the.line[2]) + as.numeric(the.line[5]))
  ans$bs.pos <- as.numeric(strsplit(bs.data[2], "[ \t]")[[1]][5])
  ans$bs.neg <- as.numeric(strsplit(bs.data[3], "[ \t]")[[1]][5])

  tot.bases <- 0
  err.bases <- 0

  # GS: line 4 is the headers
  for (i in 5:length(bs.data)) {
    the.line <- strsplit(bs.data[i], "\t")[[1]]
    tot.bases <- tot.bases + as.integer(the.line[12])
    err.bases <- err.bases + as.integer(the.line[11])
  }
  ans$err.rate <- err.bases/tot.bases

  return(ans)
}

# get useful information from levels
parse.levels <- function(levels.path) {
  if(!file.exists(levels.path))
    stop(paste0("bad levels file: ", levels.path))

  suppressWarnings({
    lev <- yaml::yaml.load_file(levels.path)
  })

  ans <- list()
  ans$mean.cytosine <- lev$cytosines$mean_meth
  ans$weigh.cytosine <- lev$cytosines$mean_meth_weighted
  ans$frac.cytosine <- lev$cytosines$fractional_meth

  ans$mean.cpg<- lev$cpg$mean_meth
  ans$weigh.cpg<- lev$cpg$mean_meth_weighted
  ans$frac.cpg<- lev$cpg$fractional_meth

  ans$mean.cpg.sym <- lev$cpg_symmetric$mean_meth
  ans$weigh.cpg.sym <- lev$cpg_symmetric$mean_meth_weighted
  ans$frac.cpg.sym <- lev$cpg_symmetric$fractional_meth

  ans$mut.cytosine       <- lev$cytosines$mutations
  ans$mut.cpg            <- lev$cpg$mutations
  ans$mut.cpg.sym        <- lev$cpg_symmetric$mutations

  ans$frac.cytosine      <- lev$cytosines$fractional_meth
  ans$frac.cpg           <- lev$cpg$fractional_meth
  ans$frac.cpg.sym       <- lev$cpg_symmetric$fractional_meth

  ans$coverage.cytosine  <- lev$cytosines$sites_covered_fraction
  ans$coverage.cpg       <- lev$cpg$sites_covered_fraction
  ans$coverage.cpg.sym   <- lev$cpg_symmetric$sites_covered_fraction

  ans$depth.cytosine     <- lev$cytosines$mean_depth_covered
  ans$depth.cpg          <- lev$cpg$mean_depth_covered
  ans$depth.cpg.sym      <- lev$cpg_symmetric$mean_depth_covered
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
    map.total = NA,
    map.unique = NA,

    sam.err.rate = NA,

    bs.conv.tot = NA,
    bs.conv.rand = NA,
    bs.pos = NA,
    bs.neg = NA,
    bs.err.rate = NA,

    c.meth.mean = NA,
    c.meth.weigh = NA,
    c.meth.frac = NA,

    c.frac = NA,
    c.mut = NA,
    c.cov = NA,
    c.depth = NA,

    cpg.meth.mean = NA,
    cpg.meth.weigh = NA,
    cpg.meth.frac = NA,

    cpg.frac = NA,
    cpg.mut = NA,
    cpg.cov = NA,
    cpg.depth = NA,
    
    cpg.sym.meth.mean = NA,
    cpg.sym.meth.weigh = NA,
    cpg.sym.meth.frac = NA,

    cpg.sym.frac = NA,
    cpg.sym.mut = NA,
    cpg.sym.cov = NA,
    cpg.sym.depth = NA,
    
    row.names = paste0(srr,"_",species)
  )
  # path to outputs directory
  base.path <- paste0(BENCHMARK_PATH, "/", mapper, "/")

  # fetch mapping data
  path <- paste0(base.path, srr, "_", species, mapstats.suffix)
  tryCatch({
    mapdata <- f(path)
    ans$total.reads = mapdata$total.reads
    ans$map.total = mapdata$pct.total
    ans$map.unique = mapdata$pct.unique
  }, error = function(e){
    write(paste0("failure opening file: ", path), stderr());
  })

  # samtools stats statistics
  # bsrate statistics
  path <- paste0(base.path, srr, "_", species, ".samstats")
  tryCatch({
    ans$sam.err.rate = parse.samstats(path)
  }, error = function(e) {
    write(paste0("failure opening file: ", path), stderr());
  })


  # bsrate statistics
  path <- paste0(base.path, srr, "_", species, ".bsrate")
  tryCatch({
    bsdata <- parse.bsrate(path)
    ans$bs.conv.tot = bsdata$bs.conv.tot
    ans$bs.conv.rand = bsdata$bs.conv.rand
    ans$bs.pos <- bsdata$bs.pos
    ans$bs.neg <- bsdata$bs.neg
    ans$bs.err.rate <- bsdata$err.rate
  }, error = function(e) {
    write(paste0("failure opening file: ", path), stderr());
  })

  # levels statistics
  path <- paste0(base.path, srr, "_", species, ".levels")
  tryCatch({
    levelsdata <- parse.levels(path)
    ans$c.meth.mean <- levelsdata$mean.cytosine
    ans$c.meth.weigh <- levelsdata$weigh.cytosine
    ans$c.meth.frac <- levelsdata$frac.cytosine

    ans$c.frac <- levelsdata$frac.cytosine
    ans$c.mut <- levelsdata$mut.cytosine
    ans$c.cov <- levelsdata$coverage.cytosine
    ans$c.depth <- levelsdata$depth.cytosine

    ans$cpg.meth.mean <- levelsdata$mean.cpg
    ans$cpg.meth.weigh <- levelsdata$weigh.cpg
    ans$cpg.meth.frac <- levelsdata$frac.cpg

    ans$cpg.frac <- levelsdata$frac.cpg
    ans$cpg.mut <- levelsdata$mut.cpg
    ans$cpg.cov <- levelsdata$coverage.cpg
    ans$cpg.depth <- levelsdata$depth.cpg

    ans$cpg.sym.meth.mean <- levelsdata$mean.cpg.sym
    ans$cpg.sym.meth.weigh <- levelsdata$weigh.cpg.sym
    ans$cpg.sym.meth.frac <- levelsdata$frac.cpg.sym

    ans$cpg.sym.frac <- levelsdata$frac.cpg.sym
    ans$cpg.sym.mut <- levelsdata$mut.cpg.sym
    ans$cpg.sym.cov <- levelsdata$coverage.cpg.sym
    ans$cpg.sym.depth <- levelsdata$depth.cpg.sym

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
        write(paste0("[!!!] inconsistent # of reads between bismark and abismal: ",
                  srr, ". bismark = ", k.data$total.reads, ",
                  abismal = ", tot.reads), stderr())
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
    write.log("...reading row ", datasets[i,1])
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

# get the elapsed time value from /get.time.from.snakemake -v
# given the file name as input
get.time.from.snakemake <- function(filename) {
  # GS: use readlines here because of mappers
  # that say "1 day, ..." on it
  line <- readLines(filename)[2]
  time <- strsplit(line, "\t")[[1]][1]
  as.numeric(time)
}

# in gigabytes
get.mem.from.snakemake <- function(filename) {
  line <- readLines(filename)[2]
  mem <- strsplit(line, "\t")[[1]][3]
  as.numeric(mem)
}

get.times.for.srr <- function(srr, species) {
  mappers <- c("abismal", mappers, "minimap2")
  times <- c()
  for (i in mappers) {
    file <- paste0(BENCHMARK_PATH, "/", i, "/snakemake_time_",
                   srr, "_", species, ".txt")
    if (!file.exists(file)) {
      if (i != "minimap2")
        write(paste0("failure opening file: ", file), stderr())
      times <- c(times, NA)
    }
    else {
      times <- c(times, get.time.from.snakemake(file))
    }
  }
  names(times) <- mappers
  times
}

get.mems.for.srr <- function(srr, species) {
  mappers <- c("abismal", mappers, "minimap2")
  mems <- c()
  for (i in mappers) {
    file <- paste0(BENCHMARK_PATH, "/", i, "/snakemake_time_",
                   srr, "_", species, ".txt")

    if (!file.exists(file))  mems <- c(mems, NA)
    else mems <- c(mems, get.mem.from.snakemake(file))

  }
  names(mems) <- mappers
  mems
}

make.index.times <- function(times) {
  cur.abismal.times <- unlist(times[,"abismal"])
  species <- gsub("^.*_","",names(cur.abismal.times))
  index.times <- list()

  # get indexing times for all species
  ref.path <- "results/indexing_times"
  for (i in unique(species)) {
    line <- readLines(paste0(ref.path, "/", i, "_usr_bin_time.txt"))[5]
    line <- gsub("^.* ", "", line)
    minutes <- as.numeric(gsub(":.*$","",line))
    seconds <- as.numeric(gsub("^.*:","",line))
    index.times[[i]] <- 60*minutes + seconds
  }
  return(do.call(c, index.times))
}

make.times <- function(tbl) {
  ans <- NULL
  ans.two <- NULL
  srrs <- tbl$srr
  ids <- tbl$id
  species <- tbl$species

  for (i in 1:length(srrs))
    ans <- rbind(ans, get.times.for.srr(srrs[i], species[i]))

  # add the "num reads" column from tbl which will be the x axis
  ans.two <- tbl$total.reads
  ans.two <- cbind(ans.two, tbl$protocol)
  
  # convert to a format that can be manipulated
   
  rownames(ans) <- rownames(ans.two) <- paste0(ids, "_", species)
  colnames(ans) <- c("abismal", mappers, "minimap2")
  colnames(ans.two) <- c("total.reads", "protocol")
  list(times = ans, meta = ans.two)
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

plot.map <- function(tbl, protocols, main, what, ylab, fig.label, scale.x = T, pct = T) {
  tbl <- tbl[unlist(tbl$protocol) %in% protocols,]
  tbl <- get.data(tbl, what)
  num.mappers <- ncol(tbl)
  tbl[is.na(tbl)] <- 0
  tbl <- t(readable.by.barplot(tbl))

  if (pct)
    tbl <- tbl*100;

  # remove species name
  colnames(tbl) <- gsub(" .*$", "", colnames(tbl))

  # put a space instead of underscore for style purposes
  colnames(tbl) <- gsub("_", " ", colnames(tbl))

  palette(plot.good.colors)
  cols <- seq(1, num.mappers)

  # set the x-axis scale
  max.val <- ifelse(pct, 100, 1)
  if (!scale.x)
    max.val <- max(as.vector(tbl))

  barplot(tbl, beside = T, horiz = T, xlab = ylab,
          col = cols, main = main, xlim = c(0, max.val),
          cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.5, cex.main = 1.5)

  fig_label(fig.label, cex = 2)
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
  
  nr <- nrow(a)
  nc <- ncol(a)
  rn <- rownames(a)
  cn <- colnames(a)
  a <- unlist(a)
  a <- matrix(a, nrow = nr, ncol = nc)
  rownames(a) <- rn
  colnames(a) <- cn

  data.frame(a)
}

##################################
# MAKE FIGURES
##################################
make.accuracy.figure <- function(sim.tbl, tbl, minimap2) {
  pdf("results/figures/accuracy.pdf", width = 15, height = 6)
  layout(mat = matrix(c(8,8,8,1,2,3,4,5,6,
                        9,9,9,1,2,3,4,5,6,
                        7,7,7,7,7,7,7,7,7), nrow = 3, byrow = T),
         heights = c(.425, .425, .15))

  tbl.w.min2 <- cbind(tbl, minimap2)
  tbl.w.min2 <- tbl.w.min2[tbl$protocol == "wgbs_single",]

  par(mar = rep(2, 4), family = "Times")
  plot.map(tbl.w.min2, "wgbs_single", "Traditional (SE)",
           "map.total", "% reads mapped", "B")
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "map.total", "% pairs mapped", "")
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "map.total", "% pairs mapped", "")
  plot.map(tbl.w.min2, "wgbs_single", "Traditional (SE)",
           "sam.err.rate", "error rate (%)", "C", scale.x = F)
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "sam.err.rate", "error rate (%)", "", scale.x = F)
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "sam.err.rate", "error rate (%)", "", scale.x = F)

  # common legend to all plots
  num.mappers <- 5
  cols <- seq(1, num.mappers)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", legend = c("abismal", mappers, "minimap2"), cex = 1.5, pt.cex = 3, pch = 22,
         pt.bg = cols, col = rep("black", num.mappers + 1), horiz = T, bty = "n")

  ############### simulation figure #############
  the.lengths <- sim.tbl$lengths

  # accuracy
  sens <- sim.tbl$err_3$sensitivity
  # hack to put abismal up front
  sens <- sens[, seq(ncol(sens), 1, -1)]

  the.xlim <- c(min(the.lengths), max(the.lengths))
  the.ylim <- c(min(sens, na.rm=T), max(sens, na.rm = T))
  plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
       xlab = "read length", ylab = "correct reads / total reads", main = "sensitivity",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  for (i in 1:ncol(sens))
    points(the.lengths, sens[,i], type = 'o', pch = 22, cex = 2, bg = ncol(sens) - i + 1)
  fig_label("A", cex = 2)

  # specificity
  spec <- sim.tbl$err_3$specificity

  # hack to put abismal up front
  spec <- spec[, seq(ncol(spec), 1, -1)]

  the.ylim <- c(min(spec, na.rm=T), max(spec, na.rm = T))
  plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
       xlab = "read length", ylab = "correct reads / reported reads", main = "specificity",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  for (i in 1:ncol(spec))
    points(the.lengths, spec[,i], type = 'o', pch = 22, cex = 2, bg = ncol(spec) - i + 1)

  ########### FINISHED PLOTTING FIGURE ###########
  dev.off()
}

make.sim.comparison.figure <- function(sim.tbl) {
  pdf("results/figures/supp_2_simulations.pdf", width = 18, height = 9)
  layout(mat = matrix(c(1,3,5,7,9,11,
                        2,4,6,8,10,12,
                        13,13,13,13,13,13), nrow = 3, byrow = T),
         heights = c(.425, .425, .15))

  the.lengths <- sim.tbl$lengths
  the.errs <- sim.tbl$error.rates

  for (err in the.errs) {
    # accuracy
    sens <- sim.tbl[[paste0("err_",err)]]$sensitivity
    # hack to put abismal up front
    sens <- sens[, seq(ncol(sens), 1, -1)]

    the.xlim <- c(min(the.lengths), max(the.lengths))
    the.ylim <- c(min(sens, na.rm=T), max(sens, na.rm = T))
    plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
         xlab = "read length", ylab = "correct reads / total reads",
         main = paste0("sensitivity (error = ", err, "%)"),
         cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
    for (i in 1:ncol(sens))
      points(the.lengths, sens[,i], type = 'o', pch = 22, cex = 2, bg = ncol(sens) - i + 1)

    # specificity
    spec <- sim.tbl[[paste0("err_",err)]]$specificity

    # hack to put abismal up front
    spec <- spec[, seq(ncol(spec), 1, -1)]

    the.ylim <- c(min(spec, na.rm=T), max(spec, na.rm = T))
    plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
         xlab = "read length", ylab = "correct reads / reported reads",
         main = paste0("specificity (error = ", err, "%)"),
         cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
    for (i in 1:ncol(spec))
      points(the.lengths, spec[,i], type = 'o', pch = 22, cex = 2, bg = ncol(spec) - i + 1)
  }
  ########### FINISHED PLOTTING FIGURE ###########
  # common legend to all plots
  num.mappers <- 5
  cols <- seq(1, num.mappers)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", legend = c("abismal", mappers), cex = 1.5, pt.cex = 3, pch = 22,
         pt.bg = cols, col = rep("black", num.mappers + 1), horiz = T, bty = "n")


  dev.off()
}

# aux function to turn times to reads per sec
time.to.reads.per.sec <- function(times) {
  tot.reads <- as.numeric(unlist(times$meta[, "total.reads"]))
  #for(i in c("abismal", mappers, "minimap2")) {
  for(i in c("abismal", mappers, "minimap2")) {
    col <- unlist(times$times[,i])
    times$times[,i] <- (tot.reads / col) * (3600 / 1000000)
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

replace.species.id.with.name <- function(x) {
  x <- gsub("hg38", "H. sapiens", x)
  x <- gsub("mm10", "M. musculus", x)
  x <- gsub("danre11", "D. rerio", x)
  x <- gsub("galgal6", "G. gallus", x)
  x <- gsub("tair10", "A. thaliana", x)
}

make.resources.figure <- function(times, mems, minimap2) {
  times <- time.to.reads.per.sec(times)
  
  # we'll use this for mem distribution
  single.datasets <- which(unlist(times$meta[,"protocol"]=="wgbs_single"))

  # convert to readable form for barplot
  times$times <- t(readable.by.barplot(times$times[,1:(length(mappers) + 2)]))

  # convert mems to GB
  mems <- mems / (1024)

  # for plotting purposes, cap memory values to 32GB
  max.mem.to.report <- 10
  for (i in 1:nrow(mems))
    for (j in 1:ncol(mems))
      mems[i,j] = min(mems[i,j], max.mem.to.report)

  cn <- colnames(times$times)
  cn <- gsub("_", " (", cn)
  cn <- gsub("-", " ", cn)
  cn <- gsub("$", ")", cn)
  cn <- replace.species.id.with.name(cn)

  palette(plot.good.colors)
  cols.paired <- seq(1, 1 + length(mappers))
  cols.single <- seq(1, 2 + length(mappers))

  pdf("results/figures/resources.pdf", width = 12, height = 5)
  layout(matrix(c(1,2,3,4,5,16,16,
                  6,7,8,9,10,16,16,
                  11,12,13,14,15,16,16), nrow = 3, byrow = T))
  par(family = "Times")
  par(mar = c(1.5, 3, 1.5, 3))
  for (i in 1:ncol(times$times)) {
    cols.to.plot <- cols.paired
    if (i %in% single.datasets)
      cols.to.plot <- cols.single

    barplot(times$times[,i], beside = T, horiz = F,
            col = cols.to.plot, names.arg = NA,
            main = cn[i], ylab = "reads per hour (M)", cex.lab = 1.5,
            font.main = 1, cex.axis = 1.5)

    if (i == 1) fig_label("A", cex = 2)
  }

  mems.to.plot <-mems[single.datasets,]
  # transforms mems so that only the species shows up
  rownames(mems.to.plot) <- gsub("^.*_", "", rownames(mems.to.plot))
  mems.to.plot <- t(mems.to.plot)

  # hack to make abismal appear on top
  mems.to.plot <- mems.to.plot[seq(nrow(mems.to.plot), 1, -1),]
  colnames(mems.to.plot) <- replace.species.id.with.name(colnames(mems.to.plot))

  dotchart(mems.to.plot, las = 2, ylab = "Memory (GB)", pch = 21, col = "black",
           bg = rev(cols.single))

  fig_label("B", cex = 2)
  dev.off()
}

wt <- function(x, file) {
  write.table(as.matrix(x), file, quote = F, sep="\t", row.names = T, col.names = T)
}

make.minimap2 <- function(datasets) {
  ans <- NULL
  datasets <- datasets[datasets$protocol == "wgbs_single" & datasets$primary == "yes",]
  for (i in 1:nrow(datasets)) {
    dat <- list()
    dat$map.total <- NA
    dat$map.unique <- NA
    dat$time <- NA
    dat$meim <- NA

    tryCatch({
      srr <- datasets[i,"srr"]
      species <- datasets[i,"species"]
      stats.path <- paste0(BENCHMARK_PATH, "/minimap2/",srr,"_",species,".samstats")
      times.path <- paste0(BENCHMARK_PATH, "/minimap2/snakemake_time_",srr,"_",species, ".txt")
      dat <- parse.samstats.minimap2(stats.path)
      dat$time <- get.time.from.snakemake(times.path)
      dat$mem <- get.mem.from.snakemake(times.path)
    }, error = function(e) {
      write(paste0("cannot open minimap2 file: ", stats.path), stderr())
    })

    ans <- rbind(ans, dat)
  }
  colnames(ans) <- paste0(colnames(ans), ".minimap2")
  rownames(ans) <- paste0(datasets$id , " (", datasets$species, ")")
  ans
}

make.sim.tbl <- function(read.lengths, num.reads = 2000000) {
  cols <- c("abismal", mappers)

  errs <- seq(0, 5)

  ans <- list()
  ans$lengths <- read.lengths
  ans$error.rates <- errs
  for (err in errs) {
    sens <- data.frame(NA, nrow = length(cols), ncol = length(read.lengths))
    spec <- data.frame(NA, nrow = length(cols), ncol = length(read.lengths))
    for (i in 1:length(cols)) {
      for (j in 1:length(read.lengths)) {
        tryCatch({
          the.path <- paste0("simulation_analysis/accuracy-out/",
                             cols[i], "-err-", 
                             err, "-len-",
                             read.lengths[j], ".txt")
          the.data <- yaml::yaml.load_file(the.path)
          sens[j, i] <- the.data$correct_reads / num.reads
          spec[j, i] <- the.data$specificity
        }, error = function(e) {
          sens[j, i] <- NA
          spec[j, i] <- NA
       });
      }
    }
    rownames(sens) <- rownames(spec) <- paste0("length_", read.lengths)
    colnames(sens) <- colnames(spec) <- cols
    ans[[paste0("err_",err)]]$sensitivity <- sens
    ans[[paste0("err_",err)]]$specificity <- spec
  }

  return(ans)
}

################# MAIN ################
do.sim <- T
do.tbl <- T
do.times <- T
do.mems <- T
plot <- T


sim.tbl <- NULL
tbl <- NULL
times <- NULL
mems <- NULL

if (do.sim) {
  write.log("reading simulation output...")
  read.lengths <- seq(from = 50, to = 190, by = 10)
  sim.tbl <- make.sim.tbl(read.lengths)
} else {
  write.log("SIMULATION OUTPUT SKIPPED")
}

if (do.tbl) {
  write.log("reading master table: methpipe/samtools outputs...")
  tbl <- make.table(datasets)
  minimap2 <- make.minimap2(datasets)
} else {
  write.log("MASTER TABLE SKIPPED")
}

if (do.times) {
  write.log("reading snakemake time outputs...")
  times <- make.times(tbl)
} else {
  write.log("TIME DATA SKIPPED")
}

if (do.mems) {
  write.log("reading snakemake memory outputs...")
  mems <- make.mems(tbl)
} else {
  write.log("MEMORY DATA SKIPPED")
}

if (plot) {
  if (do.tbl) {
 
    # Supp table 1 is accuracy
    conc.pairs <- round(100*get.data(tbl, "map.total"), 3)
    sam.err.rate <- round(100*get.data(tbl, "sam.err.rate", F), 3)
    bs.err.rate <- round(100*get.data(tbl, "bs.err.rate", F), 3)
    supp.tab.1 <- cbind(conc.pairs, sam.err.rate, bs.err.rate)

    write.log("SUPP TABLE 1: ACCURACY...")
    wt(supp.tab.1,
       "results/tables/supp_table_1_accuracy.tsv")

    # Supp table 2 is bsrate and coverage
    c.cov <- round(100*get.data(tbl, "c.cov"), 3)
    c.depth <- get.data(tbl, "c.depth")
    c.meth.mean <- round(100*get.data(tbl, "c.meth.mean"), 3)
    c.meth.weigh <- round(100*get.data(tbl, "c.meth.weigh"), 3)
    c.meth.frac <- round(100*get.data(tbl, "c.meth.frac"), 3)
    cpg.cov <- round(100*get.data(tbl, "cpg.cov"), 3)
    cpg.depth <- get.data(tbl, "cpg.depth")
    cpg.meth.mean <- round(100*get.data(tbl, "cpg.meth.mean"), 3)
    cpg.meth.weigh <- round(100*get.data(tbl, "cpg.meth.weigh"), 3)
    cpg.meth.frac <- round(100*get.data(tbl, "cpg.meth.frac"), 3)
    supp.tab.2 <- cbind(c.cov, c.depth, c.meth.mean, c.meth.weigh, c.meth.frac,
             cpg.cov, cpg.depth, cpg.meth.mean, cpg.meth.weigh, cpg.meth.frac)

    write.log("SUPP TABLE 2: METH ACCURACY...")
    wt(supp.tab.2, "results/tables/supp_table_2_downstream.tsv")
  }

  if (do.times & do.mems) {
    # Supp table 3 are the times
    mm <- c("abismal", mappers, "minimap2")
    a <- times$times[,mm]
    b <- mems[,mm]
    colnames(a) <- paste0("time_", colnames(a))
    colnames(b) <- paste0("mem_", colnames(b))
    supp.tab.3 <- cbind(a, b)

    write.log("SUPP TABLE 3: RESOURCES...")
    wt(supp.tab.3, "results/tables/supp_table_3_resources.tsv")

    index.times <- make.index.times(times$times)
    wt(unlist(index.times), "results/tables/supp_table_4_index_times.tsv")
  }
  ########### FIGURES ###########
  if (do.tbl) {
    write.log("FIGURE S2: SIMULATION RESULTS...")
    primary <- unlist(datasets$primary) == "yes"
    tbl <- tbl[primary,]

    # accuracy
    make.accuracy.figure(sim.tbl, tbl, minimap2)
  }
  if (do.sim) {
    write.log("FIGURE 2: ACCURACY...")
    primary <- unlist(datasets$primary) == "yes"
    make.sim.comparison.figure(sim.tbl)
  }
  if (do.times) {
    write.log("FIGURE 3: RESOURCES...")
    times$times <- times$times[primary,]
    times$meta <- times$meta[primary,]
    mems <- mems[primary,]

    # for the main figure we filter only the primary tests
    datasets <- datasets[as.character(datasets$primary_test) == "yes",]

    # resources
    make.resources.figure(times, mems, minimap2)
  }
}
write.log("DONE!")

