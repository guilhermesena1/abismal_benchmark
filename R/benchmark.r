source("R/fig_label.r")
source("R/supplementary.r")

# print warnings as they appear
options(warn = 1)

BENCHMARK_PATH = "results/benchmark-files"

datasets <- read.table("metadata/tests.txt", header = 1,row.names = NULL,stringsAsFactors = F)
mappers <- sort(c("abismal", "bismark", "bsmap", "bwa", "hisat_3n", "walt"))
mappers.rpbat <- sort(c("abismal", "bismark", "bsmap"))

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
  suppressWarnings({
    data <- yaml::yaml.load_file(path)
  })
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

  ans$mapped.unique <- data$mapped$unique
  ans$unmapped <- data$unmapped
  ans$mapped.total <- ans$mapped.unique + data$mapped$ambiguous

  ans <- get.pct(ans)
  ans

}

# get useful information from walt/abismal outputs
parse.abismal <- function(path) {
  if(!file.exists(path))
    stop(paste0(c("bad file: " , path)))
  suppressWarnings({
    data <- yaml::yaml.load_file(path)
  })
  ans <- list()

  # paired ended
  if (!is.null(data$pairs))  {
    data <- data$pairs
    ans$total.reads <- NULL
    if (!is.null(data$total_read_pairs))
      ans$total.reads <- data$total_read_pairs
    else if (!is.null(data$total_pairs))
      ans$total.reads <- data$total_pairs
  }

  # single end
  else {
    ans$total.reads <- data$total_reads
  }

  ans$mapped.unique <- data$mapped$num_unique
  ans$unmapped <- NULL

  if (!is.null(data$num_unmapped))
    ans$unmapped <- data$num_unmapped
  else if (!is.null(data$unmapped))
    ans$unmapped <- data$unmapped

  if (!is.null(data$mapped$ambiguous))
    ans$mapped.total <- ans$mapped.unique + data$mapped$ambiguous
  else if (!is.null(data$mapped$num_ambiguous))
    ans$mapped.total <- ans$mapped.unique + data$mapped$num_ambiguous

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

parse.hisat <- function(path) {
  if(!file.exists(path))
    stop(paste0("bad hisat_3n file: ", path))
  ans <- list()
  data <- readLines(path)

  # works for both PE and SE
  suppressWarnings({
    tot <- as.integer(strsplit(data[1], " ")[[1]])[1]
    uniq <- as.integer(strsplit(data[4], " ")[[1]])[5]
    ambig <- as.integer(strsplit(data[5], " ")[[1]])[5]
  })

  ans$total.reads <- tot
  ans$mapped.total <- ambig + uniq
  ans$mapped.unique <- uniq
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

parse.samstats.custom <- function(samstats.path) {
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
      else if (grepl("reads duplicated:", dat[i]))
        supp.reads <- as.numeric(explode[4])
      else if (grepl("raw total sequences:", dat[i]))
        tot.reads <- as.numeric(explode[5])
    }
  }
  ans$total.reads <- tot.reads
  ans$mapped.total <- mapped.reads
  ans$mapped.unique <- mapped.reads - supp.reads
  ans <- get.pct(ans)

  # GS: hack
  ans$map.total <- ans$pct.total
  ans$map.unique <- ans$pct.unique
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

  if (protocol == "wgbs_rpbat" & !(mapper %in% mappers.rpbat)) {

    colnames(ans) <- paste0(colnames(ans), ".", mapper)
    return(ans)
  }

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

    if (!is.null(w.data$total.reads.walt)) {
      if (is.na(tot.reads))
        tot.reads <- w.data$total.reads.walt

      if (!is.na(tot.reads) & !is.na(w.data$total.reads.walt))
        if (w.data$total.reads.walt != tot.reads)
          stop(paste0("inconsistent # of reads between walt and abismal: ",
                    srr))
      w.data$total.reads.walt <- NULL
    }
  }
  ################# BISMARK #######################
  if ("bismark" %in% mappers) {
    k.data <- get.row.for.mapper(srr, "bismark",
                                 ".mapstats",
                                 parse.bismark, species, protocol)

    if (is.na(tot.reads))
      tot.reads <- k.data$total.reads.bismark

    if (!is.na(tot.reads) & !is.na(k.data$total.reads.bismark))
      if (k.data$total.reads.bismark != tot.reads)
        write(paste0("[!!!] inconsistent # of reads between bismark and abismal: ",
                  srr, ". bismark = ", k.data$total.reads.bismark, ",
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

  if ("bwa" %in% mappers) {
    bw.data <- get.row.for.mapper(srr, "bwa", ".samstats", parse.samstats.custom, species, protocol)
    bw.data$total.reads.bwa <- NULL
  }

  if ("hisat_3n" %in% mappers) {
    hisat.data <- get.row.for.mapper(srr, "hisat_3n", ".mapstats", parse.hisat, species, protocol)
    if (is.na(tot.reads))
      tot.reads <- hisat.data$`total.reads.hisat_3n`

    if (!is.na(tot.reads) & !is.na(hisat.data$`total.reads.hisat_3n`))
      if(hisat.data$`total.reads` != tot.reads) {
        stop(paste0("inconsistent # of reads between hisat_3n and abismal: ",
                  srr))
      }
    hisat.data$`total.reads.hisat_3n` <- NULL
  }

  # flatten information from all mappers
  row <- a.data
  if ("bismark" %in% mappers) row <- cbind(row, k.data)
  if ("bsmap" %in% mappers) row <- cbind(row, p.data)
  if ("walt" %in% mappers) row <- cbind(row, w.data)
  if ("bwa" %in% mappers) row <- cbind(row, bw.data)
  if ("hisat_3n" %in% mappers) row <- cbind(row, hisat.data)
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

get.times.for.srr <- function(srr, species, protocol) {
  times <- c()
  for (i in mappers) {
    file <- paste0(BENCHMARK_PATH, "/", i, "/snakemake_time_",
                   srr, "_", species, ".txt")
    if (!file.exists(file)) {
      if (!(protocol == "wgbs_rpbat" & !(i %in% mappers.rpbat))) {
        write(paste0("failure opening file: ", file), stderr())
      }
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
  protocol <- tbl$protocol

  for (i in 1:length(srrs))
    ans <- rbind(ans, get.times.for.srr(srrs[i], species[i], protocol[i]))

  # add the "num reads" column from tbl which will be the x axis
  ans.two <- tbl$total.reads
  ans.two <- cbind(ans.two, tbl$protocol)
  
  # convert to a format that can be manipulated
   
  rownames(ans) <- rownames(ans.two) <- paste0(ids, "_", species)
  colnames(ans) <- mappers
  colnames(ans.two) <- c("total.reads", "protocol")
  list(times = ans, meta = ans.two)
}



make.mems <- function(tbl) {
  #tbl <- tbl[tbl$protocol =="wgbs_single",]
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

  # only put rpbat for mappers that have rpbat mode
  if ("wgbs_rpbat" %in% protocols)
    tbl <- tbl[, paste0(what, ".", mappers.rpbat)]

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


get.data <- function(tbl, what, erase.prefix = F, do.round = F) {
	a <- tbl[, grepl(what, colnames(tbl))]

  if (erase.prefix)
    colnames(a) <- gsub("^.*\\.", "", colnames(a))
  
  nr <- nrow(a)
  nc <- ncol(a)
  rn <- rownames(a)
  cn <- colnames(a)
  a <- as.numeric(unlist(a))
  a <- matrix(a, nrow = nr, ncol = nc)
  rownames(a) <- rn
  colnames(a) <- cn

  a <- data.frame(a)
  if (do.round)
    a <- round(100*a, 3)
  a
}

##################################
# MAKE FIGURES
##################################
make.accuracy.figure <- function(sim.tbl, tbl) {
  pdf("results/figures/accuracy.pdf", width = 12, height = 6)
  layout(mat = matrix(c(8,8,8,1,2,3,4,5,6,
                        9,9,9,1,2,3,4,5,6,
                        7,7,7,7,7,7,7,7,7), nrow = 3, byrow = T),
         heights = c(.425, .425, .15))

  
  par(mar = rep(2, 4), family = "Times")
  plot.map(tbl, "wgbs_single", "Traditional (SE)",
           "map.total", "% reads mapped", "B")
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "map.total", "% pairs mapped", "")
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "map.total", "% pairs mapped", "")
  plot.map(tbl, "wgbs_single", "Traditional (SE)",
           "sam.err.rate", "error rate (%)", "C", scale.x = F)
  plot.map(tbl, "wgbs_paired", "Traditional (PE)",
           "sam.err.rate", "error rate (%)", "", scale.x = F)
  plot.map(tbl, "wgbs_rpbat", "RPBAT (PE)",
           "sam.err.rate", "error rate (%)", "", scale.x = F)

  # common legend to all plots
  num.mappers <- 6
  cols <- seq(1, num.mappers)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", legend = mappers, cex = 1.5, pt.cex = 3, pch = 22,
         pt.bg = cols, col = rep("black", num.mappers + 2), horiz = T, bty = "n")

  ############### simulation figure #############
  the.lengths <- sim.tbl$pe$lengths

  # accuracy
  sens <- sim.tbl$pe$avg$sensitivity
  sens <- sens[, colnames(sens) != "walt"]

  # hack to put abismal up front
  sens <- sens[, seq(ncol(sens), 1, -1)]

  the.xlim <- c(min(the.lengths), max(the.lengths))
  the.ylim <- c(min(sens), max(sens))
  plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
       xlab = "read length", ylab = "correct reads / total reads", main = "sensitivity",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  for (i in 1:ncol(sens))
    points(the.lengths, sens[,i], type = 'o', pch = 22, cex = 2, bg = ncol(sens) - i + 1)
  fig_label("A", cex = 2)

  # specificity
  spec <- sim.tbl$pe$avg$specificity

  # hack to put abismal up front
  spec <- spec[, seq(ncol(spec), 1, -1)]
  spec <- spec[, colnames(spec) != "walt"]

  the.ylim <- c(min(spec, na.rm=T), max(spec, na.rm = T))
  plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
       xlab = "read length", ylab = "correct reads / reported reads", main = "specificity",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  for (i in 1:ncol(spec))
    points(the.lengths, spec[,i], type = 'o', pch = 22, cex = 2, bg = ncol(spec) - i + 1)

  ########### FINISHED PLOTTING FIGURE ###########
  dev.off()
}

make.sim.comparison.figure <- function(sim.tbl, outfile) {
  pdf(outfile, width = 16, height = 10)
  par(family = "Times")
  layout(mat = matrix(c(1,4,7,10,13,16,
                        2,5,8,11,14,17,
                        3,6,9,12,15,18,
                        19,19,19,19,19,19), nrow = 4, byrow = T),
         heights = c(.3, .3, .3, .1))
  par(mar = c(1,1,1,1))
  palette(plot.good.colors)

  the.lengths <- sim.tbl$lengths
  the.errs <- sim.tbl$error.rates

  for (err in the.errs) {

    ############## sensitivity ###############

    sens <- sim.tbl[[paste0("err_",err)]]$sensitivity
    # hack to put abismal up front
    sens <- sens[, seq(ncol(sens), 1, -1)]

    the.xlim <- c(min(the.lengths, na.rm = T), max(the.lengths, na.rm = T))
    the.ylim <- c(min(sens, na.rm=T), max(sens, na.rm = T))
    plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
         xlab = "read length", ylab = "correct reads / total reads",
         main = paste0("sensitivity (error = ", err, "%)"),
         cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
    for (i in 1:ncol(sens))
      points(the.lengths, sens[,i], type = 'o', pch = 22, cex = 2, bg = ncol(sens) - i + 1)

    ############## specificity #################

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

    ############## specificity #################

    time <- sim.tbl[[paste0("err_",err)]]$time

    # hack to put abismal up front
    time <- time[, seq(ncol(time), 1, -1)]

    the.ylim <- c(min(time, na.rm=T), max(time, na.rm = T))
    plot(1, type = 'n', xlim = the.xlim, ylim = the.ylim,
         xlab = "read length", ylab = "mapping time (s)",
         main = paste0("time (error = ", err, "%)"),
         cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
    for (i in 1:ncol(spec))
      points(the.lengths, time[,i], type = 'o', pch = 22, cex = 2, bg = ncol(time) - i + 1)


  }
  ########### FINISHED PLOTTING FIGURE ###########
  # common legend to all plots
  num.mappers <- length(mappers)
  cols <- seq(1, num.mappers)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("top", legend = mappers, cex = 1.5, pt.cex = 3, pch = 22,
         pt.bg = cols, col = rep("black", num.mappers), horiz = T, bty = "n")


  dev.off()
}

# aux function to turn times to reads per sec
time.to.reads.per.sec <- function(times) {
  tot.reads <- as.numeric(unlist(times$meta[, "total.reads"]))
  for(i in mappers) {
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

make.resources.figure <- function(times, mems) {
  times <- time.to.reads.per.sec(times)
  
  # we'll use this for mem distribution
  single.datasets <- which(unlist(times$meta[,"protocol"]=="wgbs_single"))
  rpbat.datasets <- which(unlist(times$meta[,"protocol"]=="wgbs_rpbat"))

  # convert to readable form for barplot
  times$times <- t(readable.by.barplot(times$times))

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
  cols.paired <- seq(1, length(mappers))
  cols.single <- seq(1, length(mappers))
  cols.rpbat <- seq(length(mappers.rpbat))

  pdf("results/figures/resources.pdf", width = 12, height = 5)
  layout(matrix(c(1,2,3,4,5,16,16,
                  6,7,8,9,10,16,16,
                  11,12,13,14,15,16,16), nrow = 3, byrow = T))
  par(family = "Times")
  par(mar = c(1.5, 3, 1.5, 3))
  for (i in 1:ncol(times$times)) {
    to.plot <- times$times[,i]
    cols.to.plot <- cols.paired
    if (i %in% single.datasets)
      cols.to.plot <- cols.single
    
    if (i %in% rpbat.datasets) {
      cols.to.plot <- cols.rpbat
      to.plot <- to.plot[mappers.rpbat]
    }

    barplot(to.plot, beside = T, horiz = F,
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

  dotchart(mems.to.plot, las = 2, ylab = "Memory (GB)", pch = 22, col = "black",
           bg = rev(cols.single))

  fig_label("B", cex = 2)
  dev.off()
}

wt <- function(x, file) {
  write.table(as.matrix(x), file, quote = F, sep="\t", row.names = T, col.names = T)
}


make.sim.tbl <- function(read.lengths, num.reads = 2000000) {
  cols <- mappers
  cols.r <- mappers.rpbat

  errs <- seq(0, 5)

  ans <- list()
  ans$pe <- ans$rpbat <- list()
  ans$pe$lengths <- ans$rpbat$lengths <- read.lengths
  ans$pe$error.rates <- ans$rpbat$error.rates <- errs
  for (err in errs) {
    sens <- data.frame(NA, ncol = length(cols), nrow = length(read.lengths))
    spec <- data.frame(NA, ncol = length(cols), nrow = length(read.lengths))
    f1   <- data.frame(NA, ncol = length(cols), nrow = length(read.lengths))
    time <- data.frame(NA, ncol = length(cols), nrow = length(read.lengths))

    sens.r <- data.frame(NA, ncol = length(cols.r), nrow = length(read.lengths))
    spec.r <- data.frame(NA, ncol = length(cols.r), nrow = length(read.lengths))
    f1.r   <- data.frame(NA, ncol = length(cols.r), nrow = length(read.lengths))
    time.r <- data.frame(NA, ncol = length(cols.r), nrow = length(read.lengths))
 
    for (i in 1:length(cols)) {
      for (j in 1:length(read.lengths)) {
        index.time <- 0
        index.time.r <- 0
        tryCatch({
          the.path <- paste0("simulation_analysis/accuracy-out-pe/",
                             cols[i], "-err-",  err, "-len-", read.lengths[j], ".txt")
          the.time <- paste0("simulation_analysis/times-pe/",
                             cols[i], "-err-",  err, "-len-", read.lengths[j], ".txt")

          if (cols[i] %in% mappers.rpbat) {
            the.path.r <- paste0("simulation_analysis/accuracy-out-rpbat/",
                               cols[i], "-err-",  err, "-len-", read.lengths[j], ".txt")
            the.time.r <- paste0("simulation_analysis/times-rpbat/",
                               cols[i], "-err-",  err, "-len-", read.lengths[j], ".txt")
          }
          the.data <- yaml::yaml.load_file(the.path)
          the.data.r <- yaml::yaml.load_file(the.path.r)

          sens[j, i] <- the.data$correct_reads / num.reads
          spec[j, i] <- the.data$specificity
          f1[j, i] <- the.data$F1
          time[j, i] <- get.time.from.snakemake(the.time) - index.time

          if (cols[i] %in% mappers.rpbat) {
            sens.r[j, i] <- the.data.r$correct_reads / num.reads
            spec.r[j, i] <- the.data.r$specificity
            f1.r[j, i] <- the.data.r$F1
            time.r[j, i] <- get.time.from.snakemake(the.time.r) - index.time.r
          }


        }, error = function(e) {
          sens[j, i]   <- spec[j, i]   <- f1[j,i]   <- time[j, i] <- 
          sens.r[j, i] <- spec.r[j, i] <- f1.r[j,i] <- time.r[j, i] <- NA
       });
      }
    }
    rownames(sens) <- rownames(spec) <- rownames(f1) <- rownames(time) <-
      rownames(sens.r) <- rownames(spec.r) <- rownames(f1.r) <- rownames(time.r) <- paste0("length_", read.lengths)
    colnames(sens) <- colnames(spec) <- colnames(f1) <- colnames(time) <- cols
    colnames(sens.r) <- colnames(spec.r) <- colnames(f1.r) <- colnames(time.r) <- cols.r

    ans$pe[[paste0("err_",err)]]$sensitivity <- sens
    ans$pe[[paste0("err_",err)]]$specificity <- spec
    ans$pe[[paste0("err_",err)]]$f1 <- f1
    ans$pe[[paste0("err_",err)]]$time <- time

    ans$rpbat[[paste0("err_",err)]]$sensitivity <- sens.r
    ans$rpbat[[paste0("err_",err)]]$specificity <- spec.r
    ans$rpbat[[paste0("err_",err)]]$f1 <- f1.r
    ans$rpbat[[paste0("err_",err)]]$time <- time.r
  }

  sens.pe <- ans$pe$err_0$sensitivity
  spec.pe <- ans$pe$err_0$specificity
  f1.pe <- ans$pe$err_0$f1
  
  sens.rpbat <- ans$rpbat$err_0$sensitivity
  spec.rpbat <- ans$rpbat$err_0$specificity
  f1.rpbat <- ans$rpbat$err_0$f1

  for (err in errs)
    if (err != 0) {
      sens.pe <- sens.pe + ans$pe[[paste0("err_",err)]]$sensitivity
      spec.pe <- spec.pe + ans$pe[[paste0("err_",err)]]$specificity
      f1.pe <- f1.pe + ans$pe[[paste0("err_",err)]]$f1

      sens.rpbat <- sens.rpbat + ans$rpbat[[paste0("err_",err)]]$sensitivity
      spec.rpbat <- spec.rpbat + ans$rpbat[[paste0("err_",err)]]$specificity
      f1.rpbat <- f1.rpbat + ans$rpbat[[paste0("err_",err)]]$f1
    }
  sens.pe <- sens.pe/length(errs)
  spec.pe <- spec.pe/length(errs)
  f1.pe <- f1.pe/length(errs)
  sens.rpbat <- sens.rpbat/length(errs)
  spec.rpbat <- spec.rpbat/length(errs)
  f1.rpbat <- f1.rpbat/length(errs)

  ans$pe$avg$sensitivity <- sens.pe
  ans$pe$avg$specificity <- spec.pe
  ans$pe$avg$f1 <- f1.pe
  ans$rpbat$avg$sensitivity <- sens.rpbat
  ans$rpbat$avg$specificity <- spec.rpbat
  ans$rpbat$avg$f1 <- f1.rpbat

  return(ans)
}

################# MAIN ################
do.sim <- T
do.tbl <- T
do.times <- T
do.mems <- T
do.collision <- T
plot <- T

sim.tbl <- NULL
tbl <- NULL
tbl.primary <- NULL
times <- NULL
mems <- NULL
primary <- unlist(datasets$primary) == "yes"

if (do.sim) {
  tryCatch({
    write.log("reading simulation output...")
    read.lengths <- seq(from = 50, to = 150, by = 10)
    sim.tbl <- make.sim.tbl(read.lengths)
  }, error = function(e) {
    write.log("failure reading simulation")
  })
} else {
  write.log("SIMULATION OUTPUT SKIPPED")
  do.sim <- F
}

if (do.tbl) {
  tryCatch({
    write.log("reading master table: methpipe/samtools outputs...")
    tbl <- make.table(datasets)
    tbl.primary <- tbl[primary,]
  }, error = function(e) {
    print(e)
    write.log("failure reading master table")
    do.tbl <- F
  })
} else {
  write.log("MASTER TABLE SKIPPED")
}

if (do.times) {
  tryCatch({
    write.log("reading snakemake time outputs...")
    times <- make.times(tbl)
  }, error = function(e) {
    write.log("failure reading time outputs")
    do.times <- F
  })
} else {
  write.log("TIME DATA SKIPPED")
}

if (do.mems) {
  tryCatch({
    write.log("reading snakemake memory outputs...")
    mems <- make.mems(tbl)
  }, error = function(e) {
    write.log("failure reading memory outputs")
    do.mems <- F
  })
} else {
  write.log("MEMORY DATA SKIPPED")
}

if (plot) {
  if (do.tbl) {

    write.log("SUPP TABLE 1: DATASETS")
    wt(datasets[datasets$primary == "yes",], "results/tables/supp_table_1_datasets.tsv")

    # Supp table 1 is accuracy
    conc.pairs <- round(100*get.data(tbl, "map.total"), 3)
    sam.err.rate <- round(100*get.data(tbl, "sam.err.rate", F), 3)
    supp.tab.1 <- cbind(conc.pairs, sam.err.rate)

    # Supp table 2 is bsrate and coverage
    bs.conv <- round(100*get.data(tbl, "bs.conv.tot"), 3)
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
    supp.tab.2 <- cbind(conc.pairs, sam.err.rate, bs.conv, c.cov, c.depth, c.meth.mean, c.meth.weigh, c.meth.frac,
             cpg.cov, cpg.depth, cpg.meth.mean, cpg.meth.weigh, cpg.meth.frac)

    write.log("SUPP TABLE 2: ACCURACY...")
    wt(supp.tab.2, "results/tables/supp_table_2_downstream.tsv")
  }
  
  
  if (do.times & do.mems) {
    # Supp table 3 are the times
    a <- times$times[,mappers]
    b <- mems[,mappers]
    colnames(a) <- paste0("time_", colnames(a))
    colnames(b) <- paste0("mem_", colnames(b))

    write.log("SUPP TABLE 3: RESOURCES...")
    wt(a, "results/tables/supp_table_3_1_times.tsv")
    wt(b, "results/tables/supp_table_3_2_mems.tsv")
  }

  ########### FIGURES ###########
  if (do.sim) {
    write.log("FIGURE S1: SIMULATION RESULTS PE...")
    make.sim.comparison.figure(sim.tbl$pe,
      "results/figures/supp_1_simulations_pe.pdf"
    )
    write.log("FIGURE S2: SIMULATION RESULTS RPBAT...")
    make.sim.comparison.figure(sim.tbl$rpbat,
      "results/figures/supp_2_simulations_rpbat.pdf"
    )

  }
  if (do.collision) {
    write.log("FIGURE 1B: COLLISION RATES...")
    make.hit.ratio.figure("results/figures/fig1b_collision.pdf")
  }

  if (do.tbl) {
    write.log("FIGURE 2: ACCURACY...")
    make.accuracy.figure(sim.tbl, tbl.primary)
  }

  if (do.times) {
    write.log("FIGURE 3: RESOURCES...")
    times.primary <- times
    mems.primary <- mems
    times.primary$times <- times$times[primary,]
    times.primary$meta <- times$meta[primary,]

    primary.single <- primary[datasets$protocol == "wgbs_single"]
    mems.primary <- mems[primary.single,]

    # for the main figure we filter only the primary tests
    datasets <- datasets[as.character(datasets$primary_test) == "yes",]

    # resources
    make.resources.figure(times.primary, mems.primary)
  }
}
if (file.exists("Rplots.pdf"))
  file.remove("Rplots.pdf")
write.log("DONE!")
