### RELEVANT LIBRARIES ### 
#library(readxl)
#library(outbreaks)
#library(incidence)
#library(epicontacts)
library(ggplot2)
library(magrittr)
library(binom)
library(ape)
library(outbreaker2)
library(igraph)
library(lubridate)

# DNA - change file and file path as desired
dna <- read.FASTA("/Users/gmoreno/Downloads/relevant_for_Gage/all_peacock_remapped_taxids.high-qual.msa.UTRs-masked.gappy-removed.fasta")
names(dna) <- gsub("/", "-", names(dna))
dna

# Metadata
linelist <- read.table(file = "/Users/gmoreno/Downloads/relevant_for_Gage/specimen_dates.tsv", sep = '\t', header = TRUE)
linelist$collection_or_symptom_date <- as.Date(linelist$collection_or_symptom_date)
linelist$barcode <- gsub("/", "-", linelist$barcode)
linelist <- subset(linelist, barcode %in% names(dna))
#linelist <- subset(linelist, percent_genome_covered_95_percent == "Yes")

### Contact Flat Pairs ----

# Contact tracing data
contacts <- read.csv("/Users/gmoreno/Downloads/relevant_for_Gage/2021_contacts_flat_pairs_with_peacock-peacock_only.tsv", sep='\t',row.names=NULL)
contacts$from <- gsub("/", "-", contacts$from)
contacts$to <- gsub("/", "-", contacts$to)

# get test dates
dates <- linelist$collection_or_symptom_date
names(dates) <- linelist$barcode
dna <- dna[which(names(dna) %in% names(dates))]


# generation interval
mean_gen <- 5.2
sd_gen <- 1.5
shape <- (mean_gen^2)/(sd_gen^2)
rate <- (mean_gen)/(sd_gen^2)

## incubtation/colonization time
mean_inc <- 6
sd_inc <- 1.5
shape_inc <- (mean_inc^2)/(sd_inc^2)
rate_inc <-(mean_inc)/(sd_inc^2)

## setup outbreaker2
w <- dgamma(1:21, shape=shape, rate=rate) # gamma distribution for serial transmission interval
i <- dgamma(1:21, shape=shape_inc, rate=rate_inc)
n_iter <- 1e4 # number of iterations to run outbreaker2
burn <- n_iter*.1 # burn-in, i.e. discard first 10%
set.seed(0) # set random state

# Mutation rate
mu <- 4.5e-6
### RUN OUTBREAKER2 ### 

# set up configuration for outbreaker2
config <- create_config(#find_import=FALSE,
  init_tree = "random", #tree_pairs,
  #init_alpha = tree_pairs,
  n_iter=n_iter,
  #init_mu = mu,
  pb=TRUE
  #move_mu=FALSE
)
data <- outbreaker_data(
  dna=dna,
  dates=dates,
  w_dens=w,
  f_dens=i,
  ctd=contacts
)

res <- outbreaker(data, config)

tree <- summary(res)$tree

# Convert names to case IDs
tree$from <- names(dates)[tree$from]
tree$to <- names(dates)[tree$to]

# note: this function partially taken from outbreaker source code
get_transmissions <- function(x, burnin = burn, min_support = 0.1, labels = NULL, ...){
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  alpha <- as.matrix(x[,grep("alpha", names(x))])
  colnames(alpha) <- seq_len(ncol(alpha))
  from <- as.vector(alpha)
  to <- as.vector(col(alpha))
  from[is.na(from)] <- 0
  out_dat <- data.frame(xyTable(from,to))
  names(out_dat) <- c("from", "to", "frequency")
  ## Calculate proportion among ancestries
  get_prop <- function(i) {
    ind <- which(out_dat$to == out_dat$to[i])
    out_dat[[3]][i]/sum(out_dat[[3]][ind])
  }
  
  #out_dat[3] <- vapply(seq_along(out_dat[[3]]), get_prop, 1)
  out_dat[,3] <- out_dat[,3] / nrow(x)
  
  out_dat
}

transmissions <- get_transmissions(res)

rename <- function(x){
  if(x==0){
    NA
  }else{
    names(dates)[x]
  }
}

transmissions$from <- sapply(transmissions$from, rename)
transmissions$to <- sapply(transmissions$to, rename)

write.csv(transmissions, "/Users/gmoreno/Downloads/relevant_for_Gage/OB2_2021_contacts_flat_pairs_1e4iterations.csv")

### Wifi Spring 2d ----
# Contact tracing data
contacts <- read.csv("/Users/gmoreno/Downloads/relevant_for_Gage/wifi_spring_2d_peacock-peacock_only.tsv", sep='\t',row.names=NULL)
contacts$from <- gsub("/", "-", contacts$from)
contacts$to <- gsub("/", "-", contacts$to)

# get test dates
dates <- linelist$collection_or_symptom_date
names(dates) <- linelist$barcode
dna <- dna[which(names(dna) %in% names(dates))]


# generation interval
mean_gen <- 5.2
sd_gen <- 1.5
shape <- (mean_gen^2)/(sd_gen^2)
rate <- (mean_gen)/(sd_gen^2)

## incubtation/colonization time
mean_inc <- 6
sd_inc <- 1.5
shape_inc <- (mean_inc^2)/(sd_inc^2)
rate_inc <-(mean_inc)/(sd_inc^2)

## setup outbreaker2
w <- dgamma(1:21, shape=shape, rate=rate) # gamma distribution for serial transmission interval
i <- dgamma(1:21, shape=shape_inc, rate=rate_inc)
n_iter <- 1e4 # number of iterations to run outbreaker2
burn <- n_iter*.1 # burn-in, i.e. discard first 10%
set.seed(0) # set random state

# Mutation rate
mu <- 4.5e-6
### RUN OUTBREAKER2 ### 

# set up configuration for outbreaker2
config <- create_config(#find_import=FALSE,
  init_tree = "random", #tree_pairs,
  #init_alpha = tree_pairs,
  n_iter=n_iter,
  #init_mu = mu,
  pb=TRUE
  #move_mu=FALSE
)
data <- outbreaker_data(
  dna=dna,
  dates=dates,
  w_dens=w,
  f_dens=i,
  ctd=contacts
)

res <- outbreaker(data, config)

tree <- summary(res)$tree

# Convert names to case IDs
tree$from <- names(dates)[tree$from]
tree$to <- names(dates)[tree$to]

# note: this function partially taken from outbreaker source code
get_transmissions <- function(x, burnin = burn, min_support = 0.1, labels = NULL, ...){
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  alpha <- as.matrix(x[,grep("alpha", names(x))])
  colnames(alpha) <- seq_len(ncol(alpha))
  from <- as.vector(alpha)
  to <- as.vector(col(alpha))
  from[is.na(from)] <- 0
  out_dat <- data.frame(xyTable(from,to))
  names(out_dat) <- c("from", "to", "frequency")
  ## Calculate proportion among ancestries
  get_prop <- function(i) {
    ind <- which(out_dat$to == out_dat$to[i])
    out_dat[[3]][i]/sum(out_dat[[3]][ind])
  }
  
  #out_dat[3] <- vapply(seq_along(out_dat[[3]]), get_prop, 1)
  out_dat[,3] <- out_dat[,3] / nrow(x)
  
  out_dat
}

transmissions <- get_transmissions(res)

rename <- function(x){
  if(x==0){
    NA
  }else{
    names(dates)[x]
  }
}

transmissions$from <- sapply(transmissions$from, rename)
transmissions$to <- sapply(transmissions$to, rename)

write.csv(transmissions, "/Users/gmoreno/Downloads/relevant_for_Gage/OB2_wifi_spring_2d_1e4iterations.csv")

### Wifi Spring 10d ----
# Contact tracing data
contacts <- read.csv("/Users/gmoreno/Downloads/relevant_for_Gage/wifi_spring_10d_peacock-peacock_only.tsv", sep='\t',row.names=NULL)
contacts$from <- gsub("/", "-", contacts$from)
contacts$to <- gsub("/", "-", contacts$to)

# get test dates
dates <- linelist$collection_or_symptom_date
names(dates) <- linelist$barcode
dna <- dna[which(names(dna) %in% names(dates))]


# generation interval
mean_gen <- 5.2
sd_gen <- 1.5
shape <- (mean_gen^2)/(sd_gen^2)
rate <- (mean_gen)/(sd_gen^2)

## incubtation/colonization time
mean_inc <- 6
sd_inc <- 1.5
shape_inc <- (mean_inc^2)/(sd_inc^2)
rate_inc <-(mean_inc)/(sd_inc^2)

## setup outbreaker2
w <- dgamma(1:21, shape=shape, rate=rate) # gamma distribution for serial transmission interval
i <- dgamma(1:21, shape=shape_inc, rate=rate_inc)
n_iter <- 1e4 # number of iterations to run outbreaker2
burn <- n_iter*.1 # burn-in, i.e. discard first 10%
set.seed(0) # set random state

# Mutation rate
mu <- 4.5e-6
### RUN OUTBREAKER2 ### 

# set up configuration for outbreaker2
config <- create_config(#find_import=FALSE,
  init_tree = "random", #tree_pairs,
  #init_alpha = tree_pairs,
  n_iter=n_iter,
  #init_mu = mu,
  pb=TRUE
  #move_mu=FALSE
)
data <- outbreaker_data(
  dna=dna,
  dates=dates,
  w_dens=w,
  f_dens=i,
  ctd=contacts
)

res <- outbreaker(data, config)

tree <- summary(res)$tree

# Convert names to case IDs
tree$from <- names(dates)[tree$from]
tree$to <- names(dates)[tree$to]

# note: this function partially taken from outbreaker source code
get_transmissions <- function(x, burnin = burn, min_support = 0.1, labels = NULL, ...){
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  alpha <- as.matrix(x[,grep("alpha", names(x))])
  colnames(alpha) <- seq_len(ncol(alpha))
  from <- as.vector(alpha)
  to <- as.vector(col(alpha))
  from[is.na(from)] <- 0
  out_dat <- data.frame(xyTable(from,to))
  names(out_dat) <- c("from", "to", "frequency")
  ## Calculate proportion among ancestries
  get_prop <- function(i) {
    ind <- which(out_dat$to == out_dat$to[i])
    out_dat[[3]][i]/sum(out_dat[[3]][ind])
  }
  
  #out_dat[3] <- vapply(seq_along(out_dat[[3]]), get_prop, 1)
  out_dat[,3] <- out_dat[,3] / nrow(x)
  
  out_dat
}

transmissions <- get_transmissions(res)

rename <- function(x){
  if(x==0){
    NA
  }else{
    names(dates)[x]
  }
}

transmissions$from <- sapply(transmissions$from, rename)
transmissions$to <- sapply(transmissions$to, rename)

write.csv(transmissions, "/Users/gmoreno/Downloads/relevant_for_Gage/OB2_wifi_spring_10d_1e4iterations.csv")
