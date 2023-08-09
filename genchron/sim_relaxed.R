#!/usr/bin/env Rscript

library(ape)
#devtools::install_github('dosreislab/simclock')
library(simclock)  # Requires at least version 0.1.4 to handle polytomic trees in .sim.gbm
library(stringr)

# Do not use.
poisson.clock <- function(tree, r, seqlength=3000) {
  tt <- tree
  nb <- length(tt$edge.length)
  tt$edge.length <- tt$edge.length * rpois(nb, lambda=seqlength*r)/seqlength
}

# Compute the path lengths of all pairs edge midpoints
#dist.edge.midpoints <- function(tree) {
#
#}

# Indelible cannot parse branch lengths in scientific notation... and APE does
# not give control over the branch length format.
reformat.lengths <- function(newicks, prec=10) {
  ntrees <- length(newicks)
  fmt <- paste0('%.', prec, 'f')
  pattern <- '[0-9]+(\\.[0-9]+)?[Ee][+-]?[0-9]+'  # scientific notation. Ignore case.

  matches <- str_extract_all(unlist(newicks), pattern)
  numbers <- lapply(matches, function(m) {gsub('\\.?0+$', '',
                                               sprintf(fmt, as.double(m)))})
  unmatched <- strsplit(unlist(newicks), pattern)

  sapply(1:ntrees, function(i) {paste0(unmatched[[i]], numbers[[i]], collapse='')})

}

# Unused. For manual verification only. Recompute the resulting rate deviation
effective.rate.dev <- function(simulated, phyltree) {
  rates <- sapply(simulated, function(tree) {tree$edge.length/phyltree$edge.length})
  # trees are in columns, edges in rows.
  meanrate <- mean(rates)
  meansd_per_tree <- mean(apply(rates, 2, sd))
  
  list(meanrate=meanrate, meansd_per_tree=meansd_per_tree)
}

average.midpoint.depth <- function(tree) {
  depths <- node.depth.edgelength(tree)
  edgetips <- tree$edge[,2]
  midpoint.depths <- depths[edgetips] - (tree$edge.lengths)
  mean(midpoint.depths)
}

parse_args <- function() {
  cargs <- commandArgs(trailingOnly=TRUE)
  if ((length(cargs)==1) && (cargs[1] %in% c('-h', '--help'))) {
    cat("Generate INDELible control files with trees simulated under relaxed clocks.

OPTIONS:
    -h print help
    -n number of replicates for each parameter set [500].
    -R random seeds for each set of parameters [random one for each, written out]
    -d output directory format: ['clock_%1$s/L%2$dn_t%3$s_s%4$s%5$s']. Use 1$,2$,3$,4$,$5 indices for clock, length, rate, diffu.ctrl, s2.
    -p print only but do not run (debugging).

    -c clock model (iln/gbm). Default all.
    -L sequence length. Default 300, 3000, 30000.
    -r mean rate. Default 5e-4, 15e-4, 40e-4, 80e-4.
    
    -s sigma^2 (diffusion parameter). Default ln(2) means rate std == rate mean.
    -a automatic diffusion: interpret -s as the *real* rate std and compute sigma^2 = log(1 + r^2/s^2)

Clock models:
    iln: Independent log-normal (simclock)
    gbm: Geometric Brownian Motion (simclock)
    
    NOTE: simclock implements 'clk', which just multiplies lengths by the rate (no variation)
")
    #pois: Poisson variation (constant mean rate).

    quit(save='no', status=0)
  }
  argvalues <- list('-n'=500,
                    '-d'='clock_%1$s/L%2$dn_t%3$s_%4$s%5$s',
                    '-R'=NULL,
                    '-c'=c('iln', 'gbm'),
                    '-L'=c(300, 3000, 30000),
                    '-r'=c('5e-4', '15e-4', '40e-4', '80e-4'),
                    '-s'=c('ln(2)'),
                    '-p'=FALSE,
                    '-a'=FALSE)
  if (!length(cargs))
    return(argvalues)

  keynames <- names(argvalues)
  extendable <- c('-R', '-c', '-L', '-r', '-s')
  key <- NULL
  for (k in seq_along(cargs)) {
    if (cargs[k] %in% keynames) {
      if ( !is.null(key) && length(argvalues[[key]])==0 )
        stop(paste('Missing argument for option', key))
      key <- cargs[k]
      if (key %in% c('-p', '-a')) {
        argvalues[[cargs[k]]] <- TRUE
        key <- NULL
      } else {
        argvalues[[key]] <- c()
      }
      next
    }
    if (is.null(key)) {
      stop(paste('Extra argument', cargs[k]))
    }
    argvalues[[key]] <- c(argvalues[[key]], cargs[k])
    if ( ! (key %in% extendable) ) {
      key <- NULL
    }
  }
  if ( !is.null(key) && length(argvalues[[key]])==0 )
    stop(paste('Missing argument for option', key))

  return(argvalues)
}

argvalues <- parse_args()

debugging <- argvalues[['-p']]

# Tree with branch lengths in million years.
#FIXME: use any species tree. The following is available at DOI:10.5281/zenodo.8235767.
phyltree <- read.tree('PhylTree-Primates.TimeTree201901.Ensembl93-like.goodQual.outgroupsMic-Pco.seqnames.nwk')

# Load the tree without internal labels or set node.labels to NULL
# Then the written newick will be in the expected format for INDELible.
if ( !is.null(phyltree$node.label) )
  phyltree$node.label <- NULL
if ( !is.null(phyltree$root.edge) )
  phyltree$root.edge <- NULL

for (key in c('-n', '-L', '-R')) {
  argvalues[[key]] <- as.integer(argvalues[[key]])
}
dirfmt       <- argvalues[['-d']]
nrep         <- argvalues[['-n']]
clock_models <- argvalues[['-c']]
seq_lengths  <- argvalues[['-L']]
rate_txt     <- argvalues[['-r']]
rates <- as.double(rate_txt)
diffu_txt <- argvalues[['-s']]
diffu_txt[diffu_txt=='ln(2)'] <- "ln2"
diffu <- as.double(ifelse(diffu_txt=='ln2', 0.693147180559945, diffu_txt))

if (argvalues[['-a']]) {
  if ( 'ln2' %in% diffu_txt )
    warning('Unexpected combination -a and -s ln(2)')
  diffu <- log(1 + rates^2 / diffu^2)
  # NOTE that this changes the iteration: values of rates,diffu are iterated jointly (instead of product)
  auto <- 'auto'
  diffu.ctrl <- 'std'
} else {
  auto <- ''
  diffu.ctrl <- 'diffu'
}

nparams <- length(clock_models) * length(seq_lengths) * length(rates)
if ( !argvalues[['-a']] )
  nparams = nparams * length(diffu)


if ( length(argvalues[['-R']]) == 0 ) {
  sim_seeds <- sample.int(65535, nparams)
} else {
  sim_seeds <- rep_len(argvalues[['-R']], nparams)
}
cat('Will simulate', nrep, 'replicates for', nparams, 'sets of parameters.', length(argvalues[['-R']]), 'random seeds.\n', file=stderr())

sim <- 1
for (model in clock_models) {
for (k in seq_along(rates)) {
  rk <- rates[k]

  # Insertion-deletion rates
  #FIXME: 0.00177 substitutions/codon/My was obtained on my 21 Primates tree from Ensembl 93.
  # By dividing by 'rk' (mean subst rate), we ensure that all simul will generate alignments that
  # have the same quantities of gaps.
  irate <- (1/14)/3 * 0.00177/rk
  drate <- (2/14)/3 * 0.00177/rk

  seq_diffu <- if( argvalues[['-a']] ) { k } else { seq_along(diffu) }

  for (j in seq_diffu) {
  s2_j <- diffu[j]
  for (len in seq_lengths) {

    set.seed(sim_seeds[sim])

    #NOTE: it could be acceptable to use the same trees for all lengths.
    # Parameter s2 : sigma^2 (subst/My)^2
    simulated <- lapply(1:nrep,
      function(i) {relaxed.tree(phyltree, model=model, r=rk, s2=s2_j)}
      )
    indel_seed <- sample.int(65535, 1)

    outdir <- sprintf(dirfmt, model, len, rate_txt[k], diffu.ctrl, diffu_txt[j])
    filename <- paste0(outdir, '/control.txt')
    write(paste('#', sim, filename, 'seed:', sim_seeds[sim]), file=stderr())
    sim <- sim + 1
    if (debugging) next
    dir.create(outdir, recursive=TRUE)

    #out <- stdout()
    out <- file(filename, 'wt')

    indices <- seq(0, nrep-1)

    write(file=out, c(
      sprintf("//  INDELible V1.03 control file
//  length = %d nucleotides; t rate (subst/nucl): mean=%s %s=%s;
//  clock model '%s' generated with R simclock (by dosreislab) with seed %d

[TYPE] CODON 1

[SETTINGS]
    [randomseed] %d
    [fastaextension] fa

[MODEL]    model0            //  Evolutionary models are defined in [MODEL] blocks.
  //  Substitution model is discrete gamma with kappa=4 and 4 omega categories
  //  Mean omega is 0.256
  // p:   0.25000  0.25000  0.25000  0.24700 0.003
  // w:   0.00000  0.00298  0.07966  0.88577 5.53
  [submodel]    4
                0.25000  0.25000  0.25000  0.24700
                0.00000  0.00298  0.07966  0.88577 5.53
  [insertrate]  %.2g
  [deleterate]  %.2g
  [indelmodel]  LAV 2 300
  [statefreq]", len, rate_txt[k], diffu.ctrl, diffu_txt[j], model, sim_seeds[sim], indel_seed, irate, drate),
  readLines('codonsfreq_m0ns5.txt'),
  "",
  sprintf("[TREE] tree%d %s",
          indices,
          #sub('Primates:0;$', ';', write.tree(do.call(c, simulated)))
          reformat.lengths(
            write.tree(do.call(c, simulated))  # uses c.phylo which creates a multiPhylo object.
          )
         ),
  "\n// [PARTITIONS] blocks say which models go with which trees and define",
  "// the length of the  sequence generated at the root.",
  sprintf("[PARTITIONS] partition%d
  [tree%d model0 %4d]", indices, indices, len/3),
  "\n// Evolve 1 replicate of each partition",
  "[EVOLVE]",
  sprintf("  partition%d 1 out%03d", indices, indices)
  ))
  close(out)
  }
  }
}
}

warnings()

# Alternative parameters:

# For omega constant per site, use:
# [submodel]  4  0.175  //  Substitution model is M0 with kappa=4, omega=0.175
