#!/usr/bin/env Rscript

library(ape)
#devtools::install_github('dosreislab/simclock')
library(simclock)


parse_args <- function() {
  cargs <- commandArgs(trailingOnly=TRUE)
  if ((length(cargs)==1) && (cargs[1] %in% c('-h', '--help'))) {
    cat("Generate INDELible control files with trees simulated under relaxed clocks.

OPTIONS:
    -h print help
    -n number of simulations [1000].
    -c clock model (clk/iln/gbm). Default all.
    -r mean rate. Default 5e-4, 15e-4, 40e-4, 80e-4.
    -L sequence length. Default 300, 3000, 30000.

Clock models:
    clk: Clock-like (constant rate)
    iln: Independent log-normal
    gbm: Geometric Brownian Motion
")
    quit(save='no', status=0)
  }
  argvalues <- list('-n'=1000,
                    '-c'=c('clk', 'iln', 'gbm'),
                    '-r'=c('5e-4', '15e-4', '40e-4', '80e-4'),
                    '-L'=c(300, 3000, 30000))
  if (!length(cargs))
    return(argvalues)

  key <- NULL
  for (k in seq_along(cargs)) {
    if (cargs[k] %in% c('-n', '-c', '-r', '-L')) {
      if (!is.null(key))
        stop(paste('Missing argument for option', key))
      key <- cargs[k]
    } else if (is.null(key)) {
      stop(paste('Extra argument', cargs[k]))
    } else if (key %in% c('-n', '-L')) {
      argvalues[key] <- as.integer(cargs[k])
      key <- NULL
    } else {
      argvalues[key] <- cargs[k]
      key <- NULL
    }
  }
  if (!is.null(key))
    stop(paste('Missing argument for option', key))

  return(argvalues)
}

argvalues <- parse_args()

# Tree with branch lengths in million years.
phyltree <- read.tree('PhylTree-Primates.TimeTree201901.Ensembl93-like.goodQual.outgroupsMic-Pco.seqnames.nwk')

# Load the tree without internal labels or set node.labels to NULL
# Then the written newick will be in the expected format for INDELible.
if ( !is.null(phyltree$node.label) )
  phyltree$node.label <- NULL
if ( !is.null(phyltree$root.edge) )
  phyltree$root.edge <- NULL

set.seed(55555)

nsim         <- argvalues[['-n']]
clock_models <- argvalues[['-c']]
seq_lengths  <- argvalues[['-L']]
rate_txt     <- argvalues[['-r']]
rates <- as.double(rate_txt)
#diffu <- c(5e-4, 15e-4, 40e-4, 80e-4)  #FIXME: random guess

for (model in clock_models) {
for (k in seq_along(rates)) {
  rk <- rates[k]
  #dk <- diffu[k]

  # Insertion-deletion rates
  irate <- (1/14)/3 * 0.00177/rk
  drate <- (2/14)/3 * 0.00177/rk
  #TODO: write an INDELIBLE control file
  for (len in seq_lengths) {

    #NOTE: it could be acceptable to use the same trees for all lengths.
    # Parameter s2 : sigma^2 (subst/My)^2
    simulated <- lapply(1:nsim,
      function(i) {relaxed.tree(phyltree, model='iln', r=rk, s2=50*rk)}
      )

    dir.create(sprintf('clock_%s/L%dn_t%s', model, len, rate_txt[k]), recursive=TRUE)
    filename <- sprintf('clock_%s/L%dn_t%s/control.txt', model, len, rate_txt[k])
    write(paste('#', filename), file=stderr())

    #out <- stdout()
    out <- file(filename, 'wt')

    indices <- seq(0, nsim-1)

    write(file=out, c(
      sprintf("//  INDELible V1.03 control file
//  length = %d nucleotides; t rate (subst/nucl) = %s;
//  clock model '%s' generated with R simclock (by dosreislab)

[TYPE] CODON 1

[SETTINGS]
    [randomseed] 9342
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
  [statefreq]", len, rk, model, irate, drate),
  readLines('codonsfreq_m0ns5.txt'),
  "",
  sprintf("[TREE] tree%d %s",
          indices,
          #sub('Primates:0;$', ';', write.tree(do.call(c, simulated)))
          write.tree(do.call(c, simulated))  # uses c.phylo which creates a multiPhylo object.
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

warnings()

# Alternative parameters:

# For omega constant per site, use:
# [submodel]  4  0.175  //  Substitution model is M0 with kappa=4, omega=0.175
