#!/usr/bin/env Rscript

usage="
date_dup.R <datasetname> <calibfile> <treefile> <ncores>
"


library(ape)
library(parallel)
source('~/scripts/dendro/chronos-function.R')  # Own fixed version


# Helper function for interactive investigation:
display_edges <- function(tree) {
  data.frame(tree$edge,
        node.label=c(tree$tip.label, tree$node.label)[tree$edge[,2]],
        edge.length=tree$edge.length)
}
display_calibrated_edges <- function(tree, calib) {
  dedges <- display_edges(tree)
  # NOTE: this won't show the root node age.
  cbind(dedges,
        calib=calib[match(dedges$node.label, row.names(calib)),"age.min"]
        )
}

load_calibration <- function(agefile) {
  #ages <- read.delim(agefile, row.names=1)
  #calib <- ages[ages$type == "spe", c("age_dist", "type")]
  #names(calib) <- c("name", "age.min", "type")
  #calib$age.max <- calib$age.min
  # cbind
  # calib$soft.bounds
  # calib$node
  #return(calib)

  # Load without row.names, then drop duplicate row names.
  ages <- read.delim(agefile, row.names=NULL, header=TRUE)
  dup_rownames <- duplicated(ages[,'name'])
  if( any(dup_rownames) ) {
    cat('WARNING: drop', sum(dup_rownames), 'duplicate row.names from ages:\n',
        paste(head(ages[dup_rownames,'name'], 20), collapse='\n '), '...\n', file=stderr())
    #TODO: remove, based on the 'is_outgroup' flag.
    #if('is_outgroup' %in% colnames(ages)) {
    #  ages <- ages[(!dup_rownames & !ages$is_outgroup),]
    #} else {
      ages <- ages[!dup_rownames,]
  }
  rownames(ages) <- ages[,'name']
  ages <- ages[,-1]  # Assuming 'name' is the first column
  return(ages)
}

subset_calibration <- function(ages, tree, age_col="age_dist", leaf_age=0) {
  # return a "calibration" dataframe formatted for ape::chronos.
  
  ### Does not retrieve leaf ages, because chronos doesn't accept them.
  ntips = Ntip(tree)
  subcalib_names <- tree$node.label
  subcalib_nodes <- ntips + seq_along(subcalib_names)
  subcalib_agemin <- ages[match(subcalib_names, rownames(ages)), age_col] - leaf_age
  
  calibrated <- ages[match(subcalib_names, rownames(ages)), "calibrated"] == 1
  # handle NA type nodes (not found in ages: duplication as root)
  # Or it is the outgroup # Better to export the subtrees with generate_dNdStable.py

  #max_age <- max(subcalib_agemin, na.rm=TRUE)
  #subcalib_agemin[is.na(calibrated)] <- max_age
  #calibrated[is.na(calibrated)] <- TRUE  # Need those nodes before the first calibration for MPL stats.
  calibrated[is.na(calibrated)] <- FALSE

  subcalib <- data.frame(node=subcalib_nodes,
                         age.min=subcalib_agemin,
                         age.max=subcalib_agemin,
                         soft.bounds=FALSE,
                         row.names=subcalib_names)[calibrated,]
}


# the MPL algo works only with purely dichotomic trees. So I need to (randomly)
# convert multifurcations to bifurcations, but then convert them back with the
# correct length
delete_edges <- function(tree, removed_edges) {
  #tree <- reorder(tree, "cladewise")
  # attr(tree, 'order') == 'cladewise'
  
  e <- tree$edge
  el <- tree$edge.length

  # iterate in reverse order (normally from leaves to root)
  removed_edges <- sort(removed_edges, decreasing=TRUE)
  removed_nodes <- e[removed_edges,2]
  removed_int_nodes <- removed_nodes - Ntip(tree)

  tree$node.label <- tree$node.label[-removed_int_nodes]

  for (i in removed_edges) {
    # Find the edges *starting* with this node.
    removed_node <- e[i,2]
    children_edges <- which(e[,1] == removed_node)

    # Update the branch start to the base of the deleted edge
    e[children_edges, 1] <- e[i,1]

    # Add the length of the current edge to the descendant edges
    el[children_edges] <- el[children_edges] + el[i]
  }

  tree$edge.length <- el[-removed_edges]

  # Update the node numbering.
  for (i in removed_nodes) {
    # All edges should be numbered consecutively, to be assigned to an edge.length
    e[e>i] <- as.integer(e[e>i] - 1)
  }

  tree$edge <- e[-removed_edges,]
  tree$Nnode <- as.integer(tree$Nnode - length(removed_edges))

  a <- attributes(tree)

  if( "stderr" %in% names(a) )
    attr(tree, 'stderr') <- attr(tree, 'stderr')[-removed_int_nodes]
  if( "Pval" %in% names(a) )
    attr(tree, 'Pval') <- attr(tree, 'Pval')[-removed_int_nodes]

  #checkValidPhylo(tree)
  return(tree)
}


mark_extra_edges <- function(tree, condition, value=-10^6) {
  # Old condition: function(nodenb) nodenb>maxnb
  e <- tree$edge
  el <- tree$edge.length

  removed_edges <- which(condition(e[,2]))
  
  #print(length(removed_edges))
  for(r in removed_edges) {
    # extend the lengths of edges **starting** with that node.
    descendants_of_removed <- which(e[,1] == e[r,2])
    tree$edge.length[descendants_of_removed] <- el[descendants_of_removed] + el[r]
    tree$edge.length[r] <- value
  }
  return(tree)
}

restore_multi <- function(dicho_tree, orig_pattern='\\.orig$') {
  ntips <- Ntip(dicho_tree)
  all_labels <- c(dicho_tree$tip.label, dicho_tree$node.label)

  nnb <- dicho_tree$edge[,2]

  return(delete_edges(dicho_tree,
                      which((nnb > ntips)
                            & !grepl(orig_pattern, all_labels[nnb]))))
}

# Available dating methods:
#  - chronos: penalized likelihood
#  - chronos(phy, lambda=Inf): NPRS
#  - chronosMPL: Mean Path Length, Britton et al

## First, the arguments controling the optimisation
mychronos_control <- chronos.control(epsilon=1e-15)#tol=1e8, iter.max=1e5, eval.max=1e5, dual.iter.max=40)
mychronos_control_C <- mychronos_control
mychronos_control_C[['nb.rate.cat']] <- 1

chronos_PL1 <- function(tree, calibration, lambda=1, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
chronos_PL100 <- function(tree, calibration, lambda=100, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}
chronos_PL10000 <- function(tree, calibration, lambda=10000, ...) {
  return(chronos(tree, lambda=lambda, calibration=calibration, ...))
}

chronos_NPRS <- function(tree, calibration, ...) {
  return(chronos(tree, lambda=10^16, calibration=calibration, ...))
}

chronos_MPL <- function(subtree, calibration, ...) {
  # TODO: split the input tree into subtrees with calibrated leaves
  if ( !is.binary(subtree) ) subtree <- multi2di(subtree)
  return(chronoMPL(subtree, ...))
}

date_PL0 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=0, calibration=subcalib, quiet=TRUE,
                 control=mychronos_control))
}
date_PL1 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=1, calibration=subcalib, quiet=TRUE,
                 control=mychronos_control))
}
date_PL100 <- function(subtree, subcalib, quiet=TRUE) {
  return(chronos(subtree, lambda=100, calibration=subcalib, quiet=quiet,
                 control=mychronos_control))
}
date_PL10000 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=10000, calibration=subcalib, quiet=TRUE,
                 control=mychronos_control))
}
date_NPRS <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=1e16, calibration=subcalib, quiet=TRUE,
                 control=mychronos_control))
}
date_R1 <- function(subtree, subcalib) {
  return(chronos(subtree, model='relaxed', calibration=subcalib, quiet=TRUE,
                 control=mychronos_control))
}
date_R100 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=100, model='relaxed', calibration=subcalib,
                 quiet=TRUE, control=mychronos_control))
}
date_R10000 <- function(subtree, subcalib) {
  return(chronos(subtree, lambda=10000, model='relaxed', calibration=subcalib,
                 quiet=TRUE, control=mychronos_control))
}
date_C <- function(subtree, subcalib) {
  return(chronos(subtree, model='discrete', calibration=subcalib,
                 control=mychronos_control_C, quiet=TRUE))
}

date_MPL <- function(subtree, subcalib) {
  is_multi <- !is.binary(subtree)
  if(is_multi) {
    # Mark original nodes
    subtree$node.label <- paste0(subtree$node.label, '.orig')
    subtree <- multi2di(subtree)
  }
  chronogram <- chronoMPL(subtree, se=TRUE, test=TRUE)
  #neg_branchlen <- chronogram$edge.length < 0
  if( is_multi ) {
    chronogram <- restore_multi(chronogram, '\\.orig$')
    chronogram$node.label <- sub('\\.orig$', '', chronogram$node.label)
  }
  # And rescale edge lengths
  MPL_depth <- max(node.depth.edgelength(chronogram))
  real_depth <- subcalib[chronogram$node.label[1], "age.min"]
  if( is.na(real_depth) )
    stop(paste("While rescaling MPL: calibration gives no age for root node", chronogram$node.label[1]))
  # If multiple calibration, do that separately for each which.edge(rowname of subcalib)
  chronogram$edge.length <- chronogram$edge.length * real_depth/MPL_depth
  attr(chronogram, "stderr") <- attr(chronogram, "stderr") * real_depth/MPL_depth
  return(chronogram)
}


KL.div <- function(freqs1, freqs2) {
  # Kullback-Leibler divergence of freqs1 compared to freqs2
  return(sum(freqs1 * log(freqs1/freqs2), na.rm=TRUE))  # 0 * 1/0 returns NaN
}


# Measure the divergence from the initiation point.
date_PL100_stability <- function(subtree, subcalib, amin, amax) {
  #age.start.bounds) {
  # age.start.bounds: data.frame(amin=, amax=) with rownames in subtree$node.label
  
#add_starting_age <- function(subtree, subcalib, a) {
  ntips <- Ntip(subtree)
  unknown <- which( !(subtree$node.label %in% rownames(subcalib)) )
  if( length(unknown)>1 ) {
    stop("More than 1 unknown nodes. You have to skip this stability test")
    #unknown <- unknown[-1] ## Assuming this is because the root is a dup (in robusts only)
  }

  newsubcalib <- merge(subcalib,
                       data.frame(node=(ntips+unknown),
                                  age.min=amin,  # Because nlminb fails on NA in LOWER/UPPER
                                  age.max=amax,
                                  soft.bounds=FALSE,
                                  age.start=amin),
                       all=TRUE)  # This resets row.names.
  unknown_row <- which( newsubcalib$node == (ntips+unknown) )  # Optional...
  cat("Unknown row: nb", unknown_row, "\n")

  starting <- seq(amin, amax, length.out=100)
  dates <- numeric(100)
  i <- 1
  for(a in starting) {
    cat(a, "")
  
    newsubcalib[unknown_row, "age.start"] <- a
    dated <- date_PL100(subtree, newsubcalib, quiet=F)
    dates[i] <- branching.times(dated)[unknown]
    i = i+1
  }
  cat('\n')
  #chronograms <- lapply(seq(amin, amax, length.out=100), function(a) {
  #                        date_PL100()
  #                      })
  return(dates)
}

date_PL100_stability_summary <- function(subtree, subcalib, amin, amax) {
  dates <- date_PL100_stability(subtree, subcalib, amin, amax)
  bins <- seq(amin, amax, length.out=11)
  return(KL.div(hist(dates, bins, plot=FALSE)$density,  # density
                0.1))
}
### NOTE: the best lambda value (penalty roughness) must be chosen by
###       cross-validation. Testing multiple subtrees to check the prediction on
###       the removed nodes
#date_bestPL <- function(tree, calibration, minl=0, maxl=2, lstep=0.1) {
#  all_chro <- lapply(seq(minl, maxl, lstep), function(lambda) {
#                        chronos(tree, lambda=lambda, calibration=calibration)})
#  #all_lambda <- sapply(all_chro, function(chro) {attr(chro, "PHIIC", exact=T)$lambda})
#  all_ploglik <- sapply(all_chro, function(chro) {attr(chro, "ploglik", exact=T)})
#
#  return(all_chro[which.min(all_ploglik)])
#
#}

# Function that fills an age at one row (i.e one edge). To use in `apply`.
fill.dates <- function(edgedatum, datation.env=new.env(), scale=1) {
  edge    <- as.character(edgedatum[1:2])
  elength <- edgedatum[3]
  datation.env$dated[edge[2], "age"] <- datation.env$dated[edge[1], "age"] - scale * elength
}

chronogram2table <- function(chronogram,
                             root.age,
                             leaf.age=0,
                             datation.env=new.env(), 
                             rescale=FALSE)
{
  # also see ape::branching.times
  ntips <- Ntip(chronogram)
  rootnb <- ntips + 1
  nnodes <- Nnode(chronogram)
  tot_nodes <- ntips + nnodes
  # either write calibrated nodes too, or remove them from the edges beforehand
  #dated <- data.frame(row.names=seq(1, ntips+nnodes),
  #                    name=c(chronogram$tip.label, chronogram$node.label), 
  #                    age=numeric(ntips+nnodes))
  dated <- data.frame(row.names=seq(rootnb, tot_nodes),
                      name=chronogram$node.label,
                      age=numeric(nnodes))
  dated[as.character(rootnb), "age"] <- root.age
  assign("dated", dated, envir=datation.env)
  edge_data <- cbind(chronogram$edge, chronogram$edge.length)
  scale <- ifelse(rescale, root.age - leaf.age, 1)
  edge_data <- edge_data[edge_data[,2] > ntips,,drop=FALSE]
  apply(edge_data, 1, fill.dates, datation.env, scale)
  return(datation.env$dated)
}


### DEPRECATED
process_line <- function(line, calibration, outfile,
                         date_func=date_PL,
                         #date_func_list=list(
                         age_col="age_dist", datation.env=new.env(),
                         col.names=TRUE, ...) {
  count_iter <<- ifelse(exists("count_iter"), count_iter + 1, 0)
  rescale <- identical(date_func, date_MPL)
  tree <- read.tree(text=line)
  leaf_ages <- calibration[tree$tip.label, age_col]
  # In case tips are not in the table, ignore (assume leaf age is zero)
  leaf_ages[is.na(leaf_ages)] <- 0
  test1 <- all(leaf_ages == leaf_ages[1])
  subcalib <- subset_calibration(calibration, tree, age_col=age_col, leaf_age=leaf_ages[1])
  test2 <- nrow(subcalib) < tree$Nnode  # There are remaining nodes to date.
  test3 <- nrow(subcalib) > 0
  test4 <- all(tree$edge.length > 0)
  cat("\r", count_iter, tree$node.label[1], "     ")
  if (test1 && test2 && test3 && test4) {

    dtree <- tryCatch(date_func(tree, calibration=subcalib, ...),
                            error=function(e) {
                                    cat(" Got Error:", file=stderr())
                                    if (e$message == "NA/NaN gradient evaluation") {
                                      cat("known error", e$message, file=stderr())
                                      warning(e)
                                      return(e)
                                    } else {
                                      cat("unknown", file=stderr())
                                      stop(e)
                                    }})
    if (data.class(dtree) != "chronos") {
      # return the error
      return(c(as.character(as.vector(dtree)), tree$node.label[1], test1, test2, test3))
    }
    root.age <- subcalib[tree$node.label[1], "age.min"]
    
    dated <- chronogram2table(dtree, root.age, leaf_ages[1],
                              datation.env=datation.env,
                              rescale=rescale)
    dated <- cbind(dated,
                   calibration[as.character(dated$name),
                               c("calibrated", "type", "parent", "taxon", "subgenetree")])
    write.table(dated, outfile, sep='\t', quote=FALSE,
                row.names=FALSE, col.names=col.names)
    return(c("calibrate", tree$node.label[1], test1, test2, test3))
  } else {
    return(c("skip", tree$node.label[1], test1, test2, test3))
  }

}


#date_all <- function(treefile, agefile, outfilename, date_func=date_PL,
#                     age_col="age_dist", n=-1, ...) {
#  require(ape)
#  calib <- load_calibration(agefile)
#  #alltrees <- read.tree(treefile)
#  allnewick <- scan(treefile, what=character(), n=n, sep="\n")
#}

extract_MPL_nodeinfo <- function(dtree) {
  # For outputs of the `chronoMPL` function
  nodeinfo <- do.call(cbind, attributes(dtree)[c("stderr", "Pval")])
  rownames(nodeinfo) <- dtree$node.label
  return(nodeinfo)
}

extract_rate_edgeinfo <- function(dtree, input_tree) {
  # For outputs of the `chronos` function
  edgeinfo <- attr(dtree, "rates")
  #if( data.class(dtree) %in% c("chronos", "phylo") ) {
  if( is.null(edgeinfo) ) {
    edgeinfo <- rep(NA, Nedge(input_tree))
    # Add terminal node name as edge label
    #if( length(edgeinfo)>1 )
    names(edgeinfo) <- c(input_tree$tip.label, input_tree$node.label)[input_tree$edge[,2]]
  } else {
    names(edgeinfo) <- c(dtree$tip.label, dtree$node.label)[dtree$edge[,2]]
  }
  return(edgeinfo)
}

extract_runinfo <- function(dtree) {
  # For outputs of the `chronos` function
  if( inherits(dtree, 'error') ) {
    return(list(message=trimws(paste0(gsub('\n', ' ', dtree, fixed=TRUE), collapse=':'), 'right'),
                ploglik=NA, loglik=NA, PHIIC=NA, niter=NA))
  } else {
    dattr <- attributes(dtree)
    return(list(message=dattr$message,
                ploglik=dattr$ploglik,
                loglik=dattr$PHIIC$logLik,
                PHIIC=dattr$PHIIC$PHIIC,
                niter=dattr$niter))
  }
}

date_methods <- list(PL1=date_PL1, PL100=date_PL100, PL10000=date_PL10000, #NPRS=date_NPRS,
                     C=date_C, R1=date_R1, R100=date_R100, R10000=date_R10000, MPL=date_MPL)

make_error_catcher <- function(date_func) {
  #date_func_name <- paste0(substitute(date_func), collapse='')
  date_func_name <- match.call()[2]
  catch_date_func <- function(tree, calibration, ...) {
    #withCallingHandlers
    # See <https://github.com/cran/ape/blob/master/R/chronos.R> search 'stop('
    # for possible errors.
    dtree <- tryCatch(
                      date_func(tree, calibration, ...),
                      error=function(e) {
                              if (grepl(paste('NA/NaN gradient evaluation',
                                              'cannot find reasonable starting dates',
                                              sep='|'),
                                        e$message)) {
                                #cat("known error", e$message, file=stderr())
                                warning(paste0('Caught expected:',
                                               paste0(gsub('\n', ' ', e, fixed=TRUE),
                                                      collapse=':')))
                                return(e)
                              } else {
                                #traceback()  # No traceback available
                                e$message <- paste0(gsub('\n', ' ', e$message,
                                                         fixed=TRUE),
                                                    '. In function ',
                                                    #substitute(date_func),
                                                    date_func_name,
                                                    ': ', tree$node.label[1])
                                cat('Unexpected error:', as.character(e),
                                    file=stderr())
                                #stop(e)
                                return(e)
                              }})
    #if (data.class(dtree) != "chronos") {
    #  # return the error
    #  return(c(as.character(as.vector(dtree)), tree$node.label[1], test1, test2, test3))
    #}
  }
  return(catch_date_func)
}

branching.times.error <- function(dtree, input_tree) {
    if (! (data.class(dtree) %in% c("chronos", "phylo")) ) {
      return(setNames(rep(NA, Nnode(input_tree)), input_tree$node.label))
    } else {
      return(branching.times(dtree))
    }
}


date_all_methods <- function(line, calibration,
                         date_methods=list(
                                           PL1=date_PL1,
                                           PL100=date_PL100,
                                           PL10000=date_PL10000,
                                           #NPRS=date_NPRS,
                                           C=date_C,
                                           R1=date_R1,
                                           R100=date_R100,
                                           R10000=date_R10000,
                                           MPL=date_MPL),
                         age_col="age_dist") {
  count_iter <<- count_iter + 1  # Not correct with parallel.
  tree <- read.tree(text=line)
  cat("\n", count_iter, tree$node.label[1], "     ")

  if( !is.rooted(tree) ) stop(paste('Tree not rooted!', tree$node.label[1]))

  leaf_ages <- calibration[tree$tip.label, age_col]

  # In case tips are not in the table, ignore (assume leaf age is zero)
  leaf_ages[is.na(leaf_ages)] <- 0
  subcalib <- subset_calibration(calibration, tree, age_col=age_col, leaf_age=leaf_ages[1])
  tests <- c(
             leaves_same_age=all(leaf_ages == leaf_ages[1]),
             any_uncalibrated=(nrow(subcalib) < tree$Nnode),  # There are remaining nodes to date.
             any_calibrated=(nrow(subcalib) > 0),
             not_neg_length=all(tree$edge.length >= 0),  # Might give weird results with 0 lengths.
             no_consecutive_zeros=FALSE)
  #check for 2 consecutive 0 lengths.
  el <- tree$edge.length
  e <- tree$edge
  consecutive_el <- cbind(el[match(e[,1], e[,2])], el)
  tests[5] <- all(rowSums(consecutive_el) > 0, na.rm=TRUE)  # NA for the root distance.

  if ( all(tests[1:4]) ) {

    dtrees <- lapply(date_methods, do.call, list(tree, subcalib))

    dated <- sapply(dtrees, branching.times.error, tree)
    if( is.null(dim(dated)) )
      stop("Not able to simplify2array the branching.times (different number of nodes)")
    # Standard error of node ages and pvalue
    dating_nodeinfo <- extract_MPL_nodeinfo(dtrees[["MPL"]])
    # Extend dated with the tip data.
    dating_nodeinfo <- rbind(dating_nodeinfo,
                             matrix(NA, Ntip(tree), ncol(dating_nodeinfo)))
    dated <- rbind(dated,
                   matrix(0, Ntip(tree), ncol(dated)))
    rownames(dated)[(Nnode(tree)+1):(Ntip(tree)+Nnode(tree))] <- tree$tip.label
    
    # Molecular rates
    dating_edgeinfo <- sapply(dtrees[c("PL1", "PL100", "PL10000", #"NPRS",
                                       "R1", "R100", "R10000")],
                              extract_rate_edgeinfo, tree)
    rateC <- attr(dtrees[["C"]], "rates", exact=TRUE)
    if ( is.null(rateC) ) rateC <- NA
    dating_edgeinfo <- cbind(dating_edgeinfo, C=rateC)
    dating_edgeinfo <- rbind(NA, dating_edgeinfo)
    rownames(dating_edgeinfo)[1] <- tree$node.label[1]
    # Reorder based on the nodes then the tips, as in the other data.frames
    #print(dating_edgeinfo)
    dating_edgeinfo <- dating_edgeinfo[c(tree$node.label, tree$tip.label), ]

    info <- calibration[c(tree$node.label, tree$tip.label),
                        c("calibrated", "type", "parent", "taxon", "subgenetree",
                          "is_outgroup")]
    rownames(info) <- c(tree$node.label, tree$tip.label)
    info[tree$tip.label,"subgenetree"] <- info[tree$node.label[1],"subgenetree"]
    #print(info[tree$node.label[1],"subgenetree"])

    #if ( !all(sapply(list(dating_nodeinfo, dating_edgeinfo, info),
    #                 function(d1, d2) identical(rownames(d1), rownames(d2)),
    #                 dated)
    #         ) )
    #  stop("NOT equal rownames\n")

    dated_info <- data.frame(age=dated, age.MPL=dating_nodeinfo, rate=dating_edgeinfo, info)

    runinfo <- unlist(lapply(dtrees[c("PL1", "PL100", "PL10000",
                                      "R1", "R100", "R10000", "C")],
                             extract_runinfo),
                      recursive=FALSE)#, use.names=FALSE)

    loginfo <- c(list(name=tree$node.label[1], action="calibrate"), tests)
    return(list(log=loginfo, run=runinfo, ages=dated_info))
  } else {
    loginfo <- c(list(name=tree$node.label[1]), action="skip", tests)
    return(list(log=loginfo, run=NULL, ages=NULL))
  }

}

### DEPRECATED
date_all_trees <- function(allnewicks, calib, outfilename, date_func=date_PL,
                           age_col="age_dist", ...) {
  datation.env <- new.env()
  count_iter  <<- 0
  outfile <- file(outfilename, "w")
  cat("name\tage\tcalibrated\ttype\tparent\ttaxon\tsubgenetree\n", file=outfile)
  #cat("\n")
  r <- sapply(allnewicks, process_line, calib, outfile, date_func, age_col,
              datation.env, col.names=FALSE, ..., USE.NAMES=FALSE)
  close(outfile)
  skipped <- sum(r[1,] == "skip")
  calibrated <- sum(r[1,] == "calibrate")
  failed <- ncol(r) - skipped - calibrated
  cat("Calibrated:", calibrated, "Skipped:", skipped, "Failed:", failed, "\n",
      file=stderr())
  return(t(r))
}


### DEPRECATED
run <- function(datasetname, agefile, treefile) {
                #outdir="~/ws2/DUPLI_data85/alignments_analysis/ages",
                #subtrees_src=".dSsubtrees.nwk") {

  ## Penalized likelihood
  #agefile <- "ages/Rodentia_M1_agesdNdSdist.subtreesNoEdit-um2.tsv"
  #treefile <- "trees/Rodentia_M1_dist.subtreesNoEdit.fulltree.nwk"
  #datasetbase <- paste0(outdir, "/", dataset)
  #agefile <- paste0(datasetbase, ".tsv")
  #treefile <- paste0(outdir_trees, datasetbase, subtrees_src)

  calib <- load_calibration(agefile)
  allnewicks <- scan(treefile, what=character(), n=n, sep="\n")

  outfile <- paste0(datasetname, "-C.tsv")
  x_c <- date_all_trees(treefile, agefile, outfile, model='discrete',
                  control=chronos.control(nb.rate.cat=1), quiet=TRUE)

  outfile <- paste0(datasetname, "-PL1.tsv")
  x1 <- date_all_trees(treefile, agefile, outfile, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL100.tsv")
  x100 <- date_all_trees(treefile, agefile, outfile, lambda=100, quiet=TRUE)

  outfile <- paste0(datasetname, "-PL10000.tsv")
  x10000 <- date_all_trees(treefile, agefile, outfile, lambda=10000, quiet=TRUE)
  
  # Penalized relaxed
  outfile <- paste0(datasetname, "-PLR.tsv")
  x_r <- date_all_trees(treefile, agefile, outfile, model='relaxed', quiet=TRUE)

  ## MPL
  outfile <- paste0(datasetname, "-MPL.tsv")
  xMPL <- date_all_trees(treefile, agefile, outfile, date_MPL,
                         se=FALSE, test=FALSE)
}


date_all_trees_all_methods <- function(datasetname, agefile, treefile, ncores=6, n=-1) {
  calib <- load_calibration(agefile)
  allnewicks <- scan(treefile, what=character(), n=n, sep="\n")
  outfile_ages <- paste0(datasetname, ".chronos-ages.tsv")
  outfile_runs <- paste0(datasetname, ".chronos-runs.tsv")
  outfile_logs <- paste0(datasetname, ".chronos-logs.tsv")
  count_iter <<- 0

  # Name each tree by its root node name.
  rootname_pattern <- "^.*\\)([A-Za-z.0-9]+):[0-9.-][^(),]+;$"
  names(allnewicks) <- strcapture(rootname_pattern,
                                  allnewicks,
                                  proto=data.frame(root=character()))$root

  date_methods <- list(PL1=date_PL1,
                       PL100=date_PL100,
                       PL10000=date_PL10000,
                       #NPRS=date_NPRS,
                       C=date_C,
                       R1=date_R1,
                       R100=date_R100,
                       R10000=date_R10000,
                       MPL=date_MPL)
  try_date_methods <- lapply(date_methods, make_error_catcher)
  try_date_all <- function(newickline, ...) {
    tryCatch(date_all_methods(newickline, ...),
             error=function(e) {
               rootname <- strcapture(rootname_pattern,
                                      newickline,
                                      proto=data.frame(root=character()))$root
               e$message <- paste0(gsub('\n', ' ', e$message, fixed=TRUE),
                                   ': ', rootname)
               #traceback()
               cat('Error in at least one dating method:', as.character(e),
                   file=stderr())
               stop(e)
             })
  }

  if( ncores == 1 ) {
    out <- lapply(allnewicks, try_date_all, calib, try_date_methods)
  } else {
    cl <- makeForkCluster(ncores, outfile="")
    out <- parLapply(cl, allnewicks, try_date_all, calib, try_date_methods)
    stopCluster(cl)
  }

  out_logs <- do.call(rbind.data.frame, lapply(out, `[[`, "log"))
  out_runs <- do.call(rbind, lapply(out, `[[`, "run"))
  #names(out) <- NULL  # Otherwise, those names will prefix each row name.
  out_ages <- do.call(rbind, c(lapply(unname(out), `[[`, "ages"), deparse.level=0))  # Do not check row names
  #TODO: Do not append 'age' before 'MPL.Pval' and 'MPL.stderr'

  # Just for a shorter representation and lighter files.
  out_logs[3:ncol(out_logs)] <- apply(out_logs[3:ncol(out_logs)], 2,
                                      factor,
                                      levels=c("TRUE", "FALSE"),
                                      labels=c('T', 'F'))

  cat("Total:", nrow(out_logs),
      "Calibrated:", sum(out_logs$action == "calibrate"),
      "Skipped:", sum(out_logs$action == "skip"),
      "Failed:", sum(apply(out_runs, 1,
                             function(r){any(grepl('^Error', r), na.rm=TRUE)})
                      ),
      "\n",
      file=stderr())

  write.table(out_logs, outfile_logs, quote=FALSE, sep="\t", na="", row.names=FALSE)
  write.table(out_runs, outfile_runs, quote=FALSE, sep="\t", na="")
  write.table(out_ages, outfile_ages, quote=FALSE, sep="\t", na="")
  #return(list(out_ages, out_runs, out_logs))
}

if(!interactive()) {
  #options(error=function(){warnings();traceback();return(1)})

  args <- c("Simiiformes_m1w04_ages.subtreesGoodQualO2-ci",
            "Simiiformes_m1w04_ages.subtreesGoodQualO2-ci-um2.tsv",
            "Simiiformes_m1w04_ages.subtreesGoodQualO2-ci.dSsubtrees.nwk",1,35)
  args <- c("ages/Simiiformes_m1w04_ages.subtreesGoodQualO2-cCatarrhini",
            "ages/Simiiformes_m1w04_ages.subtreesGoodQualO2-cCatarrhini-um1.new.tsv",
            "trees/Simiiformes_m1w04_ages.subtreesGoodQualO2-cCatarrhini.dSfulltrees.nwk",
            1, 35)
  args <- commandArgs(trailingOnly=TRUE)
  if( !(length(args) %in% 3:5) ) {
    stop("Wrong number of arguments", usage)
  }
  #datasetname <- "Simiiformes_m1w04_ages.subtreesCleanO2-um2-withSG"
  #run(datasetname)
  datasetname <- args[1]  # base of the outputs
  agefile <- args[2]      # calibrations
  treefile <- args[3]
  ncores <- ifelse(length(args) >= 4, as.integer(args[4]), 1)
  nlines <- ifelse(length(args) >= 5, as.integer(args[5]), -1)

  date_all_trees_all_methods(datasetname, agefile, treefile, ncores, nlines)
  warnings()
}
