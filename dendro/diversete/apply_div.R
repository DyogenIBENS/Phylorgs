#!/usr/bin/env Rscript

# RUN THIS ANALYSIS FROM: ~/ws2/DUPLI_data90/div-VS-dup

library(magrittr)
library(parallel)

library(nlme)

library(ape)
library(geiger)
library(RPANDA)
library(diversitree)

#param_name <- "listens90size10"

backup_cleanup <- function() {
  # BACK UP AND CLEAN UP CURRENT WORKSPACE
  if( exists("param_str", envir=.GlobalEnv) ){
    # Back up
    if( !exists("analyses", envir=.GlobalEnv) ) analyses <<- list()
    analyses[[param_str]] <- new.env()

    important_var <- c("param_str", "all_stats", "div_stats")

    for( var in important_var){
      if( exists(var, envir=.GlobalEnv) )
        assign(var, get(var, envir=.GlobalEnv), envir=analyses[[param_str]])
    }
    # Clean up
    allvars <- ls(all.names=TRUE, envir=.GlobalEnv)
    allvars <- allvars[allvars != "analyses"]
    rm(list=allvars, envir=.GlobalEnv)
  }
}

# SINGLE CONFIGURATION PARAMETER
param_str <- "age150-size10"

param_suffix <- paste0("listens90", param_str)
source_dir <- "~/ws2/"
source_dir2 <- "~/ws7/"
div_path   <- "databases/timetree/Opisthokonta-"
dup_path   <- "DUPLI_data90/"

#workdir <- paste0(source_dir, dup_path, "div-VS-dup")
#setwd(workdir)

if( interactive() ) {
  subtreefile  <- paste0(source_dir, div_path, param_suffix, ".subtrees.nwk")
  divtablefile <- paste0(source_dir, div_path, param_suffix, ".tsv")
  duptablefile <- paste0(source_dir, dup_path, "event_rates-", param_str, ".tsv")
  args <- c(subtreefile, divtablefile, duptablefile)
} else {
  args <- commandArgs(trailingOnly=TRUE)
  subtreefile <- args[1]
}

if(length(args)>1) {
  divtablefile  <- args[2]
  duptablefile  <- args[3]
  clade.div.data <- read.delim(divtablefile, row.names=1, header=TRUE)
  clade.dup.data <- read.delim(duptablefile, row.names=1, header=TRUE)
} else {
  #divtablefile <- NULL
  clade.div.data <- NULL
  clade.dup.data <- NULL
}

fix_single_nodes <- function(lines){
  # newick line with a single node and no parenthesis can't be read by read.tree
}

leaf_heights <- function(tree) {
  cum_edge_length <- tree$edge.length
  for(e in seq_along(cum_edge_length)){
      edge <- tree$edge[e,]
      cur_dist <- tree$edge.length[e]
      prev_edge_pos <- tree$edge[,2] == edge[1]
      prev_height <- ifelse(any(prev_edge_pos),
                            cum_edge_length[prev_edge_pos],
                            0)
      cum_edge_length[e] <- prev_height + cur_dist
  }
  leaf_edge_pos <- which(tree$edge[,2] %in% 1:Ntip(tree))
  return(cum_edge_length[leaf_edge_pos])
  # similar to phylobase::nodeHeights:
  #return(phylobase::nodeHeights(tree)[tree$edge[,2] %in% 1:Ntip(tree), 2])
}

extend_to_ultrametric <- function(tree) {
  # Extend tree leaves to the maximum of the leaf-to-root distance.
  # Return the maximum extension performed
  lheights <- leaf_heights(tree)
  maxh <- max(lheights)
  cat(paste0("Max extension: ", maxh - min(lheights), "\n"), file=stderr())
  leaf_edge_pos <- which(tree$edge[,2] %in% 1:Ntip(tree))
  tree$edge.length[leaf_edge_pos] <- tree$edge.length[leaf_edge_pos] + (maxh - lheights)
  return(tree)
}

summary_posterior <- function(samples, paramname, probaname="p",
                              CIrange=c(0.05, 0.95), burnin_prop=0.1){
  # Computes the cumulative distribution function of the parameter, given a
  # mcmc run ('samples'). 'samples' must contain the columns 'paramname'
  # (parameter) and 'probaname' (the posterior probability of the parameter)
  # Return the value of maximum posterior proba and the 5% interval.
  samples <- samples[-(1:floor(nrow(samples) * burnin_prop)),]
  #max_posterior_pos <- which(samples[,probaname] == max(samples[,probaname]))
  #o_samples <- samples[order(samples[,paramname]),]
  
  #tot_area <- sum(diff(o_samples[,paramname]) * o_samples[-1, probaname])
  ## cumulative distribution function
  # cdf <- cumsum(diff(o_samples[,paramname]) * o_samples[-1, probaname]) / tot_area
  #CI_lower_pos <- which(cdf >= 0.05)[1]
  #CI_upper_pos <- which(cdf >= 0.95)[1]
  #CI <- o_samples[c(CI_lower_pos, CI_upper_pos)+1, paramname]
  #best <- samples[max_posterior_pos, paramname]
  post_density <- density(samples[,paramname])
  best <- post_density$x[which.max(post_density$y)]
  CI <- quantile(samples[,paramname], CIrange)
  out <- list(best=best, CI=CI)
  #names(out) <- if(outname){paste0(outname, ".", c("best", "CI"))}{outname}
  return(out)
}


get_div_stats <- function(tree, metadata) {
  stderrf <- stderr()
  if(!is.ultrametric(tree)) tree <- extend_to_ultrametric(tree)
  treename <- tree$node.label[1]
  cat("- Processing", treename, file=stderrf)
  if(Ntip(tree) <= 2) {
    cat("Ntip <= 2. Skip.\n", file=stderrf)
    return()
  }
  sampling_f <- metadata[treename, "ncbi_sp_sampling"]

  out <- list()
  # In Ape (maximum likelihood)
  cat(" Ape... ", file=stderrf)
  BD_0 <- birthdeath(tree)
  out$b0 <- as.numeric(BD_0$para["b-d"] / (1 - BD_0$para["d/b"]))
  out$d0 <- as.numeric(BD_0$para["d/b"] * out$b0)
  out[["b0-d0"]] <- as.numeric(BD_0$para["b-d"])
  out[["d0/b0"]] <- as.numeric(BD_0$para["d/b"])
  out[["CI_b0-d0"]] <- BD_0$CI["b-d",]
  out[["CI_d0/b0"]] <- BD_0$CI["d/b",]

  out$gamma <- gammaStat(tree)
  # In Geiger (moment estimators)
  cat("Geiger... ", file=stderrf)
  tree$root.edge <- 0         # otherwise will compute the stem estimate
  out$net_d_ms <- bd.ms(tree)
  out$net_d_ms_missing <- bd.ms(tree, missing=round(Ntip(tree) * (1 - sampling_f)))
  out$net_d_km <- bd.km(tree) # missing not implemented for Kendall-Moran
  # In RPANDA
  cat("RPANDA... ", file=stderrf)
  cst.rate <- function(t,r){r[1]}
  tot_time <- max(leaf_heights(tree))
  BD_1 <- fit_bd(tree, tot_time, cst.rate, cst.rate, out$b0, out$d0,
                 cst.lamb=TRUE, cst.mu=TRUE)
  # With incomplete sampling
  BD_2 <- fit_bd(tree, tot_time, cst.rate, cst.rate, out$b0, out$d0, sampling_f,
                 cst.lamb=TRUE, cst.mu=TRUE)
  out$ml.b1 <- BD_1$lamb_par
  out$ml.d1 <- BD_1$mu_par
  out$ml.sampl.b2 <- BD_2$lamb_par
  out$ml.sampl.d2 <- BD_2$mu_par
  # In diversitree (MCMC)
  cat("Diversitree..", file=stderrf)
  bd_func <- make.bd(tree)
  bd_func_sampl <- make.bd(tree, sampling.f=sampling_f)
  # Max-likelihood
  # TODO: catch warnings (improper likelihood convergence)
  cat("ML...", file=stderrf)
  BD_3 <- find.mle(bd_func, c(out$b0, out$d0), method="optim") #condition.surv=F
  BD_4 <- find.mle(bd_func_sampl, c(out$b0, out$d0), method="optim")
  # MCMC
  cat("MCMC...\n", file=stderrf)
  samples_5 <- mcmc(bd_func, as.numeric(c(out$b0, out$d0)), nsteps=10000, print.every=0,
                    w=c(0.1, 0.1))
  samples_6 <- mcmc(bd_func_sampl, as.numeric(c(out$b0, out$d0)), nsteps=10000, print.every=0,
                    w=c(0.1, 0.1))

  out <- c(out, ml2=coef(BD_3), ml2.sampl=coef(BD_4),
           mcmc=summary_posterior(samples_5, 'lambda'),
           mcmc=summary_posterior(samples_5, 'mu'),
           mcmc.sampl=summary_posterior(samples_6, 'lambda'),
           mcmc.sampl=summary_posterior(samples_6, 'mu'))

  return(unlist(out))
}

load_subtrees <- function(subtreefile) {
  subtreelines  <- scan(subtreefile, what=character(), sep="\n")
  not_single_node_lines <- grepl('(', subtreelines, fixed=TRUE)
  subtreelines <- subtreelines[not_single_node_lines]

  subtrees <- read.tree(text=subtreelines, keep.multi=TRUE)
  subtrees <- subtrees[sapply(subtrees, Ntip) > 2]
  names(subtrees) <- sapply(subtrees, function(tree){tree$node.label[1]})
  return(subtrees)
}

analyze_div <- function(subtrees, save=TRUE) {
  #div_stats <- sapply(subtrees, get_div_stats, clade.data)
  cl <- makeForkCluster(max(1, detectCores()-2), outfile="tmp_outfile.txt")
  tmp_output <- parSapplyLB(cl, subtrees, get_div_stats, clade.div.data,
                            simplify=TRUE)
  div_stats <- data.frame(t(tmp_output))
  stopCluster(cl)
  #cat(readLines("tmp_outfile.txt"), sep="\n")

  if(save)
    write.table(div_stats, paste0("div_stats-", param_suffix, ".tsv"), sep='\t', quote=F)

  return(div_stats)
}

load_clade_converter <- function() {
  clade.converter <- read.delim(paste0(source_dir, dup_path, "ens90", param_str,
                                       ".txt"),
                                header=FALSE, col.names=c('div', 'dup'),
                                as.is=TRUE)
  # complete missing values:
  no.dup.name <- clade.converter$dup == ''
  clade.converter$dup[no.dup.name] <- gsub('_', ' ',
                                           clade.converter$div[no.dup.name])
  no.div.name <- clade.converter$div == ''
  clade.converter$div[no.div.name] <- gsub(' ', '_',
                                           clade.converter$dup[no.div.name])
  # Convert to tip labels (ape removes spaces)
  clade.converter$duptree <- gsub(' ', '', clade.converter$dup)
  return(clade.converter)
}

fuse_divdup_data <- function(div_stats, clade.dup.data, clade.converter, save=TRUE) {
  # Convert rownames to make them match.
  rownames.dup.data <- rownames(clade.dup.data)
  rownames(clade.dup.data) <- clade.converter$div[match(rownames.dup.data,
                                                        clade.converter$dup)]
  # Compare
  common.clades <- intersect(rownames(div_stats), rownames(clade.dup.data))
  cat("Unconverted div_stats rows:\n")
  cat(setdiff(rownames(div_stats), rownames(clade.dup.data)), '\n')
  cat("Unconverted dup data rows:\n")
  cat(setdiff(rownames(clade.dup.data), rownames(div_stats)), '\n')

  all_stats <- cbind(div_stats[common.clades,], clade.dup.data[common.clades,])
  all_stats$allDup <- all_stats$tandemDup + all_stats$dispDup
  all_stats$allnew <- all_stats$allDup + all_stats$birth
  if(save)
    write.table(all_stats, paste0("all_stats-", param_suffix, ".tsv"),
                sep='\t', quote=F)
  return(all_stats)
}

combine_maindata <- function(all_stats, maintree, clade.converter) {
  # Modifies inplace!
  #maintree$tip.label <- clade.converter$div[match(maintree$tip.label,
  #                                                clade.converter$duptree)]

  maintreedi <- multi2di(maintree)
  maintreedi$tip.label <- clade.converter$div[match(maintreedi$tip.label,
                                                  clade.converter$duptree)]
  # First ensure the correspondance between maintree$tip.label and
  # rownames(all_stats):
  #cat("Ignoring the following branches of the tree:\n")
  #print(setdiff(maintree$tip.label, rownames(all_stats)))
  #cat("Ignoring the following data from all_stats:\n")
  #print(setdiff(rownames(all_stats), maintree$tip.label))

  # subset both data:
  maindata <- geiger::treedata(maintreedi, all_stats, sort=TRUE)
  maindata$data <- data.frame(maindata$data)
  return(maindata)
}

if( !interactive() ) {

  #backup_cleanup()

  subtreefile %>% load_subtrees %>% analyze_div -> div_stats

  load_clade_converter() -> clade.converter
  fuse_divdup_data(div_stats, clade.dup.data, clade.converter) -> all_stats

  # Now the phylogenetic correlation
  maintree <- read.tree(paste0(source_dir, dup_path, "event_rates-", param_str,
                               ".nwk"))
  maindata <- combine_maindata(all_stats, maintree, clade.converter)

  pic.b2        <- pic(maindata$data[,"ml.sampl.b2"], maindata$phy)
  pic.d2        <- pic(maindata$data[,"ml.sampl.d2"], maindata$phy)
  pic.tandemDup <- pic(maindata$data[,"tandemDup"]  , maindata$phy)
  pic.dispDup   <- pic(maindata$data[,"dispDup"]    , maindata$phy)
  pic.allDup    <- pic(maindata$data[,"allDup"]     , maindata$phy)
  pic.allnew    <- pic(maindata$data[,"allnew"]     , maindata$phy)

  pic.cor.tests <- list(
    tandemDup=cor.test(pic.b2, pic.tandemDup),
    dispDup=cor.test(pic.b2, pic.dispDup),
    allDup=cor.test(pic.b2, pic.allDup),
    allnew=cor.test(pic.b2, pic.allnew))

  phycovar.B <- corBrownian(1, maindata$phy)
  phycovar.G <- corGrafen(1, maindata$phy)
  #phycovar.M <- corMartins
  #phycovar.P05 <- corPagel(0.5, maindata$phy)

  gls.tests.B <- list(
    tandemDup=gls(ml.sampl.b2~tandemDup, maindata$data, phycovar.B),
    dispDup=gls(ml.sampl.b2~dispDup, maindata$data, phycovar.B),
    allDup=gls(ml.sampl.b2~allDup, maindata$data, phycovar.B),
    allnew=gls(ml.sampl.b2~allnew, maindata$data, phycovar.B))

  gls.tests.G <- list(
    tandemDup=gls(ml.sampl.b2~tandemDup, maindata$data, phycovar.G),
    dispDup=gls(ml.sampl.b2~dispDup, maindata$data, phycovar.G),
    allDup=gls(ml.sampl.b2~allDup, maindata$data, phycovar.G),
    allnew=gls(ml.sampl.b2~allnew, maindata$data, phycovar.G))

  gls.tests.B.log <- list(
    tandemDup=gls(log(tandemDup)~ml.sampl.b2, maindata$data, phycovar.B),
    dispDup=gls(  log(dispDup)  ~ml.sampl.b2, maindata$data, phycovar.B),
    allDup=gls(   log(allDup)   ~ml.sampl.b2, maindata$data, phycovar.B),
    allnew=gls(   log(allnew)   ~ml.sampl.b2, maindata$data, phycovar.B))
  lapply(gls.tests.B.log, summary)  # check p-values.

  # Obviously the two fish clades are outliers: super diverse, very few dup.
  # Let's remove them and try again.

  tetrapods <- !(rownames(all_stats) %in% c('Percomorphaceae', 'Otophysi'))
  maindata.tetra <- combine_maindata(all_stats[tetrapods,], maintree, clade.converter)

  pic.b2.tetra        <- pic(maindata.tetra$data[,"ml.sampl.b2"], maindata.tetra$phy)
  pic.d2.tetra        <- pic(maindata.tetra$data[,"ml.sampl.d2"], maindata.tetra$phy)
  pic.tandemDup.tetra <- pic(maindata.tetra$data[,"tandemDup"]  , maindata.tetra$phy)
  pic.dispDup.tetra   <- pic(maindata.tetra$data[,"dispDup"]    , maindata.tetra$phy)
  pic.allDup.tetra    <- pic(maindata.tetra$data[,"allDup"]     , maindata.tetra$phy)
  pic.allnew.tetra    <- pic(maindata.tetra$data[,"allnew"]     , maindata.tetra$phy)
  #pic.GL.tetra        <- pic(maindata.tetra$data[,"GL"]         , maindata.tetra$phy)

  pic.cor.tests.tetra <- list(
    tandemDup=cor.test(pic.b2.tetra, pic.tandemDup.tetra),
    dispDup=cor.test(pic.b2.tetra, pic.dispDup.tetra),
    allDup=cor.test(pic.b2.tetra, pic.allDup.tetra),
    allnew=cor.test(pic.b2.tetra, pic.allnew.tetra))
  # None significant, by far.

  phycovar.B.tetra <- corBrownian(1, maindata.tetra$phy)
  phycovar.G.tetra <- corGrafen(1, maindata.tetra$phy)
  nophycovar.tetra <- corBrownian(0, maindata.tetra$phy)
  #phycovar.M <- corMartins
  #phycovar.P05 <- corPagel(0.5, maindata$phy)

  gls.tests.B.tetra <- list(
    tandemDup=gls(ml.sampl.b2~tandemDup, maindata.tetra$data, phycovar.B.tetra),
    dispDup=gls(ml.sampl.b2~dispDup, maindata.tetra$data, phycovar.B.tetra),
    allDup=gls(ml.sampl.b2~allDup, maindata.tetra$data, phycovar.B.tetra),
    allnew=gls(ml.sampl.b2~allnew, maindata.tetra$data, phycovar.B.tetra))
  lapply(gls.tests.B.tetra, summary)  # check p-values.

  gls.tests.G.tetra <- list(
    tandemDup=gls(ml.sampl.b2~tandemDup, maindata.tetra$data, phycovar.G.tetra),
    dispDup=gls(ml.sampl.b2~dispDup, maindata.tetra$data, phycovar.G.tetra),
    allDup=gls(ml.sampl.b2~allDup, maindata.tetra$data, phycovar.G.tetra),
    allnew=gls(ml.sampl.b2~allnew, maindata.tetra$data, phycovar.G.tetra))

  # TODO:
  # 1. Add the generation time;
  # 2. log transform duplication numbers.

  gls.tests.B.log.tetra <- list(
    tandemDup=gls(ml.sampl.b2~log(tandemDup), maindata.tetra$data, phycovar.B.tetra),
    dispDup=gls(ml.sampl.b2~log(dispDup), maindata.tetra$data, phycovar.B.tetra),
    allDup=gls(ml.sampl.b2~log(allDup), maindata.tetra$data, phycovar.B.tetra),
    allnew=gls(ml.sampl.b2~log(allnew), maindata.tetra$data, phycovar.B.tetra))
  # Reverse x/y
  gls.tests.B.log.tetra <- list(
    tandemDup=gls(log(tandemDup)~ml.sampl.b2, maindata.tetra$data, phycovar.B.tetra),
    dispDup=gls(  log(dispDup)~ml.sampl.b2, maindata.tetra$data, phycovar.B.tetra),
    allDup=gls(   log(allDup)~ml.sampl.b2, maindata.tetra$data, phycovar.B.tetra),
    allnew=gls(   log(allnew)~ml.sampl.b2, maindata.tetra$data, phycovar.B.tetra))
  lapply(gls.tests.B.log.tetra, summary)  # check p-values.
  # OK, after this p-hacking, I get dispDup as significant, with slope 0.07
  pvalues <- lapply(gls.tests.B.log.tetra, function(glsres) {summary(glsres)$tTable['ml.sampl.b2','p-value']})
  # But with multiple testing correction?
  p.adjust(pvalues, 'BH')  # -> Doesn't resist it.

  # Non phylogenetic regression
  ols.tests.log.tetra <- list(
    tandemDup=gls(log(tandemDup)~ml.sampl.b2, maindata.tetra$data, nophycovar.tetra),
    dispDup=gls(  log(dispDup)~ml.sampl.b2,   maindata.tetra$data, nophycovar.tetra),
    allDup=gls(   log(allDup)~ml.sampl.b2,    maindata.tetra$data, nophycovar.tetra),
    allnew=gls(   log(allnew)~ml.sampl.b2,    maindata.tetra$data, nophycovar.tetra))
  lapply(ols.tests.log.tetra, summary)  # check p-values.
  # Hmmm still seems to have used phylogenetic covariance. Much different from:
  # Control:
  summary(lm(log(dispDup)~ml.sampl.b2, maindata.tetra$data))
  
  # Double-check: PGLS from caper:
  library(caper)
  caperdata.tetra <- comparative.data(maindata.tetra$phy,
                                      cbind(maindata.tetra$data, clades=rownames(maindata.tetra$data)),
                                      'clades',
                                      vcv=TRUE)
  pglsfit <- pgls(log(dispDup)~ml.sampl.b2, caperdata.tetra)
  # Ouf, les coefs et p-values sont suffisamment proches.

  GL <- read.table(paste0(source_dir2, 'databases/Generation_Length_for_Mammals.csv'),
                                   sep='\t', header=TRUE, as.is=TRUE)
  GL$Calculated_GL_d <- as.numeric(GL$Calculated_GL_d)
  standard.clades <- clade.converter$dup[match(rownames(all_stats),
                                               clade.converter$div)]
  standard.clades %in% c(GL$Order, GL$Family, GL$Genus)
  # Need Orders for:
  Orders <- list(
    Marsupialia=c('Dasyuromorphia', 'Didelphimorphia', 'Diprotodontia', 'Microbiotheria', 'Notoryctemorphia', 'Paucituberculata', 'Peramelemorphia'),
    Atlantogenata=c('Cingulata', 'Pilosa', 'Afrosoricida', 'Hyracoidea', 'Macroscelidea', 'Proboscidea', 'Sirenia', 'Tubulidentata'),
    Insectivora=c('Eulipotyphla'))
  Families <- list(
    Haplorrhini=c('Cercopithecidae', 'Hominidae', 'Hylobatidae', 'Aotidae', 'Atelidae', 'Cebidae', 'Pitheciidae', 'Tarsiidae'),
    Strepsirrhini=c('Daubentoniidae', 'Cheirogaleidae', 'Indriidae', 'Lemuridae', 'Lepilemuridae', 'Galagidae', 'Lorisidae'),
    Hystricognathi=c('Abrocomidae', 'Bathyergidae', 'Capromyidae', 'Caviidae', 'Chinchillidae', 'Ctenodactylidae', 'Ctenomyidae', 'Cuniculidae',
                     'Dasyproctidae', 'Dinomyidae', 'Echimyidae', 'Erethizontidae', # 'Hydrochaeridae',
                     'Hystricidae', 'Diatomyidae', 'Myocastoridae',
                     'Octodontidae', 'Petromuridae', 'Thryonomyidae'),
    Murinae=c('Muridae')  # for my sampling
  )

  clade.GL <- rep(NA, length(standard.clades))
  names(clade.GL) <- standard.clades

  for (clade in standard.clades) {
    if(clade %in% GL$Order) {
      clade.GL[clade] <- mean(GL$Calculated_GL_d[GL$Order==clade], na.rm=TRUE)
    } else if(clade %in% GL$Family) {
      clade.GL[clade] <- mean(GL$Calculated_GL_d[GL$Family==clade], na.rm=TRUE)
    } else if(clade %in% names(Orders)) {
      clade.GL[clade] <- mean(GL$Calculated_GL_d[GL$Order %in% Orders[[clade]]], na.rm=TRUE)
    } else if(clade %in% names(Families)) {
      clade.GL[clade] <- mean(GL$Calculated_GL_d[GL$Family %in% Families[[clade]]], na.rm=TRUE)
    }
  }
  clade.GL['Neognathae'] <- 1000  # I'm making this up.
  all_stats$GL <- clade.GL
  # Rerun the combine_maindata()

  gls.tests.B.log.tetra.GL <- list(
    tandemDup=gls(ml.sampl.b2~GL+log(tandemDup), maindata.tetra$data, phycovar.B.tetra),
    dispDup=gls(ml.sampl.b2~GL+log(dispDup), maindata.tetra$data, phycovar.B.tetra),
    allDup=gls(ml.sampl.b2~GL+log(allDup), maindata.tetra$data, phycovar.B.tetra),
    allnew=gls(ml.sampl.b2~GL+log(allnew), maindata.tetra$data, phycovar.B.tetra))
  lapply(gls.tests.B.log.tetra.GL, summary)
  # OK, GL has no effect on diversification. Somewhat surprising.
}


