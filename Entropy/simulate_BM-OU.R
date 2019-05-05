require("optparse")
require("ouch")
require("phytools")
require("geiger")
require("reshape2")



option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="tree file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.csv", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tree)){
  print_help(opt_parser)
  stop("Tree file path should be supplied.n", call.=FALSE)
}

# extract input files
tree <- read.tree(opt$tree)
tree <- force.ultrametric(tree)

# start simulating the ltr
N = 100

chi_sqr <- vector()
simulation <- vector()

# simulate only bm based
for(i in 1:N) {
  sim_bm <- fastBM(tree, sig2 = 0.01)
  df <- data.frame(melt(sim_bm))
  nc <- name.check(tree,df)
  if (nc != "OK"){
    tree <- drop.tip(tree,nc$tree_not_data)
    df <- df[! df$node_name %in% nc$data_not_tree,]
  }
  ot <- ape2ouch(tree)
  otd <- as(ot,"data.frame")
  df$labels <- rownames(df)
  otd <- merge(otd,df,by="labels",all=TRUE)
  rownames(otd) <- otd$nodes
  
  ### now remake the ouch tree
  ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
  bm <- brown(tree=ot,data=otd['value'])
  otd$regimes <- as.factor("global")
  ou <- hansen(tree=ot, data=otd[c('value')], regimes=otd["regimes"], fit=TRUE, sqrt.alpha=1, sigma=1, maxit=10000)
  chi.square <- as.numeric(2*(logLik(ou) - logLik(bm)))
  
  chi_sqr[i] <- chi.square
  simulation[i] <- i
  
}


sim.results <- data.frame(simulation, chi_sqr)
write.table(sim.results, opt$out, sep = ",", col.names = T, append = T, row.names=FALSE)

print("Done")






