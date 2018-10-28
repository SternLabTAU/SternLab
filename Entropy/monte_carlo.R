require("optparse")
require("ouch")
require("phytools")
require("geiger")
require("reshape2")



option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="tree file name", metavar="character"),
  make_option(c("-v", "--value"), type="character", default='k5', 
              help="feature name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.csv", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input file should be supplied.n", call.=FALSE)
}
if (is.null(opt$tree)){
  print_help(opt_parser)
  stop("Tree file path should be supplied.n", call.=FALSE)
}

# extract input files
tree <- read.tree(opt$tree)
df <- read.csv(opt$file, header=TRUE)
rownames(df) <- df$node_name
df <- as.data.frame(df)
nc <- name.check(tree,df)
if (nc != "OK"){
  tree <- drop.tip(tree,nc$tree_not_data)
  df <- df[! df$node_name %in% nc$data_not_tree,]
}


# transform into ouch tree object and pre-process the tree and the trait values

ot <- ape2ouch(tree)
otd <- as(ot,"data.frame")
df$labels <- rownames(df)
otd <- merge(otd,df,by="labels",all=TRUE)
rownames(otd) <- otd$nodes

### now remake the ouch tree
ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

# BM model
b1 <- brown(tree=ot,data=otd[c(opt$value)])

### evaluate an OU model with a single, global selective regime
otd$regimes <- as.factor("global")
h1 <- hansen(tree=ot,data=otd[c(opt$value)],regimes=otd["regimes"],fit=TRUE,sqrt.alpha=1,sigma=1,maxit=10000)

# get parameters before simulating
sigma_ou = as.numeric(summary(h1)$sigma.squared)
sigma_bm = as.numeric(summary(b1)$sigma.squared)
alpha = as.numeric(summary(h1)$alpha)
theta = as.numeric(summary(h1)$optima$k5)
# now we have our own parameters, simulate both OU an BM

N = 100

lr_bm <- vector()
lr_ou <- vector()
simulation <- vector()


for(i in 1:N) {
  
  # simulate bm based
  
  sim_bm <- fastBM(tree, sig2 = sigma_bm)
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
  lr <- as.numeric(2*(logLik(ou) - logLik(bm)))
  
  lr_bm[i] <- lr
  
  # simulate according to ou:
  sim_ou <- fastBM(tree, sig2 = sigma_ou, alpha = alpha, theta=theta)
  df <- data.frame(melt(sim_ou))
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
  lr <- as.numeric(2*(logLik(ou) - logLik(bm)))
  
  lr_ou[i] <- lr
  simulation[i] <- i
  
}


sim.results <- data.frame(simulation, lr_bm, lr_ou)
write.table(sim.results, opt$out, sep = ",", col.names = T, append = T, row.names=FALSE)

print("Done")





