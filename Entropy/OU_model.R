# this script is for running the OU model on a tree as described in 
# A.A. King's "Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution"

library("optparse")
library("ouch")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="tree file name", metavar="character"),
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
  stop("Tree foile path shpuld be supplied.n", call.=FALSE)
}

# extract input files
tree <- read.tree(opt$tree)
df <- read.csv(opt$file, header=TRUE)

# transform into ouch tree object and pre-process the tree and the trait values

ot <- ape2ouch(tree)
otd <- as(ot,"data.frame")
df$labels <- df$tree.tip.label
otd <- merge(otd,dat,by="labels",all=TRUE)
rownames(otd) <- otd$nodes

### now remake the ouch tree
ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

# BM model
b1 <- brown(tree=ot,data=otd["values"])

### evaluate an OU model with a single, global selective regime
otd$regimes <- as.factor("global")
h1 <- hansen(
  tree=ot,
  data=otd["values"],
  regimes=otd["regimes"],
  sqrt.alpha=1,
  sigma=1,
  maxit=10000
)


h1_res <- summary(h1)
b1_res <- summary(b1)

b1_ll = b1_res$loglik
h1_ll = h1_res$loglik

lr = b1_ll / h1/ll

names <- c('BMloglik', 'OUloglik', 'LRT')
values <- c(b1_ll, h1_ll,lr)

result <- data.frame(statistics=names, values)
write.table(result, opt$out, sep = ",", col.names = T, append = T)

print(result)
print("Done")






