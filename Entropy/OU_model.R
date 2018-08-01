# this script is for running the OU model on a tree as described in 
# A.A. King's "Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution"

library("optparse")
library("ouch")
library("phytools")

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
tree <- drop.tip(tree,nc$tree_not_data)

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
h1 <- hansen(
  tree=ot,
  data=otd[c(opt$value)],
  regimes=otd["regimes"],
  fit=TRUE,
  sqrt.alpha=1,
  sigma=1,
  maxit=10000
)


h1_res <- summary(h1)
b1_res <- summary(b1)

b1_ll = b1_res$loglik
h1_ll = h1_res$loglik

# was in use when i stupidly tought that the result is minus 2 log likelihood
#chi.square <- b1_ll - h1_ll

# calculate 2glr statistic
chi.square <- 2*(h1_ll - b1_ll)

alpha = h1_res$alpha
#optima = h1_res$optima$value
ou_sigma = h1_res$sigma.squared[1]
bm_sigma = b1_res$sigma.squared[1]

pval = pchisq(chi.square, df=h1_res$dof - b1_res$dof, lower.tail=FALSE)

names <- c('alpha', 'chiSquare', 'OUsigma', 'BMsigma', 'Pvalue')
values <- c(alpha, chi.square, ou_sigma, bm_sigma, pval)

result <- data.frame(statistics=names, valus=values)
write.table(result, opt$out, sep = ",", col.names = T, append = T, row.names=FALSE)

print(result)
print("Done")






