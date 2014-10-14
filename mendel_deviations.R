library(plyr)

chrs <- read.csv('data/genotypes/chrom1.csv')
pedigree <- read.csv('data/F3_pedigree.csv')
families <- split(pedigree, list(pedigree$Dam, pedigree$Sire), drop=T)
ngenes = 31

multi.prob <- list(LLSS = c(0, 0, 1, 0), SSLL = c(0, 1, 0, 0),
                   SLLL = c(0, 1, 0, 1), LSLL = c(0, 1, 0, 1),
                   SLSL = c(1, 1, 1, 1), SSSL = c(1, 1, 0, 0),
                   LSSL = c(1, 1, 1, 1), SSLS = c(1, 1, 0, 0),
                   LSLS = c(1, 1, 1, 1), SLLS = c(1, 1, 1, 1),
                   SSSL = c(1, 1, 0, 0), LLLL = c(0, 0, 0, 1),
                   SLSS = c(1, 0, 1, 0), SSSS = c(1, 0, 0, 0),
                   LSSS = c(1, 0, 1, 0), LLSL = c(0, 0, 1, 1),
                   LLSL = c(0, 0, 1, 1), LLLS = c(0, 0, 1, 1))

calc.code <- function(x, ret.vector=F) {
    if (x[1] == 1) {
        if (!ret.vector) {
            return('SS')
        } else {
            return(c(1,0,0,0))
        }
    }

    if (x[1] == -1) {
        if (!ret.vector) {
            return('LL')
        } else {
            return(c(0,0,0,1))
        }
    }

    if (x[2] == 1) {
        if (x[3] == 1) {
            if (!ret.vector) {
                return('SL')
            } else {
                return(c(0,0,1,0))
            }
        } else {
            if (!ret.vector) {
                return('LS')
            } else {
                return(c(0,1,0,0))
            }
        }
    }
}

gene.cols <- function(id, gene) {
    gene.col <- c(paste("A", gene, sep=""), paste("D", gene, sep=""), paste("I", gene, sep=""))
    return (chrs[chrs$ID == id, gene.col])
}

calc_family_prob <- function(g, family, sire.id, dame.id) {
    gene.col <- c(paste("A", g, sep=""), paste("D", g, sep=""), paste("I", g, sep=""))
    sire.gene <- chrs[chrs$ID == sire.id, gene.col]
    dame.gene <- chrs[chrs$ID == dame.id, gene.col]
    sire.code <- calc.code(sire.gene)
    dame.code <- calc.code(dame.gene)
    parents.code <- paste(sire.code, dame.code, sep="")
    classes.prob <- multi.prob[[parents.code]]
    sons.classes <- Reduce("+", llply(family$ID, function (id) calc.code(gene.cols(id, g), ret.vector=T)))
    return(dmultinom(x=sons.classes, prob=classes.prob, log=T))
}

families.probs <- llply(names(families), function(n) {
    parents <- strsplit(n, "\\.")[[1]]
    sire.id <- as.numeric(parents[1])
    dame.id <-  as.numeric(parents[2])
    family <- chrs[chrs$ID %in% families[[n]]$ID,]

    if (dim(family)[1] == 0) {
        return(NA)
    }

    return (aaply(1:ngenes, 1, calc_family_prob, family, sire.id, dame.id))
})

names(families.probs) <- names(families)
