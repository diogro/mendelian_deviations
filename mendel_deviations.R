library(ggplot2)
library(reshape2)

source('./read_mouse_data.R')

#chrs <- read.csv('data/genotypes/chrom3.csv')
pedigree <- read.csv('data/F3_pedigree.csv')
families <- split(pedigree, list(pedigree$Dam, pedigree$Sire), drop=T)

multi_prob <- list(LLSS = c(0, 0, 1, 0), SSLL = c(0, 1, 0, 0),
                   SLLL = c(0, 1, 0, 1), LSLL = c(0, 1, 0, 1),
                   SLSL = c(1, 1, 1, 1), SSSL = c(1, 1, 0, 0),
                   LSSL = c(1, 1, 1, 1), SSLS = c(1, 1, 0, 0),
                   LSLS = c(1, 1, 1, 1), SLLS = c(1, 1, 1, 1),
                   SSSL = c(1, 1, 0, 0), LLLL = c(0, 0, 0, 1),
                   SLSS = c(1, 0, 1, 0), SSSS = c(1, 0, 0, 0),
                   LSSS = c(1, 0, 1, 0), LLSL = c(0, 0, 1, 1),
                   LLSL = c(0, 0, 1, 1), LLLS = c(0, 0, 1, 1))

calc_code <- function(x, ret_vector=F) {
    if (x[1] == 1) {
        if (!ret_vector) {
            return('SS')
        } else {
            return(c(1,0,0,0))
        }
    }
    if (x[1] == -1) {
        if (!ret_vector) {
            return('LL')
        } else {
            return(c(0,0,0,1))
        }
    }
    if (x[2] == 1) {
        if (x[3] == 1) {
            if (!ret_vector) {
                return('SL')
            } else {
                return(c(0,0,1,0))
            }
        } else {
            if (!ret_vector) {
                return('LS')
            } else {
                return(c(0,1,0,0))
            }
        }
    }
}

calc_locus_cols <- function(id, locus, chrs) {
    locus_col <- c(paste0("A", locus), paste0("D", locus), paste0("I", locus))
    return (chrs[chrs$ID == id, locus_col])
}

calc_litter_prob <- function(g, litter, sire_id, dame_id, chrs) {
    locus_col <- c(paste0("A", g), paste0("D", g), paste0("I", g))
    sire_locus <- chrs[chrs$ID == sire_id, locus_col]
    dame_locus <- chrs[chrs$ID == dame_id, locus_col]
    sire_code <- calc_code(sire_locus)
    dame_code <- calc_code(dame_locus)
    parents_code <- paste(sire_code, dame_code, sep="")
    classes_prob <- multi_prob[[parents_code]]
    sons_classes <- Reduce("+", llply(litter$ID, function (id) calc_code(calc_locus_cols(id, g, chrs), ret_vector=T)))
    return(dmultinom(x=sons_classes, prob=classes_prob, log=T))
}

get_family <- function(current_family, chrs){
    parents <- strsplit(current_family, "\\.")[[1]]
    sire <- as.numeric(parents[1])
    dame <- as.numeric(parents[2])
    litter  <- chrs[chrs$ID %in% families[[current_family]]$ID,]
    list(litter = litter, sire = sire, dame = dame)
}

get_family_prob <- function(current_family, chrs, nlocus) {
    family <- get_family(current_family, chrs)
    if (dim(family$litter)[1] == 0) {
        return(NA)
    }
    return (aaply(1:nlocus, 1, calc_litter_prob, family$litter, family$sire, family$dame, chrs))
}

re_calc_impossible <- function(fuck_up_id, nlocus){
    family = get_family(fuck_up_id, chrs)
    good_litter <- aaply(1:dim(family$litter)[1], 1, function(ind) all(is.finite(aaply(1:nlocus, 1, calc_litter_prob, family$litter[ind,], family$sire, family$dame, chrs))))
    aaply(1:nlocus, 1, calc_litter_prob, family$litter[good_litter,], family$sire, family$dame, chrs)
}

runCromossome <- function(cromossome){
    chrs <- mouse_gen[[cromossome]]
    nlocus = (length(chrs)-1)/3
    families_probs <- llply(names(families), get_family_prob, chrs, nlocus)
    names(families_probs) <- names(families)
    families_probs <- families_probs[!is.na(families_probs)]
    fuck_up <- names(families_probs[!laply(families_probs, function(x) all(is.finite(x)))])
    #TODO: if() nao tenha ninguem zoado
    families_probs[fuck_up] <- llply(fuck_up, re_calc_impossible, nlocus)
    return(families_probs)
}

families_probs <- runCromossome(3)
chrs <- mouse_gen[[3]]
fuck_up_id <- fuck_up[[1]]

df_probs <- ldply(families_probs)
m_probs <- melt(df_probs)
#names(m_probs)
ggplot(m_probs, aes(variable, value, group = variable, color = .id)) + geom_jitter() 
