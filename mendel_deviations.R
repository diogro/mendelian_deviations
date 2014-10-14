library(ggplot2)
library(reshape2)

source('./read_mouse_data.R')

#chrs <- read.csv('data/genotypes/chrom3.csv')
pedigree <- read.csv('data/F3_pedigree.csv')
families <- split(pedigree, list(pedigree$Dam, pedigree$Sire), drop=T)

mouse_gen[[1]] <- mouse_gen[[1]][!(mouse_gen[[1]]$ID == 3311 | mouse_gen[[1]]$ID == 2833),]

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

calc_family_prob <- function(g, family, sire_id, dame_id, chrs) {
    locus_col <- c(paste0("A", g), paste0("D", g), paste0("I", g))
    sire_locus <- chrs[chrs$ID == sire_id, locus_col]
    dame_locus <- chrs[chrs$ID == dame_id, locus_col]
    sire_code <- calc_code(sire_locus)
    dame_code <- calc_code(dame_locus)
    parents_code <- paste(sire_code, dame_code, sep="")
    classes_prob <- multi_prob[[parents_code]]
    sons_classes <- Reduce("+", llply(family$ID, function (id) calc_code(calc_locus_cols(id, g, chrs), ret_vector=T)))
    return(dmultinom(x=sons_classes, prob=classes_prob, log=T))
}

runCromossome <- function(cromossome){
    chrs <- mouse_gen[[cromossome]]
    nlocus = (length(chrs)-1)/3
    families_probs <- llply(names(families), function(current_family, chrs) {
                            parents <- strsplit(current_family, "\\.")[[1]]
                            sire_id <- as.numeric(parents[1])
                            dame_id <- as.numeric(parents[2])
                            family  <- chrs[chrs$ID %in% families[[current_family]]$ID,]
                            if (dim(family)[1] == 0) {
                                return(NA)
                            }
                            return (aaply(1:nlocus, 1, calc_family_prob, family, sire_id, dame_id, chrs))
                   }, chrs)
    names(families_probs) <- names(families)
    families_probs <- families_probs[!is.na(families_probs)]
    return(families_probs)
}
families_probs <- runCromossome(2)

df_probs <- ldply(families_probs)
m_probs <- melt(df_probs)
#names(m_probs)
ggplot(m_probs, aes(variable, value, group = variable, color = .id)) + geom_jitter() 
