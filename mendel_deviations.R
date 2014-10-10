chrs <- read.csv('data/genotypes/chrom1.csv')

multi_prob <- list(LLSS = c(0, 0, 1, 0), SSLL = c(0, 1, 0, 0),
                   SLLL = c(0, 1, 0, 1), LSLL = c(0, 1, 0, 1),
                   SLSL = c(1, 1, 1, 1), SSSL = c(1, 1, 0, 0),
                   LLLL = c(0, 0, 0, 1), SSSS = c(1, 0, 0, 0))

calc_code <- function(x) {
    if (x[1] == 1) {
        return('SS')
    }

    if (x[1 == -1]) {
        return('LL')
    }

    if (x[2] == 1) {
        if (x[3] == 1) {
            return('SL')
        } else {
            return('LS')
        }
    }
}

prob_family <- function(parent, mother, sons) {
    parent_code <- calc_code(parent)
    mother_code <- calc_code(mother)
    parents_code <- paste(parent_code, mother_code, sep="")
    classes_prob <- multi_prob[parents_code]

#   sons_classes <- MAGIC

    return dmultinom(x=sons_classes, prob=classes_prob)
}
