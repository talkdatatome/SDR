#' PermArr - Permutation Array for Row Assignment
#' 
#' Permutation array generated with an initial set of indices of 2-40 and 42-80, in the first column. Every subsequent column is the sequence determined by a stride of size equal to the column index, indexing starting at 1 for the first column.
#' 
#' @references
#' Slud, E., & Trudell, T. (2026). SDR Variance Estimates in Small Domains. arXiv preprint arXiv:2608.17353. doi:10.48550/arXiv.2608.17353
"PermArr"

#' PrmArr0 - Permutation Array used for defaults 
#' 
#' Permutation array used as a default value for functions when not otherwise specified. 
#' 
#' @format Similar to `PermArr`, a matrix with 78 rows and 10 columns. The indices are simply mapped from 2:40,42:80 to 1:39,42:80.
"PrmArr0"

