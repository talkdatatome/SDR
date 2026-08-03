#' PermArr - Permutation Array from Supplement S.1
#' 
#' Permutation array generated as the code description in the first part of Supplement S.1.
#' 
#' @format A matrix with 78 rows and 10 columns. The first column is a subset of possible rows from the Hadamard matrix, defined here as \code{c(2:40, 42:80)}. Each additional column a permutation of the first column.
#' @source Supplement S.1 from Slud and Trudell (in preparation).
"PermArr"

#' PrmArr0 - Permutation Array used for default values
#' 
#' Permutation array used as a default value for functions when not otherwise specified.
#' 
#' @format Similar to `PermArr`, a matrix with 78 rows and 10 columns. The first column is a subset of possible rows from the Hadamard matrix. In this case, the first column is \code{1:78}. Each additional column a permutation of the first column.
#' @source From Slud and Trudell (in preparation TODO).
"PrmArr0"

