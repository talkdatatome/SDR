#' Utility program to calculate the sparse representation of the M matrix
#'
#' @description Calculates sparse representation of M matrix for all choices of m, D, and n.
#' 
#'
#' @param n Total Sample Size
#' @param m Number of replicates used
#' @param D D cycling parameter
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @returns Returns a matrix with rows consisting of three elements, an index-pair from 1:n along with the entry for M, for all index-pairs with non-zero entries.
#' @examples
#' mcalc_out <- Mcalc(200)
#' Mmat = array(0, dim=c(200, 200))
#' Mmat[mcalc_out[,1:2]] = mcalc_out[,3]
#' Mmat
#' @export 
Mcalc <- function(n, m=78, D=1, permarr=PrmArr0[,1:D]) {
  ## matrix consisting of all nonzero values j, j', b(j,j')
  ## for j, j' in 1:n and b(j,j') the curly-bracket value in (16)
  ## permarr is m x D matrix of permutations used in D-cycling with D>1
  ## with each column a permutation of 1:m starting with 1.
  nrem = n %% (D*m)
  ncyc = n  %/% (D*m)
  gam = if(n < m) c(1:n, 1) else { 
    if(D==1) c(rep(1:m, ncyc), 1:(1+nrem)) else {
      prmvec = c(permarr)
      c(rep(prmvec, ncyc), prmvec[1:(1+nrem)]) }}
  Mout = NULL
  for(j in 1:n) for (k in 1:n)  {
    jv = gam[j:(j+1)]
    kv = gam[k:(k+1)]           
    if(length(intersect(jv,kv)))  {
      a = sum(jv==kv)
      b = sum(jv==kv[2:1])  
      mjk = 0.5*(a-b)
      Mout = rbind(Mout, c(j,k,mjk)) } }
  Mout  }

#' Produces Expected Relative SDR Bias with permutation tabulation
#' 
#' @description
#' Implementing the same basic formula as in `RelSDRBias` but adds tabulating results over multiple (either all or a random subset) permutations of strata. If maxiter set to 1, then should be identical to `RelSDRBias`.
#'
#' @param parmlst list of parameter vectors for each stratum. Required order is stratum size vector, sample size vector, stratum mean vector, stratum level standard deviation vector, and vector of stratum level domain percentages.
#' @param maxRB upper-bound of relative bias for tabulation
#' @param maxiter if number of strata is above 10, set maximum sampled permutations.
#' @returns list of relative bias prior to permutation, a count tabulation of permutation bins, the true variance, expected domain sample size, finally a vector of the permutation relative bias.
#' @examples
#' parmlst <- list(Nh=rep(1000, 5), 
#'                 nh=rep(25, 5),
#'                 muh=rep(c(1,2), length=5),
#'                 sigmah=rep(1, 5),
#'                 alphah=rep(c(0.1, 0.4), length=5))
#' ans <- PermRelBias(parmlst)
#' print(ans)
#' @export
PermRelBias <-
function(parmlst, maxRB=0.6, maxiter=5e4) {
        ## Cover all permutations if H <= 10 and H! <= maxiter, 
        ##       otherwise sample  maxiter of them.
      ##  maxRB is the upper-bound of relative-bias, below which we tabulate
      ##    the fraction of permuted-strata VSDR estimates in equall spaced intervals
        Nh = parmlst[[1]]
        N = sum(Nh)
        nh = parmlst[[2]]
        H = length(Nh) 
        ph = Nh/N
        muh = parmlst[[3]]
        sgsqh = parmlst[[4]]^2
        alphh = parmlst[[5]]
        Vtru = sum( (ph^2*alphh/nh)*(1-nh/Nh)*(sgsqh + (1-alphh)*muh^2) )             
        alth0 = ph*alphh*muh/nh
        relBpt2 = sum(ph*alphh*(sgsqh + (1-alphh)*muh^2*
                 (1+((nh-1)*Nh/(nh^2*(Nh-1))))))/N
        RelBias = (sum(diff(c(alth0,alth0[1]))^2)/2 + relBpt2)/Vtru
        nD = sum(alphh*nh)
      if(H <= 10 & gamma(H+1) <= maxiter ) {
        tmp = combinat::permn(H)
        RB.num = numeric(gamma(H+1))
        for (i in 1:gamma(H+1)) {
            alth = alth0[tmp[[i]]]
            RB.num[i] = sum(diff(c(alth,alth[1]))^2) }        
      } else {
        unifarr = array(runif(maxiter*H), c(maxiter,H))
        RB.num = numeric(maxiter)
        for (i in 1:maxiter) {
            alth = alth0[order(unifarr[i,])]
            RB.num[i] = sum(diff(c(alth,alth[1]))^2) }        
      }
        RB.vec = (RB.num+ 2*relBpt2)/(2*Vtru)
        countvec = table(cut(RB.vec, breaks=(0:10)*maxRB/10))
      ## count-tbl tabulated VSDR estimates in equally spaced intervals
     ##   and the vector of all Relative Biases from permuted strata is RB.vec
        list(RelBias=RelBias, count.tbl=countvec, Vtru=Vtru, nD=nD, RB.vec=RB.vec) }

#' Calculates SDR variance estimate based on a given stratum permutation
#' 
#' @description
#' Given weighted y vector and a sparse representation of the non-zero entries of the SDR matrix, and an optional stratum permutation, generate the SDR variance estimate.
#' 
#' @param ywt A list of stratum design weighted outcome vectors. Length is equal to number of stratum.
#' @param Mmat An length(ywt) by 3 matrix representing the non-zero indices of the SDR matrix given a D cycling configuration. Can use output from `Mcalc`.
#' @param nhvec vector representing the size of each stratum
#' @param Str.perm an integer vector that provides the permutation order to use in calculation
#' @returns Total estimate under permutation given in Str.perm.
#' @examples
#' nh <- c(3,5,7,9)
#' ywt <- do.call(c, lapply(nh, rnorm, mean=5, sd=1))
#' print(PermSDR(ywt, Mcalc(sum(nh)), nhvec=nh))
#' @export
PermSDR <- function(ywt, Mmat, nhvec, Str.perm = 1:length(nhvec)) {                    
  ### calculates SDR var estimate based on stratum permutation,
  ##  fixing order of y values within stratum
  ## NB Assumes ywt values appear in blocks of size nhvec[h]
  ## Mmat is sparse-matrix representation found in Mcalc
  H = length(Str.perm)        
  cumnh = c(0,cumsum(nhvec))
  ywprm = NULL
  for(h in 1:H) ywprm = c(ywprm, ywt[cumnh[Str.perm[h]] + 
                                       (1:nhvec[Str.perm[h]])])
  ## now use Mmat to find  ysdr; recall Yest does not use SDR info
  Ysdr = sum(ywprm[Mmat[,1]]*ywprm[Mmat[,2]]*Mmat[,3])
  Ysdr   }


#' Calculates correct variance of VSDR variance-estimate components
#'
#' @description For general m, D, and n, calculates coefficients to calculate correct variance for VSDR under simple random sampling from an i.i.d. superpopulation. The coefficients calculated are summed against expressions involving only the first four moments of t he population variables $y_i$. 
#' 
#' @param n total sample size
#' @param m number of replicates used. Defaults to 78 assuming an 80 x 80 Hadamard matrix
#' @param D D cycling parameter. Defaults to 1.
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @returns A list with vectors Coef and Trm. TODO
#' @export
VVsdrEst <- function(n, m=78, D=1, permarr=PrmArr0[,1:D, drop=F]) {
  ### begin by setting n x n Mmat matrix using Mcalc
  Mmat = Mcalc(n, m, D, permarr)
  ## calculate n* variance for SDR-var estimator 
  ## using wt*y from sample of n scaled by 1/n so that sum(wt)=1
  Coef = numeric(4)
  ## coef's apply for D-cycle SDR design respectively to:
  ##  mu4           4th central moment minus var^2 
  ##  mu3*mu        3rd central moment times mean
  ##  sig^2*mu^2     var times mean^2
  ##  sig^4          var^2
  Trm = numeric(4)
  ## coef's calculated (using M_{j,j}=1) in terms of moments of y and
  ## sums of M_{j,j'}*M_{k,k'} over all j,j'k,k' quadruples of 4 types:
  ## type 1: j=j'=k=k'
  ## type 2: 3 indices the same, one different
  ## type 3: 2 pairs (j,j')(j,j') or (j,j')(j',j), with j != j'
  ## type 4: 1 pair (j,j')(j,k) or (j,j')(k,j) or (j,j')(j',k) 
  ##                or (j,j')(k,j')with j,j',k different
  Trm[1] = n
  ## Trm[2] = 4*sum of M_{j,j'}   over all j != j'
  ## Trm[3] = 2*sum of M_{j,j'}^2 over all j != j'
  ## Trm[4] = 4*sum of M_{j,j'}*M_{j,k} over j,j',k all different
  ind = which(Mmat[,1]!=Mmat[,2])
  Mv = Mmat[ind,3]
  M1 = aggregate(Mmat[,3], by=list(Mmat[,1]), sum)$x
  Trm[2] = 4*sum(Mv)
  Trm[3] = 2*sum(Mv^2)
  Trm[4] = 4*(n - sum(Mv^2) - 2*sum(Mmat[,3]) + sum(M1^2))
  Coef[1] = Trm[1]
  Coef[2] = 4*Trm[1] + Trm[2]
  Coef[3] = 4*Trm[1] + 2*Trm[2] + 2*Trm[3] + Trm[4]
  Coef[4] = -Trm[1]+Trm[3]
  ### Final formula for n^2 * Var of SDR var-estimator is (1/n^2)*
  # Coef[1]*mu4 + Coef[2]*mu3*mu + Coef[3]*(sig*mu)^2+Coef[4]*sig^4
  list(Coef=Coef, Trm = Trm) }

#' Function to calculate variance of SDR under stratified SRSWOR.
#' Written for clarity. If speed for larger n is an issue, use `VSDR.Gen`.
#' 
#' @param parmlst list of parameter vectors for each stratum. Required order is stratum size vector, sample size vector, stratum mean vector, stratum level standard deviation vector, and vector of stratum level domain percentages.
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @export
VSDR.exp <- function(parmlst, permarr) {
  ## NOTE parmlst has components in the order "Nh" "nh" "mean" "sig" "alpha"
  ## permarr is m x D, permarr[1,] constant, & each col a permutation
  ##    of same index-set, not necessarily 1:m
  muh = parmlst[[3]]
  sgsqh = parmlst[[4]]^2
  alphh = parmlst[[5]]
  nh = parmlst[[2]]
  Nh = parmlst[[1]]
  ph = Nh/sum(Nh)
  m = nrow(permarr)
  D = ncol(permarr)
  R = max(permarr[,1])
  n = sum(nh)
  allinds = 1:R
  useinds = sort(permarr[,1])
  altvec = numeric(R)
  altvec[useinds] = 1:m ## typically altvec = c(0,1:39,0,40  :78)
  gamvec = altvec[c(permarr[,1:D])]
  gamvec = rep(gamvec, (n+m*D)%/% (m*D))[1:(n+1)]
  if(n < m) gamvec[n+1]=gamvec[1]
  H = length(nh)
  Vtru = sum( (ph^2*alphh/nh)*(1-nh/Nh)*(sgsqh + (1-alphh)*muh^2) )
  sampvec = rep(1:H, nh) 
  EVSDR = sum( (ph^2/nh)*alphh*(muh^2+sgsqh))
  Mmat = Mcalc(n,m,D,permarr)
  MmatB = Mmat[ Mmat[,1] < Mmat[,2], ]
  SameHfac = rep(1,nrow(MmatB))
  ## indices in first 2 cols are 1:n, convert them to strata
  strmat = cbind(sampvec[MmatB[,1]],sampvec[MmatB[,2]])
  ind = which(strmat[,1]==strmat[,2])
  qh = muh*alphh*ph/nh
  SameHfac[ind] = SameHfac[ind] * ((Nh-1/alphh)/(Nh-1))[strmat[ind,1]]
  EVSDR = EVSDR + 2*sum(MmatB[,3]*qh[strmat[,1]]*
                          qh[strmat[,2]]*SameHfac)
  list(Vtru = Vtru, EVSDR=EVSDR, RelBias = EVSDR/Vtru - 1)  }


#' Function to calculate variance of SDR under stratified SRSWOR, for larger n.
#'
#' @description Preferred for speed. Similar to VSDRbias but is built to run faster for larger n.
#' 
#' @param parmlst list of parameter vectors for each stratum. Required order is stratum size vector, sample size vector, stratum mean vector, stratum level standard deviation vector, and vector of stratum level domain percentages.
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @export
VSDR.Gen <- function(parmlst, permarr) {
  ## NOTE parmlst has components in the order "Nh" "nh" "mean" "sig" "alpha"
  ## permarr is m x D, permarr[1,] constant, & each col a permutation
  muh = parmlst[[3]]
  sgsqh = parmlst[[4]]^2
  alphh = parmlst[[5]]
  nh = parmlst[[2]]
  Nh = parmlst[[1]]
  ph = Nh/sum(Nh)
  R = max(permarr[,1])
  m = nrow(permarr)
  D = ncol(permarr)
  n = sum(nh)
  allinds = 1:R
  useinds = sort(permarr[,1])
  altvec = numeric(R)
  altvec[useinds] = 1:m ## typically altvec = c(0,1:39,0,40  :78)
  gamvec = altvec[c(permarr[,1:D])]        ## corrected 7/11/26
  gamvec = rep(gamvec, (n+m*D)%/% (m*D))[1:(n+1)]
  if(n < m) gamvec[n+1]=gamvec[1]
  H = length(nh)
  Vtru = sum( (ph^2*alphh/nh)*(1-nh/Nh)*(sgsqh + (1-alphh)*muh^2) )
  sampvec = rep(1:H, nh) 
  mu.alp = muh*alphh*ph/nh
  EVSDR = sum( (ph^2/nh)*alphh*(muh^2+sgsqh))
  Sam.h = outer(mu.alp,mu.alp)
  diag(Sam.h) = alphh*(ph*muh/nh)^2*(Nh*alphh-1)/(Nh-1)
  
  baslst = split(1:(n+1), factor(gamvec))
  for(k in 1:min(m,n)) {        ## min(m,n) should = length(baslst)
    inds = baslst[[k]]
    indm = setdiff(inds, n+1)
    indg = setdiff(inds,1)
    if(length(indm)) {
      vecm = sampvec[indm] 
      mat1 = Sam.h[vecm,vecm]*outer(indm,indm, function(x,y) 
        (x!=y)*(1 + (gamvec[x+1]==gamvec[y+1]))) }  else mat1=0 
    if(length(indg)) {
      vecg = sampvec[indg-1]
      mat3 = Sam.h[vecg,vecg]*outer(indg,indg, function(x,y)
        gamvec[x-1]!=gamvec[y-1])  }  else mat3=0
    mat2 = if(min(length(indg),length(indm))) Sam.h[vecm,vecg] else 0
    EVSDR= EVSDR- sum(c(mat2))+ 0.5*(sum(c(mat1))+sum(c(mat3))) }
  list(Vtru = Vtru, EVSDR=EVSDR, RelBias = EVSDR/Vtru - 1)  }


#' Function to calculate correct expectation of the variance of SDR under stratified SRSWOR
#'
#' @description Calculates the expected variance of SDR under stratified SRSWOR for all combinations of m, D, and sample size, n.
#'
#' @param parmlst list of parameter vectors for each stratum. Required order is stratum size vector, sample size vector, stratum mean vector, stratum level standard deviation vector, and vector of stratum level domain percentages.
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @export
VSDRbias <- function(parmlst, permarr) {
  ## NOTE parmlst has components "Nh" "nh" "mean" "sig" "alpha"
  ## permarr is m x D, permarr[1,] constant, & each col a permutation
  muh = parmlst[[3]]
  sgsqh = parmlst[[4]]^2
  alphh = parmlst[[5]]
  nh = parmlst[[2]]
  Nh = parmlst[[1]]
  ph = Nh/sum(Nh)
  R = max(permarr[,1])
  m = nrow(permarr)
  D = ncol(permarr)
  n = sum(nh)
  allinds = 1:R
  useinds = sort(permarr[,1])
  altvec = numeric(R)
  altvec[useinds] = 1:m ## typically altvec = c(0,1:39,0,40:78)
  gamvec = altvec[c(permarr)]
  gamvec = rep(gamvec, (n+m*D)%/% (m*D))[1:(n+1)]
  if(n < m) gamvec[n+1]=gamvec[1]
  H = length(nh)
  sampvec = rep(1:H, nh) 
  
  tmp = (muh*alphh*ph/nh)[sampvec]
  SameStrat = outer(tmp,tmp)
  cumnh = cumsum(c(0,nh))
  for(h in 1:H) {
    indv = (cumnh[h]+1):cumnh[h+1]
    SameStrat[indv,indv] = alphh[h]*(ph[h]*muh[h]/nh[h])^2*
      (Nh[h]*alphh[h]-1)/(Nh[h]-1) }
  SameStrat[cbind(1:n,1:n)] = 0  
  Mat1 = outer(1:n, 1:n, function(x,y) 
    0.5*(gamvec[x]==gamvec[y] & x!=y)*(1+(gamvec[x+1]==gamvec[y+1])) -
      (gamvec[x]==gamvec[y+1]) + 0.5*(gamvec[x+1]==gamvec[y+1])*
      (gamvec[x]!=gamvec[y]))
  
  EVSDR = sum( (ph^2/nh)*alphh*(muh^2+sgsqh)) + 
    sum(c(Mat1*SameStrat))
  Vtru = sum( (ph^2*alphh/nh)*(1-nh/Nh)*(sgsqh + (1-alphh)*muh^2) )
  list(Vtru = Vtru, EVSDR=EVSDR, RelBias = EVSDR/Vtru-1) }

#' Computes SDR relbias when n less than or equal to m
#' 
#' @description This function computes the relative bias of SDR when n < m.
#'
#' @param parmlst list of parameter vectors for each stratum. Required order is stratum size vector, sample size vector, stratum mean vector, stratum level standard deviation vector, and vector of stratum level domain percentages.
#' @returns A list with relative bias, the true variance, and the original paramter list.
#' @export
RelSDRBias <- function(parmlst) {   ## the one to use when n <= m
       ## parmlst has components $Nh, $nh, $muh, $sigmah, $alphah
      ## where "alpha" in the formula is  p(A | N_h), h=1,...,H
    ## This function now modified to reflect new formula (18) in Prop.1 of paper
          Nh = parmlst[[1]]
          N = sum(Nh)
          nh = parmlst[[2]]
          ph = Nh/N
          muh = parmlst[[3]]
          sgsqh = parmlst[[4]]^2
          alphh = parmlst[[5]]
          Vtru = sum( (ph^2*alphh/nh)*(1-nh/Nh)*(sgsqh + (1-alphh)*muh^2) )
          alth = ph*alphh*muh/nh
          RelBias = ( sum(diff(c(alth,alth[1]))^2)/2 + 
                 sum(ph*alphh*(sgsqh + (1-alphh)*muh^2*
                     (1+((nh-1)*Nh/(nh^2*(Nh-1))))))/N )/Vtru
          list(RelBias=RelBias, Vtru=Vtru, parmlst=parmlst) }

#' A function that defines a list of G groups assciated with unique gamma-pairs
#' 
#' @description Illustrates how different D cycling choices, replicate weight-factor values and sample sizes can lead to different group sizes when cases are grouped by common factor values.
#' @param n Total Sample Size
#' @param m Number of replicates used
#' @param D D cycling parameter
#' @param permarr A matrix of indices from 1 to m with dimensions m x D, each column being a permutation of the first column.
#' @returns An ordered list of groups, ordered by the earliest index in the group from the original sample index order, with each list containing the indicies 1 to n that in that group. 
#' @examples
#' GpCalcT(200, D=3)
#' @export
GpCalcT <- function(n, m=78, D=1, permarr=PrmArr0[,1:D, drop=F]) {
     ## Splitting of observations 1:n into unique groups by
     ##   distinct pairs (gamma[j],gamma[j+1])
       gam = if(n<m) c(permarr[1:n,1],permarr[1,1]) else {
            K = (n+D*m-1) %/% (D*m)
            c(rep(permarr,K), permarr[1,1])[1:(n+1)]  }
     ## length(gam) is always  n+1
       gmat = cbind(gam[-(n+1)],gam[-1])
       idvals = c(gmat %*% runif(2))    ## hash coding for unique pairs
       gpinds = match(idvals, unique(idvals))
     ## output is a list of index-sets j corresponding to
     ## G unique (gam[j],gam[j+1]) values      
     ## Groups are now ordered by original order of sample obs
       split(1:n, gpinds) }


#' Generates random parameter list
#' 
#' @description A convenience function to generate somewhat random parmlsts for use in demonstration functions
#' @returns a parmlst with components ready to use in other functions
#' @examples
#' p <- ParmSim()
#' print(p)
#' @export
ParmSim <- function() {
       H = sample(5:15,1)
       list(Nh = 1e4*sample(5:25, H, replace=T),
         nh = sample(2:25, H, replace=T),
         muh = stats::rnorm(H),
         sigmah = stats::rgamma(H,3,3),
         alphah = runif(1, 0.5,0.8)*stats::rbeta(H,1,3)) }
