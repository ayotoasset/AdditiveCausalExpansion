norm_variables <- function(y,X,moments){
  px = ncol(X); pz = ncol(Z);

  #use function for both in and out-of-sample normalization
  if(missing(moments)){
    moments = list(meanX = rep(0,px),
                   varX  = rep(1,px),
                   meanZ = rep(0,pz),
                   varZ  = rep(1,pz),
                   meanY  = 0*mean(y), # we do not center the outcome surface
                   varY  = stats::var(y))

    #only normalize non-binary variables in X
    nonbinary_varbs = unique(1:px * (apply(X,2,function(x) { all(stats::na.omit(x) %in% 0:1) })==0))
    nonbinary_varbs = nonbinary_varbs[nonbinary_varbs!=0]
    within_unit_circle = unique(1:px * apply(X,2,function(x) { all(stats::na.omit(x)<=1 & stats::na.omit(x)>=0 ) }) )
    within_unit_circle = within_unit_circle[within_unit_circle!=0]

    normalize_varbs = setdiff(nonbinary_varbs,within_unit_circle)

    tmp = data.frame(X[,normalize_varbs])

    moments$meanX[normalize_varbs] = apply(tmp,2,mean)
    moments$varX[normalize_varbs]  = apply(abs(tmp),2,max)^2

    #only normalize non-binary variables in Z
    nonbinary_varbs = unique(1:pz * (apply(Z,2,function(x) { all(stats::na.omit(x) %in% 0:1) })==0))
    nonbinary_varbs = nonbinary_varbs[nonbinary_varbs!=0]
    within_unit_circle = unique(1:pz * apply(Z,2,function(x) { all(stats::na.omit(x)<=1 & stats::na.omit(x)>=0 ) }) )
    within_unit_circle = within_unit_circle[within_unit_circle!=0]

    normalize_varbs = setdiff(nonbinary_varbs,within_unit_circle)

    tmp = data.frame(X[,normalize_varbs])

    moments$meanZ[normalize_varbs] = apply(tmp,2,mean)
    moments$varZ[normalize_varbs]  = apply(abs(tmp),2,max)^2

  }
  mynorm <- function(i){ (X[,i]-moments$meanX[i]) / sqrt(moments$varX[i]) }
  X = data.frame(sapply(1:p,mynorm))

  mynorm <- function(i){ (Z[,i]-moments$meanZ[i]) / sqrt(moments$varZ[i]) }
  Z = data.frame(sapply(1:p,mynorm))

  if(missing(y)){
    out = list(X=X,Z=Z)
  } else {
    y = (y - moments$meanY ) / sqrt(moments$varY)
    out = list(moments=moments, y=y, X=X, Z=Z)
  }
  out
}


check_inputs <- function(y,X,Z){
  if(length(y) != nrow(X) ){
    warning("y and X: nr of observations mismatch", call. = FALSE)
  }
  if(length(y) != nrow(Z) ){
    warning("y and z: nr of observations mismatch", call. = FALSE)
  }
  if(class(X) != "data.frame"){
    warning("X is not of class data.frame", call. = FALSE)
  }
}
