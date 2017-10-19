#' Plot marginal effect.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param marginal A matrix with new data. If not presented, using the training data.
#' @param plotly A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param Xlim A logical
#' @param Xstep A logical
#' @param Zlim A logical
#' @param Zstep A logical

#' @return Returns the plotly object


plot.GPspline <- function(object,marginal=TRUE,plotly=TRUE,Xlim = c(-1,1),Xstep = 0.1,Zlim = c(-1,1),Zstep = 0.1){

  #if (requireNamespace("rgl", quietly = TRUE)) {
  #  rgl::plot3d(...)
  #} else {
  #  ## do something else not involving rgl.
  #}

    #one-dimensionalX
    if(!missing(Xlim)) { Xlim = (Xlim - object$moments[2,1])   / object$moments[2,2]; }
    if(!missing(Xstep)){ Xstep = (Xstep - object$moments[2,1]) / object$moments[2,2]; }
    if(!missing(Zlim)) { Zlim = (Zlim - object$moments[3,1])   / object$moments[3,2]; }
    if(!missing(Zstep)){ Zstep = (Zstep - object$moments[3,1]) / object$moments[3,2]; }


    #Xval = apply(object$train_data$X,2,mean)
    n <- length(seq(min(Xlim),max(Xlim),by=Xstep))
    Xgrid <- seq(min(Xlim),max(Xlim),by=Xstep) * object$moments[2,2] + object$moments[2,1]
    Zgrid <- seq(min(Zlim),max(Zlim),by=Zstep) * object$moments[3,2] + object$moments[3,1]

    grid <- base::expand.grid(X=Xgrid, Z=Zgrid)
    Xgrid2 <- grid$X
    Zgrid2 <- grid$Z
    cat("Predicting", n, "x",  n, "grid points\n")
    marginal_surface <- predict(object,Xgrid2,Zgrid2,marginal=marginal)

    grid$Y <- marginal_surface$map
    gridL = grid; gridU = grid;
    gridL$Y <- marginal_surface$ci[,1]
    gridU$Y <- marginal_surface$ci[,2]

  if(!plotly){
    # Lattice ####
    lattice::wireframe(Y ~ X*Z, data=grid,
                       shade=FALSE,drape = TRUE,pretty=TRUE)
  } else {
    ## plotly ####
    Xgrid <- matrix(Xgrid,n,n)
    Zgrid <- t(matrix(Zgrid,n,n))

    cat("Plotting with plotly. This can take several seconds.\n")
    p <- plotly::plot_ly(showscale = FALSE)
    p <- plotly::layout(p,
        title = "Marginal effect of Z conditional on X[abc] and average covariates",
        scene = list(
          camera = list(eye = list(x = -1.25, y = -1.25, z = 1.25)),
          xaxis = list(title = "Control X"),
          yaxis = list(title = "Treatment Z"),
          zaxis = list(title = "Outcome Y")
        ))
    p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(marginal_surface$map,n,n))
    p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(marginal_surface$ci[,1],n,n),opacity = 0.5)
    p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(marginal_surface$ci[,2],n,n),opacity = 0.5)
    p
  }

  }
