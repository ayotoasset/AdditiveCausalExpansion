#' Plot marginal effect.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param Xcol Either numeric or character, the latter only works if the training was done with a data.frame (i.e. not GPspline.train). If missing plots the 2D graph of Y on Z.
#' @param marginal A logical statement, default: \code{FALSE}. If \code{TRUE}, plots the marginal surface of the response surface Y with respect to the basis expanded model Z using average X.
#' @param Xlim A logical
#' @param Xstep A logical
#' @param Zlim A logical
#' @param Zstep A logical
#' @param truefun For simulation tests

#' @return Returns the plotly object


plot.GPspline <- function(object,Xcol,marginal=TRUE,Xlim = c(-1,1),Xstep = 0.1,Zlim = c(-1,1),Zstep = 0.1,truefun,contour=FALSE){
  #For the colors see http://jfly.iam.u-tokyo.ac.jp/color/#what

  Blue100    <- grDevices::rgb(0, 114, 178, max=255, alpha = (100-0)*255/100)
  SkyBlue100 <- grDevices::rgb(86, 180, 233, max=255, alpha = (100-0)*255/100)
  #cool = rainbow(10, start=rgb2hsv(col2rgb(SkyBlue100)), end=rgb2hsv(col2rgb(Blue100)),alpha=1)

  Orange100  <- grDevices::rgb(230, 159, 0, max=255, alpha = (100-0)*255/100)
  Vermillion50  <- grDevices::rgb(213, 94, 0, max=255, alpha = (100-50)*255/100)
  Vermillion100 <- grDevices::rgb(213, 94, 0, max=255, alpha = (100-0)*255/100)
  #warm = rainbow(10, start=rgb2hsv(col2rgb(Orange100)), end=rgb2hsv(col2rgb(Vermillion100)),alpha=1)

  #cols = c(rev(cool), rev(warm))
  #mypalette <- grDevices::colorRampPalette(cols)(255)

  px <- ncol(object$train_data$X)
  pz <- 1

  if(!missing(Zlim)) { Zlim = (Zlim - object$moments[1+px+pz,1])   / object$moments[1+px+pz,2]; }
  if(!missing(Zstep)){ Zstep = (Zstep - object$moments[1+px+pz,1]) / object$moments[1+px+pz,2]; }
  #for continuous Z so far!
  Zgrid <- seq(min(Zlim),max(Zlim),by=Zstep) * object$moments[1+px+pz,2] + object$moments[1+px+pz,1]

  #MIGHT NOT WORK FOR px>=2!!!
  medX <- apply(object$train_data$X * object$moments[1:(px+1),2] + object$moments[1:(px+1),1],2,median)

  if(!missing(Xcol)){
    #make numeric
    if(is.numeric(Xcol)) {
        idx <- (Xcol == seq(px))*1L
      }
    else if(is.character(Xcol) && (Xcol %in% colnames(object$train_data$X))) {
        idx <- ("abc"== c("abc","dcf","hef"))
        Xcol <- which(idx)
        idx <- idx*1L
      }
    else { stop("Supply either a numeric column index for X, or if  the object was generated with the data.frame wrapper, a valid variable name.\n") }

    #
    if(!missing(Xlim)) { Xlim = (Xlim - object$moments[1+Xcol,1])   / object$moments[1+Xcol,2]; }
    if(!missing(Xstep)){ Xstep = (Xstep - object$moments[1+Xcol,1]) / object$moments[1+Xcol,2]; }

    #Xval = apply(object$train_data$X,2,mean)
    n <- length(seq(min(Xlim),max(Xlim),by=Xstep))

    Xgrid <- seq(min(Xlim),max(Xlim),by=Xstep) * object$moments[1+Xcol,2] + object$moments[1+Xcol,1]
    grid <- base::expand.grid(X=Xgrid, Z=Zgrid)
    ngrid <- n^2

    #generate matrix with medians and replace column Xcol with the grid elements
    Xgrid2 <- t(matrix(rep(medX,n),px,ngrid))
    Xgrid2[,Xcol] <- grid$X

    cat("Predicting", n, "x",  n, "grid points\n")
    surface <- predict(object,Xgrid2,grid$Z,marginal=marginal)
    #truefun for simulation
    if(missing(truefun)){
      grid$Y <- surface$map
      gridL = grid; gridU = grid;
      gridL$Y <- surface$ci[,1]
      gridU$Y <- surface$ci[,2]
    } else {
      true <- truefun(Xgrid2,grid$Z)
      grid$Y <- surface$map - true
      gridL = grid; gridU = grid;
      gridL$Y <- surface$ci[,1] - true
      gridU$Y <- surface$ci[,2] - true
    }

    if(!requireNamespace("plotly", quietly = TRUE) && (contour==FALSE)){ #!plotly
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
      p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(grid$Y,n,n))
      p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(gridL$Y,n,n),opacity = 0.5)
      p <- plotly::add_surface(p,x = ~Xgrid, y = ~Zgrid,z = ~matrix(gridU$Y,n,n),opacity = 0.5)
      p
    } else {
      # Lattice ####
      cat("Plotly package not installed. Plotting with lattice.\n")
      #lattice::wireframe(Y ~ X*Z, data=grid,
      #                   shade=FALSE,drape = TRUE,pretty=TRUE)
      graphics::image(x=Xgrid,y=Zgrid,z=matrix(grid$Y,n,n),xlab=paste0("X[,",Xcol,"]"),ylab="Z")
      graphics::contour(x=Xgrid,y=Zgrid,z=matrix(grid$Y,n,n),add=TRUE,vfont=c("sans serif","bold"))
      obsX <- object$train_data$X[,Xcol] * object$moments[1+Xcol,2]+ object$moments[1+Xcol,1]
      obsZ <- object$train_data$Z * object$moments[1+px+pz,2] + object$moments[1+px+pz,1]
      points(obsX,obsZ,cex=0.5)

    }
  } else { #if not column of X selected
    n <- length(Zgrid)
    ngrid <- n

    Xgrid <- t(matrix(rep(medX,n),px,ngrid))

    cat("Predicting", n, "grid points\n")
    surface <- predict(object,Xgrid,Zgrid,marginal=marginal)
    Ygrid <- surface$map
    YgridL <- surface$ci[,1]
    YgridU <- surface$ci[,2]
    if(!marginal){
      obsY <- object$train_data$y * object$moments[1,2]       + object$moments[1,1]
      obsZ <- object$train_data$Z * object$moments[1+px+pz,2] + object$moments[1+px+pz,1]
      y_limit <- c(min(YgridL,obsY),max(YgridU,obsY))
    } else {
      obsY <- NULL
      obsZ <- NULL
      y_limit <- c(min(YgridL),max(YgridU))
    }
    plot(obsZ, obsY, xlim=range(Zgrid),cex=0.5,ylab="Outcome Y",xlab="Basis expanded Z",ylim = y_limit)
    if(!missing(truefun)){
      lines(Zgrid,truefun(Xgrid,Zgrid),lty=2,lwd=1.5)
    }
    lines(Zgrid,Ygrid,lty=2,lwd=2,col=Vermillion100)
    polygon(c(Zgrid, rev(Zgrid)),c(YgridL,rev(YgridU)),col = Vermillion50, border = FALSE)
  }
}
