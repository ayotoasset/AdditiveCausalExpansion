#' Plot marginal effect.
#'
#' @param object An "ace" object from the ace function
#' @param Xcol Either numeric or character, the latter only works if the training was done with a data.frame (i.e. not GPspline.train). If missing plots the 2D graph of Y on Z.
#' @param marginal A logical statement, default: \code{FALSE}. If \code{TRUE}, plots the marginal surface of the response surface Y with respect to the basis expanded model Z using average X.
#' @param plot3D A logical statemnt whether a three-dimensional plot should be produced. This is only valid if an X-column is selected and the plotly package is installed.
#' @param show.observations A logical
#' @param Xlim A two-dimensional vector
#' @param Xstep A scalar
#' @param Zlim A two-dimensional vector
#' @param Zstep A scalar
#' @param truefun A function handle for simulations.

#' @return Returns the plotly object


plot.ace <- function(object, Xcol, marginal = FALSE, plot3D = FALSE, show.observations = TRUE,
                     Xlim = c(-1,1), Xstep = 0.1, Zlim = c(-1,1), Zstep = 0.1, truefun) {

  # TODO: follow style guide
  # TODO: Doesnt look nice with discrete X variables

  graphics::par(mfrow=c(1, 1))
  #For the colors see http://jfly.iam.u-tokyo.ac.jp/color/#what
  Blue100         <- grDevices::rgb(0,   114, 178, max=255, alpha = (100 -  0)*255/100)
  Blue50          <- grDevices::rgb(0,   114, 178, max=255, alpha = (100 - 50)*255/100)
  BlueishGreen50  <- grDevices::rgb(0,   158, 115, max=255, alpha = (100 - 50)*255/100)
  BlueishGreen100 <- grDevices::rgb(0,   158, 115, max=255, alpha = (100 -  0)*255/100)
  SkyBlue100      <- grDevices::rgb(86,  180, 233, max=255, alpha = (100 -  0)*255/100)
  Orange100       <- grDevices::rgb(230, 159,   0, max=255, alpha = (100 -  0)*255/100)
  Vermillion50    <- grDevices::rgb(213,  94,   0, max=255, alpha = (100 - 50)*255/100)
  Vermillion100   <- grDevices::rgb(213,  94,   0, max=255, alpha = (100 -  0)*255/100)

  px <- ncol(object$train_data$X)
  pz <- 1

  if(!missing(Zlim)) { Zlim = (Zlim - object$moments[1+px+pz,1])   / object$moments[1+px+pz,2]; }
  if(!missing(Zstep)){ Zstep = (Zstep - object$moments[1+px+pz,1]) / object$moments[1+px+pz,2]; }

  medX= matrix(NaN,1,px)
  quantX= matrix(NaN,3,px)
  for(i in 1:px){
    quantX[,i] <- stats::quantile(object$train_data$X[,i] * object$moments[(i+1),2] + object$moments[(i+1),1],probs=c(0.25,0.5,0.75))
    medX[1,i] <- stats::median(object$train_data$X[,i] * object$moments[(i+1),2] + object$moments[(i+1),1])
  }

  if(!missing(Xcol)) {
    #stop("Supply either a numeric column index for X, or if the object was generated with the data.frame wrapper, a valid variable name.\n")
  #{ #}
    #make numeric
    if(is.numeric(Xcol)) {
      nXcol <- paste0("X[,",Xcol,"]")
      idx <- (Xcol == seq(1,px))*1L
    }
    else if(is.character(Xcol) && (Xcol %in% colnames(object$train_data$X))) {
      nXcol <- Xcol
      idx  <- (Xcol == colnames(object$train_data$X))
      Xcol <- which(idx)
      idx  <- as.integer(idx)
    }

    if(!missing(Xlim)) { Xlim = (Xlim - object$moments[1+Xcol,1])   / object$moments[1+Xcol,2]; }
    if(!missing(Xstep)){ Xstep = (Xstep - object$moments[1+Xcol,1]) / object$moments[1+Xcol,2]; }
    nx <- length(seq(min(Xlim),max(Xlim),by=Xstep))

    Xgrid <- seq(min(Xlim), max(Xlim), by=Xstep) * object$moments[1 + Xcol, 2] + object$moments[1 + Xcol, 1]
  }

  if(show.observations){
    if(!missing(Xcol)){
      obsX <- object$train_data$X[,Xcol] * object$moments[1 + Xcol, 2]+ object$moments[1 + Xcol, 1]
    }
    obsZ <- object$train_data$Z * object$moments[1 + px + pz, 2]
            + object$moments[1 + px + pz, 1]
    obsY <- object$train_data$y * object$moments[1, 2]
            + object$moments[1, 1]
  }

  if (!object$train_data$Zbinary) {
    if(marginal && show.observations){
      cat("Observations for marginal predictions use predicted value as the true marginal is not observed\n")
    }
    if (missing(truefun)) {
      if (marginal==TRUE) {
        y.label <- "Marginal Y"
      } else {
        y.label <- "Outcome Y"
      }
    } else {
      if (show.observations) {
        cat("True values subtracted from Observations\n")
      }
      if (marginal==TRUE) {
        y.label <- "Diff. to true marg. Y"
      } else {
        y.label <- "Diff. to true Y"
      }
    }

    #for continuous Z so far!
    Zgrid <- seq(min(Zlim),max(Zlim),by=Zstep) * object$moments[1+px+pz,2] + object$moments[1+px+pz,1]
    nz <- length(seq(min(Zlim),max(Zlim),by=Zstep))

    if(!missing(Xcol)){

      grid <- base::expand.grid(X=Xgrid, Z=Zgrid)
      ngrid <- nx*nz

      #generate matrix with medians and replace column Xcol with the grid elements
      Xgrid2 <- t(matrix(rep(medX,nz),px,ngrid))
      Xgrid2[,Xcol] <- grid$X

      cat("Predicting", nz, "x",  nx, "grid points\n")

      surface <- predict(object,newX = as.matrix(Xgrid2),newZ = grid$Z,marginal=marginal)

      if(marginal && show.observations){
        invisible(utils::capture.output(obsY <- predict(object ,marginal = TRUE)$map))
      }

      #truefun for simulation
      if (missing(truefun)) {
        grid$Y  <- surface$map
        gridL   <- grid
        gridU   <- grid
        gridL$Y <- surface$ci[, 1]
        gridU$Y <- surface$ci[, 2]
      } else {
        #Xgrid2_notnormalized = Xgrid2
        #only single dimensional X right now:
        true <- truefun(grid$X,# * object$moments[1+Xcol,2] + object$moments[1+Xcol,1],
                        grid$Z)# * object$moments[1+px+pz,2] + object$moments[1+px+pz,1])
        gridT   <- grid
        gridT$Y <- true
        grid$Y  <- surface$map #- true
        gridL   <- grid
        gridU   <- grid
        gridL$Y <- surface$ci[,1] #- true
        gridU$Y <- surface$ci[,2] #- true

        if(show.observations==TRUE){
          #Xgrid3 <- t(matrix(rep(medX,length(obsZ)),px,length(obsZ)))
          #Xgrid3[,Xcol] <- obsX
          obsY <- obsY #- predict(object,newX = Xgrid3,newZ = grid$Z,marginal=marginal)
        }

      }

      if(requireNamespace("plotly", quietly = TRUE) && (plot3D == TRUE)) { #!plotly
        ## plotly ####
        Xgrid <- matrix(Xgrid,nx,nz)
        Zgrid <- t(matrix(Zgrid,nz,nx))

        cat("Plotting with plotly. This can take several seconds\n")
        p <- plotly::plot_ly(showscale = FALSE, type = "surface")
        if(show.observations){
          p <- plotly::add_trace(p,x = obsX,  y = c(obsZ), z = c(obsY), mode = "markers", type = "scatter3d",
                               marker = list(size = 5, color = "red", symbol = 104),inherit=FALSE,name="Observations")
        }
        if(!missing(truefun)){
          p <- plotly::add_trace(p,type = "scatter3d", mode = "markers",x = c(Xgrid), y = c(Zgrid),z = c(matrix(gridT$Y,nx,nz)),
                                 opacity = 0.8,name="True",marker = list(size = 2, color = "black", symbol = 104),
                                 #colors = c("#000000"),
                                 inherit=FALSE)
        }

        p <- plotly::layout(p,
            title = paste0("Marginal effect of Z conditional on ",nXcol," and median covariates"),
            scene = list(
              camera = list(eye = list(x = -1.25, y = -1.25, z = 1.25)),
              xaxis = list(title = paste0("Control ",nXcol)),
              yaxis = list(title = "Treatment Z"),
              zaxis = list(title = y.label)
            ))
        p <- plotly::add_surface(p, x = ~Xgrid, y = ~Zgrid, z = ~matrix(grid$Y,nx,nz),name="MAP")
        p <- plotly::add_surface(p, x = ~Xgrid, y = ~Zgrid, z = ~matrix(gridL$Y,nx,nz),opacity = 0.5,name="Lower CI")
        p <- plotly::add_surface(p, x = ~Xgrid, y = ~Zgrid, z = ~matrix(gridU$Y,nx,nz),opacity = 0.5,name="Upper CI")

        p

      } else if (requireNamespace("ggplot2", quietly = TRUE)) {
        if (plot3D == TRUE) {
          cat("plotly package not found. Proceeding with contour plots")
        }
        cat("Plotting with ggplot2\n")

        g <- ggplot2::ggplot(data.frame(x=grid$X,y=grid$Z,z=grid$Y), ggplot2::aes(x=x,y=y)) +
          ggplot2::geom_tile(ggplot2::aes(fill=z)) +
          ggplot2::scale_fill_gradient2(low=Blue100, high=Vermillion100,name=y.label)
        if (show.observations){
          g <- g + ggplot2::geom_point(data = data.frame(xobs = obsX, yobs = c(obsZ)),
                                       ggplot2::aes(x = "xobs", y = "yobs"),
                                       pch = 21, inherit.aes = FALSE, show.legend = TRUE)
        }

        g <- g + ggplot2::scale_x_continuous(name=paste0("Control ",nXcol),expand=c(0,0)) +
          ggplot2::scale_y_continuous(name="Basis Expanded Z",expand=c(0,0)) +
          ggplot2::theme_minimal()
        if( (sign(min(grid$Y))!=sign(max(grid$Y))) && (marginal==TRUE) ){
          #plot 0-line if the range of the marginal output is zero-crossing
          g <- g + ggplot2::geom_contour(bins=1,col="black",ggplot2::aes(z=z),lty=2)
        }
        g

      } else {
        # filled.contour ####
        if(plot3D==TRUE){cat("plotly package not found. Proceeding with contour plots")}
        if(requireNamespace("plotly", quietly = TRUE)){ cat("plotly package not installed\n")}
        if(requireNamespace("ggplot2", quietly = TRUE)){ cat("ggplot2 package not installed\n")}
        cat("Plotting with the standard R package graphics\n")

        #par(mar=c(4,4,1,0))
        nlevels_quadrant=10

        n_levels_pos = floor(nlevels_quadrant*max(grid$Y)/(max(grid$Y)+abs(min(grid$Y))))
        n_levels_neg = floor(nlevels_quadrant*abs(min(grid$Y))/(max(grid$Y)+abs(min(grid$Y))))

        filled.contour.lvls <- sort(unique(c(seq(min(grid$Y), 0, length.out = n_levels_neg),0,
                                 seq(0, max(grid$Y), length.out = n_levels_pos))),decreasing = FALSE)

        filled.contour.cols <- grDevices::colorRampPalette(c("white","darkgray"))(length(filled.contour.lvls) - 1)

        graphics::filled.contour(x=Xgrid,y=Zgrid,z=matrix(grid$Y,nz,nx),axes=TRUE,
                                 xlab=paste0("Control ",nXcol),ylab="Z",
                                 col=filled.contour.cols, levels=filled.contour.lvls,
                                 plot.axes = { graphics::axis(1)
                                               graphics::axis(2)
                                               if(show.observations){
                                                 graphics::points(obsX, obsZ,cex=0.5)
                                                 }
                                               }
                                 )
      }
    } else { #if not column of X selected
      cat("No column of X selected. Proceeding with 2D plot with respect to Z\n")
      nz <- length(Zgrid)
      ngrid <- nz

      cat("Predicting", nz, "grid points\n")
      Xgrid <- t(matrix(rep(medX,nz),px,ngrid))
      surface <- predict(object,Xgrid,Zgrid,marginal=marginal)

      Ygrid <- surface$map
      YgridL <- surface$ci[,1]
      YgridU <- surface$ci[,2]
      if(!marginal && show.observations){
        #obsY <- object$train_data$y * object$moments[1,2]       + object$moments[1,1]
        #obsZ <- object$train_data$Z * object$moments[1+px+pz,2] + object$moments[1+px+pz,1]
        y_limit <- c(min(YgridL,obsY),max(YgridU,obsY))
      } else {
        obsY <- NULL
        obsZ <- NULL
        y_limit <- c(min(YgridL),max(YgridU))
      }
      graphics::plot(obsZ, obsY, xlim=range(Zgrid),cex=0.5,ylab=y.label,xlab="Basis expanded Z",ylim = y_limit)
      graphics::grid()
      if(!missing(truefun)){
        graphics::lines(Zgrid,truefun(Xgrid,Zgrid),lty=2,lwd=1.5)
      }
      graphics::lines(Zgrid,Ygrid,lty=2,lwd=2,col=Vermillion100)
      graphics::polygon(c(Zgrid, rev(Zgrid)),c(YgridL,rev(YgridU)),col = Vermillion50, border = FALSE)
    }
  } else {
    #Binary Z
    if(plot3D){cat("No 3D plot available for binary Z\n")}

    if(missing(truefun)){
      if(marginal==TRUE){
        y.label <- "Y(x,1) - Y(x,0)"
      } else {
        y.label <- "Outcome Y"
      }
    } else {
      if(marginal==TRUE){
        y.label <- "Y(x,1) - Y(x,0)"
      } else {
        y.label <- "Outcome Y"
      }
    }

    if(!missing(Xcol)){
      cat("Proceeding to plot binary Z\n")

      obsZ          <- object$train_data$Z
      Xgrid2        <- t(matrix(rep(medX,nx),px,nx))
      Xgrid2[,Xcol] <- Xgrid
      n <- nrow(Xgrid2)

      if(marginal){
        #one curve ("treatment effect")
        cat("Marginal: Predicting", n, "grid points\n")
        surface <- predict.ace(object, as.matrix(Xgrid2), 0, marginal = TRUE) #Z does not matter

        y_limit <- c(min(surface$ci[,1]),max(surface$ci[,2]))

        graphics::plot(Xgrid,surface$map,col=BlueishGreen100,lwd=2,ylim=y_limit,type="l",xlab=paste0("Control ",nXcol),ylab=y.label)
        graphics::polygon(c(Xgrid,rev(Xgrid)),c(surface$ci[,1],rev(surface$ci[,2])),col=BlueishGreen50,density=20)
        if(!missing(truefun)){
          graphics::lines(Xgrid, truefun(matrix(Xgrid2)),lty=2)
        }
        if(show.observations){
          graphics::points(obsX[obsZ==0],rep(y_limit[1],sum(obsZ==0)),pch=21,cex=0.6)
          graphics::points(obsX[obsZ==1],rep(y_limit[2],sum(obsZ==1)),pch=20,cex=0.6)
        }
      } else {
        #two curves
        surface0 <- predict(object,as.matrix(Xgrid2),0,marginal=FALSE)
        surface1 <- predict(object,as.matrix(Xgrid2),1,marginal=FALSE)

        if(show.observations){
          y_limit <- c(min(surface0$ci[,1],surface1$ci[,1],obsY),
                       max(surface0$ci[,2],surface1$ci[,2],obsY))
        } else {
          y_limit <- c(min(surface0$ci[,1],surface1$ci[,1]),
                       max(surface0$ci[,2],surface1$ci[,2]))
        }

        graphics::plot(Xgrid, surface0$map, col=Blue100, lwd=2, ylim=y_limit,
                       type="l", xlab=paste0("Control ", nXcol), ylab=paste0(y.label, " (blue: Z=1,red: Z=0)"))
        graphics::polygon(c(Xgrid, rev(Xgrid)), c(surface0$ci[, 1], rev(surface0$ci[, 2])), col=Blue50, density=20, angle = -45)
        graphics::lines(Xgrid, surface1$map,col=Vermillion100, lwd=2)
        graphics::polygon(c(Xgrid, rev(Xgrid)), c(surface1$ci[, 1], rev(surface1$ci[, 2])), col=Vermillion50, density=20)
        if(!missing(truefun)){
          graphics::lines(Xgrid,truefun((Xgrid2),rep(0,nx)),lty=2)
          graphics::lines(Xgrid,truefun((Xgrid2),rep(1,nx)),lty=2)
        }
        if(show.observations){
          graphics::points(obsX[obsZ==0],rep(y_limit[1],sum(obsZ==0)),pch=21,cex=0.6)
          graphics::points(obsX[obsZ==1],rep(y_limit[2],sum(obsZ==1)),pch=20,cex=0.6)
        }
      }
    } else {
      cat("No plot produced, need to specify an X-column for binary Z")
      FALSE
    }

  }
}
