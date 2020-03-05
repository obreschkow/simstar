#' Visualise one or seveal 3D point sets
#' 
#' @importFrom cooltools nplot rasterflip lim griddata2 kde2 quadrupole rotation3
#' @importFrom png writePNG
#' @importFrom EBImage gblur
#' @importFrom grDevices pdf dev.off col2rgb
#' @importFrom graphics axis lines par rasterImage rect text
#'
#' @description Produces a square raster image visualising particle positions in 3D N-body/SPH data.
#'
#' @param x list of n-by-3 matrices, represending the coordinates of n 3D positions; n can be different for each matrix. The particles stored in different matrices represent different types of particles and are visualised in different colors.
#' @param col vector of colors used for the different particle types, corresponding to the matrix-valued elements in the list \code{x}.
#' @param center 3-vector specifying the coordinate at the center of the plot. The default is the center of mass.
#' @param radius radius (= half horizontal/vertical diameter) of the image in the same length units as \code{x}. The default is the largest distance between the center of mass and all points, multiplied by \code{radius.scale}.
#' @param radius.scale scaling factor of the \code{radius}, if the latter is not specified by the user.
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pngfile optional png-filename to save the image as raster image.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param rotation either an integer (1-6) or a 3-vector specifying a rotationn of the 3D particle positions. For integer, 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,y)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specified the rotation axis and its norm the rotation angle in the positive geometric sense.
#' @param kde logical flag, specifying whether the particles are smoothed using an addaptive kernel density estimator.
#' @param ngrid number of grid cells per side in the output image.
#' @param title Text to be added to the figure.
#' @param lum overall luminosity scaling factor (default 1).
#' @param shadows differential luminosity scaling factor for darker regions (default 1).
#' @param sigma relative smoothing scale (default 1). 
#' @param arrows logical flag, specifying if axis-arrows are drawn.
#' @param side vector made of up to four integer values (1/2/3/4), specifying which axes with tick marks are to be drawn (see \code{\link{axis}}).
#' @param xlab label on horizontal arrow (only shown if \code{arrows=TRUE}).
#' @param ylab label on vertical arrow (only shown if \code{arrows=TRUE}).
#' @param add logical flag specifying whether the images added to an existing plot. If TRUE, no pdf-file can be produced with this routine.
#' @param add.xlim 2-vector specifing the horiziontal position (left and right) of the image on the existing plot if \code{add=TRUE}.
#' @param add.ylim 2-vector specifing the vertical position (bottom and top) of the image on the existing plot if \code{add=TRUE}.
#'
#' @seealso \code{\link{sphview4}}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Example of 2x1e5 particles representing a homogenous sphere contained in another sphere
#' x.blue = cooltools::runif3(1e5)
#' x.red = cooltools::runif3(1e5, r = c(0,0.5))
#' sphview(list(x.blue,x.red), radius.scale=1.2)
#'
#' @export

sphview = function(x, col = c('#0515ff', '#ff0010', 'green', 'orange', 'yellow', 'purple'),
                   center = NULL, radius = NULL, radius.scale = 1, screen = TRUE, pngfile = NULL, pdffile = NULL,
                   rotation = 1, kde = TRUE, ngrid = 300, title = NULL, lum = 1, shadows = 1, sigma = 1,
                   arrows = TRUE, side = NULL, xlab = NULL, ylab = NULL,
                   add = FALSE, add.xlim = NULL, add.ylim = NULL) {

  # handle x
  if (is.list(x)) {
    n = length(x)
    if (n>6) stop('The list x must not contain more than 5 arrays.')
  } else {
    x = list(x)
    n = 1
  }

  # concatenate list
  xc = x[[1]]
  if (n>1) {
    for (i in seq(2,n)) xc = rbind(xc,x[[i]])
  }
  x0 = apply(xc,2,mean) # geometric centre
  npoints = dim(xc)[1]

  # plot parameters
  d = 0.05  # distance to margin for arrows ant title
  larrow = 0.1   # length of arrows
  mar = c(0,0,0,0)
  col = t(col2rgb(col)/255)

  # rotation matrix
  if (length(rotation)==3) {
    rot = t(cooltools::rotation3(rotation))
    if (is.null(xlab)) xlab = expression(e[1])
    if (is.null(ylab)) ylab = expression(e[2])
  } else {
    if (length(rotation)!=1) stop('rotation must be an integer 1,...,6 or a real 3-vector.')
    if (rotation>3) e = eigen(cooltools::quadrupole(xc))$vectors
    if (rotation==1) {
      rot = diag(3)
      if (is.null(xlab)) xlab = 'x'
      if (is.null(ylab)) ylab = 'y'
    } else if (rotation==2) {
      rot = rbind(c(0,0,1),c(1,0,0),c(0,1,0))
      if (is.null(xlab)) xlab = 'y'
      if (is.null(ylab)) ylab = 'z'
    } else if (rotation==3) {
      rot = rbind(c(0,1,0),c(0,0,1),c(1,0,0))
      if (is.null(xlab)) xlab = 'z'
      if (is.null(ylab)) ylab = 'x'
    } else if (rotation==4) {
      rot = e[,c(1,3,2)]
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression(lambda ['min'])
    } else if (rotation==5) {
      rot = e[,c(1,2,3)]
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression(lambda ['mid'])
    } else if (rotation==6) {
      rot = e[,c(2,3,1)]
      if (is.null(xlab)) xlab = expression(lambda ['mid'])
      if (is.null(ylab)) ylab = expression(lambda ['min'])
    } else {
      stop('rotation must be an integer 1,...,6 or a real 3-vector.')
    }
  }

  # center points
  if (is.null(center)) {
    center = x0
  } else {
    if (length(center)==2) center=c(center,0)
  }
  for (i in seq(n)) {
    x[[i]] = t(t(x[[i]])-center)
  }

  # determine radius
  if (is.null(radius)) {
    radius = sqrt(max(apply(t(t(xc)-x0)^2,1,sum)))*radius.scale
  }

  # rotate and project coordinates
  for (i in seq(n)) {
    x[[i]] = (x[[i]]%*%rot)[,1:2]
  }

  # raster data
  smoothing = 0.01*sigma*ngrid
  img = array(0,c(ngrid,ngrid,3))

  for (i in seq(n)) {

    if (smoothing==0) {
      g = cooltools::griddata2(x[[i]][,1],x[[i]][,2],xlim=c(-1,1)*radius,ylim=c(-1,1)*radius,n=ngrid)
    } else {
      if (kde) {
        g = cooltools::kde2(x[[i]][,1],x[[i]][,2],xlim=c(-1,1)*radius,ylim=c(-1,1)*radius,n=ngrid,s = 0.2*smoothing^0.5, sd.max=smoothing*2)
        g$n = g$d
      } else {
        g = cooltools::griddata2(x[[i]][,1],x[[i]][,2],xlim=c(-1,1)*radius,ylim=c(-1,1)*radius,n=ngrid)
        g$n = EBImage::gblur(g$n,smoothing)
      }
    }

    for (k in seq(3)) {
      img[,,k] = img[,,k]+col[i,k]*g$n
    }

  }

  # finalize image
  img = atan(img/npoints*ngrid^2*lum)/pi*2
  f = 10^max(0,shadows)*2
  img = log10(f*img+1)/log10(f+1)
  img = cooltools::lim(img) # just to be sure

  # save raster image as png
  if (!is.null(pngfile)) {
    png::writePNG(cooltools::rasterflip(img),pngfile)
  }

  # determine margin
  mar = rep(0.25,4)
  for (i in seq(4)) {
    if (i%in%side) mar[i]=3
  }

  # show on screen and save as pdf
  for (mode in seq(2)) {

    make = FALSE

    if (mode==1 & screen) {
      make = TRUE
    }

    if (mode==2 & !is.null(pdffile) & !add) {
      make = TRUE
      grDevices::pdf(pdffile,width=7,height=7)
    }

    if (make) {

      # initialize plot
      if (!add) cooltools::nplot(xlim=c(0,1), ylim=c(0,1), pty='s', mar=mar)

      # plot raster
      if (add) {
        xleft = add.xlim[1]
        ybottom = add.ylim[1]
        dx = add.xlim[2]-add.xlim[1]
        dy = add.ylim[2]-add.ylim[1]
      } else {
        xleft = 0
        ybottom = 0
        dx = 1
        dy = 1
      }
      rasterImage(cooltools::rasterflip(img),xleft,ybottom,xleft+dx,ybottom+dy)

      # arrows
      if (arrows) {
        if (!is.null(xlab)) {
          text(xleft+d+larrow,ybottom+d,xlab,pos=4,col='white')
          arrows(xleft+d,ybottom+d,xleft+d+larrow,ybottom+d,col='white',length = 0.1,angle=20)
        }
        if (!is.null(ylab)) {
          text(xleft+d*1.2,ybottom+1.7*d+larrow,ylab,pos=4,col='white', offset=-0.4)
          arrows(xleft+d,ybottom+d,xleft+d,ybottom+d+larrow,col='white',length = 0.1,angle=20)
        }
      }
      
      # title
      if (!is.null(title)) text(d*1.2, 1-d*1.2, title, pos=4, col='white', offset=-0.4)

      # side axes
      par(xpd=T)
      rect(xleft,ybottom,xleft+dx,ybottom+dy)
      par(xpd=F)
      usr = par()$usr
      par(usr=c(-1,1,-1,1)*radius)
      for (i in seq(4)) {
        if (i%in%side) axis(side=i)
      }
      par(usr=usr)

      if (mode==2) {dev.off()}

    }

  }

}
