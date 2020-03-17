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
#' @param x n-by-3 matrix, representing the 3D-coordinates of n particles.
#' @param species optional n-element vector giving the integer species indices, which will specify the color of the particles. If not given, all particles will have the first colour given in the vector \code{col}.
#' @param col vector of colors, specifying the color of the species given in the optional species vector.
#' @param center 3-vector specifying the coordinate at the center of the plot. The default is the center of mass.
#' @param xlim 2-vector specifying the horizontal range of the image in the units of \code{x}. The default is c(-1,1) times the largest distance between the center of mass and all points, multiplied by \code{auto.scale}.
#' @param ylim 2-vector specifying the vertical range of the image in the units of \code{x}. The default is c(-1,1) times the largest distance between the center of mass and all points, multiplied by \code{auto.scale}.
#' @param auto.scale number of 2-vector specifying the scaling factor for the automatically determined xlim and ylim, if the latter are not specified by the user.
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pngfile optional png-filename to save the image as raster image.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param rotation either an integer (1-6) or a 3-vector specifying a rotationn of the 3D particle positions. In case of an integer: 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,y)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specifies the rotation axis and its norm the rotation angle in radians in the positive geometric sense.
#' @param kde logical flag, specifying whether the particles are smoothed using an addaptive kernel density estimator.
#' @param ngrid number of grid cells per side in the output image. If the image is not square ngrid is interpreted as the geometric mean between the horizontal and the vertical number of pixels, such that the total number of pixels remains about ngrid^2.
#' @param title Text to be added to the figure.
#' @param title.origin optional 2-vector specifying the position of the title
#' @param lum overall luminosity scaling factor (default 1).
#' @param shadows differential luminosity scaling factor for darker regions (default 1).
#' @param sigma relative smoothing scale (default 1). 
#' @param arrows logical flag, specifying if axis-arrows are drawn.
#' @param arrow.origin optional 2-vector specifying the origin of the arrows
#' @param arrow.length optional number specitying the length of the arrows
#' @param arrow.lwd line width of arrows
#' @param side vector made of up to four integer values (1/2/3/4), specifying which axes with tick marks are to be drawn (see \code{\link{axis}}).
#' @param xlab label on horizontal arrow (only shown if \code{arrows=TRUE}).
#' @param ylab label on vertical arrow (only shown if \code{arrows=TRUE}).
#' @param cex text size
#' 
#' @return Returns a list of items
#' \item{rgb}{ngrid-by-ngrid-by-3 array of containig the rgb colour data between 0 and 1}
#' \item{radius}{positive real number specifying the radius of the image along the horizontal and vertical axes, in the units of the particle positions}
#' \item{rotationmatrix}{rotation matrix applied to the input data, i.e. the projection was x = (x\%*\%rot)[,1:2]}
#'
#' @seealso \code{\link{sphview4}}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Example of 2x1e5 particles representing a homogenous sphere contained in another sphere
#' x.red = cooltools::runif3(1e5, r = c(0,0.5))
#' x.blue = cooltools::runif3(1e5)
#' sphview(rbind(x.red,x.blue), c(rep(1,1e5),rep(2,1e5)), auto.scale=1.2)
#'
#' @export

sphview = function(x, species, col = c('#ff0010', '#0515ff', 'green', 'orange', 'yellow', 'purple'),
                   center = NULL, xlim = NULL, ylim = NULL, auto.scale = 1, screen = TRUE, pngfile = NULL, pdffile = NULL,
                   rotation = 1, kde = TRUE, ngrid = 300, title = NULL, title.origin = NULL, lum = 1, shadows = 1, sigma = 1,
                   arrows = TRUE, arrow.origin = NULL, arrow.length = NULL, arrow.lwd = 1.5,
                   side = NULL, xlab = NULL, ylab = NULL, cex=1) {
  
  # rescale luminosity scale (values matched to typical sph halos)
  lum = 0.2*lum
  shadows = 1.5*shadows

  # handle x
  if (is.array(x)) {
    if (length(dim(x))==2) {
      n = dim(x)[1]
      if (dim(x)[2]==2) x = cbind(x,rep(0,n))
      if (dim(x)[2]!=3) stop('x must be an n-by-3 dimensional array')
    } else {
      stop('x must be an n-by-3 dimensional array')
    }
  } else {
    stop('x must be an n-by-3 dimensional array')
  }
  
  # handle species
  if (is.null(species)) {
    species = rep(1,n)
  } else {
    if (length(species)!=n) stop('species must be a vector of same length as number of rows in x.')
    if (min(species)<1) stop('species must contain only positive integers.')
    if (max(species)>length(col)) stop('not enough colours given to match the maximum species value.')
  }

  # determine geometrix centre
  x0 = apply(x,2,mean) # geometric centre

  # plot parameters
  mar = c(0,0,0,0)
  col = t(col2rgb(col)/255)

  # rotation matrix
  if (length(rotation)==3) {
    rot = t(cooltools::rotation3(rotation))
    if (is.null(xlab)) xlab = expression(e[1])
    if (is.null(ylab)) ylab = expression(e[2])
  } else {
    if (length(rotation)!=1) stop('rotation must be an integer 1,...,6 or a real 3-vector.')
    if (rotation>3) e = eigen(cooltools::quadrupole(x))$vectors
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
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
    } else if (rotation==5) {
      rot = e[,c(1,2,3)]
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['mid'])
    } else if (rotation==6) {
      rot = e[,c(2,3,1)]
      if (is.null(xlab)) xlab = expression(lambda ['mid'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
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
  x = t(t(x)-center)

  # determine xlim, ylim
  radius = sqrt(max(apply(t(t(x)-x0)^2,1,sum)))
  if (length(auto.scale)==1) auto.scale = rep(auto.scale,2)
  if (is.null(xlim)) xlim = c(-radius,radius)*auto.scale[1]
  if (is.null(ylim)) ylim = c(-radius,radius)*auto.scale[2]
  width = xlim[2]-xlim[1]
  height = ylim[2]-ylim[1]
  diagonal = sqrt(width*height)
  nx = ngrid*width/diagonal
  ny = ngrid*height/diagonal
  
  # determine arrows
  if (is.null(arrow.origin)) arrow.origin = c(xlim[1],ylim[1])+0.05*diagonal
  if (is.null(arrow.length)) arrow.length = 0.1*diagonal
  if (is.null(title.origin)) title.origin = c(xlim[1]+0.06*diagonal,ylim[2]-0.06*diagonal)

  # rotate and project coordinates
  x = (x%*%rot)[,1:2]

  # raster data
  smoothing = 0.008*sigma*ngrid
  img = array(0,c(nx,ny,3))
  s = unique(species)
  
  for (i in seq_along(s)) {
    
    sel = species==s[i]
    if (smoothing==0) {
      g = cooltools::griddata2(x[sel,1],x[sel,2],xlim=xlim,ylim=ylim,n=c(nx,ny))
    } else {
      if (kde) {
        g = cooltools::kde2(x[sel,1],x[sel,2],xlim=xlim,ylim=ylim,n=c(nx,ny),s = 0.2*smoothing^0.5, sd.max=smoothing*2)
        g$n = g$d
      } else {
        g = cooltools::griddata2(x[sel,1],x[sel,2],xlim=c(-1,1)*radius,ylim=c(-1,1)*radius,n=c(nx,ny))
        g$n = EBImage::gblur(g$n,smoothing)
      }
    }

    for (k in seq(3)) {
      img[,,k] = img[,,k]+col[s[i],k]*g$n
    }

  }

  # finalize image
  img = atan(img/n*ngrid^2*lum)/pi*2
  f = 10^max(0,shadows)*2
  img = log10(f*img+1)/log10(f+1)
  img = cooltools::lim(img) # just to be sure

  # save raster image as png
  if (!is.null(pngfile)) {
    png::writePNG(cooltools::rasterflip(img),pngfile)
  }

  # determine margin
  for (i in seq(4)) {
    if (i%in%side) mar[i]=3
  }

  # show on screen and save as pdf
  for (mode in seq(2)) {

    make = FALSE

    if (mode==1 & screen) {
      make = TRUE
    }

    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      grDevices::pdf(pdffile,width=7*width/diagonal,height=7*height/diagonal)
    }

    if (make) {

      # initialize plot
      cooltools::nplot(xlim=xlim, ylim=ylim, asp=1, mar=rep(0,4))

      # plot raster
      rasterImage(cooltools::rasterflip(img),xlim[1],ylim[1],xlim[2],ylim[2])

      # arrows
      if (arrows) {
        if (!is.null(xlab)) {
          text(arrow.origin[1]+arrow.length,arrow.origin[2],xlab,pos=4,col='white',cex=cex)
          arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1]+arrow.length,arrow.origin[2],col='white',length = 0.1,angle=20,lwd=arrow.lwd)
        }
        if (!is.null(ylab)) {
          text(arrow.origin[1],arrow.origin[2]+arrow.length,ylab,pos=3,col='white',cex=cex)
          arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1],arrow.origin[2]+arrow.length,col='white',length = 0.1,angle=20,lwd=arrow.lwd)
        }
      }
      
      # title
      if (!is.null(title)) text(title.origin[1], title.origin[2], title, pos=4, col='white', offset=-0.4)

      # side axes
      par(xpd=T)
      rect(xlim[1],ylim[1],xlim[2],ylim[2])
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
  
  # return results
  invisible(list(rgb = img, radius = radius, rotationmatrix = rot))

}
