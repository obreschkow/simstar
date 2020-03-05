#' Show one or seveal 3D point sets in four projections
#' 
#' @importFrom cooltools nplot
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics lines
#'
#' @description Produces a 2-by-2 panel figure, where each panel is made of an \code{\link{sphview}} image in a different rotation.
#'
#' @param x list of n-by-3 matrices, represending the coordinates of n 3D positions; n can be different for each matrix. The particles stored in different matrices represent different types of particles and are visualised in different colors.
#' @param rotations 4-vector or 4-list specifying four different values/vectors of the argument \code{rotation} in \code{\link{sphview}}. These different rotations are shown in the bottom left, bottom right, top left, and top right panel, respectively, in this order.
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param title Text to be added to the figure.
#' @param ... additional parameters for \code{\link{sphview}}.
#' 
#' @seealso \code{\link{sphview}}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Example of three particle species shown in four projections
#' x.red = cooltools::runif3(1e5,polarangle = c(0,pi/2), azimuth = c(0,pi*0.75))
#' x.blue = cooltools::runif3(1e5,polarangle = c(pi/2,pi), azimuth = c(pi*1.25,2*pi))
#' x.green = t(t(cooltools::runif3(2e4,r=0.5))+c(1,1/4,-1))
#' sphview4(list(x.blue , x.red, x.green), radius.scale=1.2)
#'
#' @export

sphview4 = function(x, rotations=c(2,3,4,1), screen = TRUE, pdffile = NULL, title = NULL, ...) {

  for (mode in seq(2)) {

    make = FALSE
    if (mode==1 & screen) {
      make = TRUE
    }
    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      grDevices::pdf(pdffile,width=7,height=7)
    }

    if (make) {

      cooltools::nplot(xlim=c(0,2), ylim=c(0,2), pty='s', mar=c(0,0,0,0))

      ix = iy = 0

      for (i in seq(4)) {

        sphview(x, rotation = rotations[[i]], pngfile = NULL, pdffile = NULL, add = TRUE, add.xlim = c(0,1)+ix, add.ylim = c(0,1)+iy, screen = TRUE, ...)
        ix = (ix+1)%%2
        if (ix==0) iy = iy+1

      }

      # lines between panels
      lines(c(1,1),c(0,2),col='grey')
      lines(c(0,2),c(1,1),col='grey')
      
      # title
      if (!is.null(title)) text(0.06, 1.94, title, pos=4, col='white', offset=-0.4)

      if (mode==2) {dev.off()}

    }

  }

}