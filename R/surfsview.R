#' Show VELOCIraptor groups in SURFS N-body simulation
#' 
#' @importFrom rhdf5 h5read
#'
#' @description Dedicated function to visualise data in ICRAR's SURFS simulation suite, processed by the VELOCIraptor halo-finder. Requires the external tool surfsuite, available at https://github.com/obreschkow/surfsuite.
#'
#' @param haloid halo index in the VELOCIraptor output
#' @param subhalos options specifying whether subhalos should be included (only works, if the specified VELOCIraptor output contains substructure)
#' @param snapshot simulation snapshot index
#' @param at optional snapshot index, specifying at what snapshot the particles of the halo extracted from snapshot "snapshot" should be displayed
#' @param parameterset parameterset name; must be matched by a parameterset name in the parameterfile.txt in the path.surfsuite directory. This parameterset name specifies the actual N-body simulation, as well as the VELOCIRAPTOR output.
#' @param species optional vector listing the species of particles to be shown
#' @param fourprojections logical flag specifying whether the group is visualised in a single projedction (using \code{\link{sphview}}) or in four projections (\code{\link{sphview4}}).
#' @param ... additional arguments to passed to \code{\link{sphview}} or \code{\link{sphview4}}, respectively.
#' 
#' @seealso \code{\link{sphview}} and \code{\link{sphview4}}
#'
#' @author Danail Obreschkow
#'
#' @export

surfsview = function(haloid = 1, subhalos = TRUE, snapshot = 199, at = NULL,
                     parameterset = 'L210_N1024-Hydro3D', species = NULL,
                     fourprojections = FALSE,...) {
  
  # check existence of surfsuite
  if (is.na(paths()$surfsuite)) stop('Please specify a surfsuite directory using paths(surfsuite = "...").')
  
  # load halo
  dat = surfsread(haloid = haloid, subhalos = subhalos, snapshot = snapshot, at = at,
                  parameterset = parameterset)
  
  # recast positions into matrix
  if (is.null(at)) {
    x = cbind(dat$particles$rx,dat$particles$ry,dat$particles$rz)
  } else {
    x = cbind(dat$particles$snapshot$rx,dat$particles$snapshot$ry,dat$particles$snapshot$rz)
  }
  
  # species selection
  if (!is.null(species)) {
    sel = dat$particles$species%in%species
    x = x[sel,]
    dat$particles$species = dat$particles$species[sel]
  }
  
  # visualize
  if (fourprojections) {
  
    sphview4(x, species=dat$particles$species, ...)
    invisible('no values returned when all four projections are displayed')
    
  } else {
    
    s = sphview(x, species=dat$particles$species, ...)
    invisible(s)
    
  }
  
}
