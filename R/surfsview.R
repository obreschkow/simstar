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
#' @param simulation simulation name; must be matched by a simulation name in the parameterfile.txt in the path.surfsuite directory. This simulation name specified the actual N-body simulation, as well as the VELOCIRAPTOR output.
#' @param species optional vector listing the species of particles to be shown
#' @param fourprojections logical flag specifying whether the group is visualised in a single projedction (using \code{\link{sphview}}) or in four projections (\code{\link{sphview4}}).
#' @param ... additional arguments to passed to \code{\link{sphview}} or \code{\link{sphview4}}, respectively.
#' 
#' @seealso \code{\link{sphview}} and \code{\link{sphview4}}
#'
#' @author Danail Obreschkow
#'
#' @export

surfsview = function(haloid = 1, subhalos = T, snapshot = 199, at = NULL,
                     simulation = 'L210_N1024-Hydro3D', species,
                     fourprojections = FALSE, ...) {
  
  # check existence of surfsuite
  if (is.na(paths()$surfsuite)) stop('Please specify a surfsuite directory using paths(surfsuite = "...").')
  
  # load halo
  fn = 'tmphalo.hdf'
  if (is.null(at)) {
    call = sprintf('%ssurfsuite gethalo %d -simulation %s -snapshot %d -subhalos %d -center 1 -parameterfile %sparameters.txt -outputfile %s',
                 paths()$surfsuite,haloid,simulation,snapshot,as.integer(subhalos),paths()$surfsuite,fn)
  } else {
    call = sprintf('%ssurfsuite trackhalo %d -simulation %s -snapshot %d -from %d -to %d -subhalos %d -center 1 -parameterfile %sparameters.txt -outputfile %s',
                   paths()$surfsuite,haloid,simulation,snapshot,at,at,as.integer(subhalos),paths()$surfsuite,fn)
  }
  system(call)
  dat = rhdf5::h5read(fn,'/')
  call = paste0('rm ',fn)
  system(call)
  
  # recast positions into matrix
  if (is.null(at)) {
    x = cbind(dat$particles$rx,dat$particles$ry,dat$particles$rz)
  } else {
    x = cbind(dat$particles$snapshot$rx,dat$particles$snapshot$ry,dat$particles$snapshot$rz)
  }
  
  # visualize
  if (fourprojections) {
  
    sphview4(x, dat$particles$species, ...)  
    
  } else {
    
    sphview(x, dat$particles$species, ...)
    
  }
  
}
