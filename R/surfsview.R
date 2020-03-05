#' Show VELOCIraptor groups in SURFS N-body simulation
#' 
#' @importFrom rhdf5 h5read
#'
#' @description Dedicated function to visualise data in ICRAR's SURFS simulation suite, processed by the VELOCIraptor halo-finder. Requires the external tool surfsuite, available at https://github.com/obreschkow/surfsuite.
#'
#' @param haloid halo index in the VELOCIraptor output
#' @param subhalos options specifying whether subhalos should be included (only works, if the specified VELOCIraptor output contains substructure)
#' @param snapshot simulation snapshot index
#' @param simulation simulation name; must be matched by a simulation name in the parameterfile.txt in the path.surfsuite directory. This simulation name specified the actual N-body simulation, as well as the VELOCIRAPTOR output.
#' @param species optional vector listing the species of particles to be shown
#' @param lum overall luminosity scaling factor (default 1).
#' @param shadows differential luminosity scaling factor for darker regions (default 1).
#' @param fourprojections logical flag specifying whether the group is visualised in a single projedction (using \code{\link{sphview}}) or in four projections (\code{\link{sphview4}}).
#' @param ... additional arguments to passed to \code{\link{sphview}} or \code{\link{sphview4}}, respectively.
#' 
#' @seealso \code{\link{sphview}} and \code{\link{sphview4}}
#'
#' @author Danail Obreschkow
#'
#' @export

surfsview = function(haloid = 1, subhalos = T, snapshot = 199,
                     simulation = 'L210_N1024-Hydro3D', species,
                     lum = 1, shadows = 1, fourprojections = FALSE, ...) {
  
  # check existence of surfsuite
  if (substr(paths()$surfsuite,1,1)=='[') stop('Please specify a surfsuite directory using paths(surfsuite = "...").')
  
  # load halo
  fn = 'tmphalo.hdf'
  call = sprintf('%ssurfsuite gethalo %d -simulation %s -snapshot %d -subhalos %d -center 1 -parameterfile %sparameters.txt -outputfile %s',
                 paths()$surfsuite,haloid,simulation,snapshot,as.integer(subhalos),paths()$surfsuite,fn)
  system(call)
  dat = rhdf5::h5read(fn,'/')
  call = paste0('rm ',fn)
  system(call)
  
  # extract positions into matrix
  x = cbind(dat$particles$rx,dat$particles$ry,dat$particles$rz)
  
  # separate two species
  xlist = list()
  species = sort(unique(dat$particles$species))
  for (i in seq(length(species),1)) {
    xlist[[i]] = x[dat$particles$species==species[i],]
  }
  
  # visualize
  if (fourprojections) {
  
    sphview4(xlist, lum = 0.2*lum, shadows=1.5*shadows, ...)  
    
  } else {
    
    sphview(xlist, lum = 0.2*lum, shadows=1.5*shadows, ...)
    
  }
  
}
