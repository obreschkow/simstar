#' Load VELOCIraptor groups in SURFS N-body simulation
#' 
#' @importFrom rhdf5 h5read
#'
#' @description Dedicated function to load halos in cosmological N-body simulations, identified by the VELOCIraptor halo-finder. Requires the external tool surfsuite, available at https://github.com/obreschkow/surfsuite.
#'
#' @param haloid halo index in the VELOCIraptor output
#' @param subhalos options specifying whether subhalos should be included (only works, if the specified VELOCIraptor output contains substructure)
#' @param snapshot simulation snapshot index
#' @param at optional snapshot index, specifying at what snapshot the particles of the halo extracted from snapshot "snapshot" should be displayed
#' @param parameterset parameterset name; must be matched by a parameterset name in the parameterfile.txt in the path.surfsuite directory. This parameterset name specifies the actual N-body simulation, as well as the VELOCIRAPTOR output.
#' 
#' @seealso \code{\link{surfsview}} and \code{\link{surfsmovie}}
#'
#' @author Danail Obreschkow
#'
#' @export

surfsread = function(haloid = 1, subhalos = TRUE, snapshot = 199, at = NULL,
                     parameterset = 'L210_N1024-Hydro3D') {
  
  # check existence of surfsuite
  if (is.na(paths()$surfsuite)) stop('Please specify a surfsuite directory using paths(surfsuite = "...").')
  
  # load halo
  fn = 'tmphalo.hdf'
  if (is.null(at)) {
    call = sprintf('%ssurfsuite gethalo %d -parameterset %s -snapshot %d -subhalos %d -center 1 -parameterfile %sparameters.txt -outputfile %s',
                   paths()$surfsuite,haloid,parameterset,snapshot,as.integer(subhalos),paths()$surfsuite,fn)
  } else {
    call = sprintf('%ssurfsuite trackhalo %d -parameterset %s -snapshot %d -from %d -to %d -subhalos %d -center 1 -parameterfile %sparameters.txt -outputfile %s',
                   paths()$surfsuite,haloid,parameterset,snapshot,at,at,as.integer(subhalos),paths()$surfsuite,fn)
  }
  system(call)
  dat = rhdf5::h5read(fn,'/')
  call = paste0('rm ',fn)
  system(call)
  
  # convert data
  if (is.null(at)) {
    dat$particles = data.frame(dat$particles)
  }
  
  # return data
  return(dat)
  
}
