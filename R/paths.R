#' Get and set paths to external routines
#' 
#' @param surfsuite path of the tool surfsuite, available at https://github.com/obreschkow/surfsuite
#' @param procorr path of the tool procorr, available at https://github.com/obreschkow/procorr
#' @param temporary path to save temporary data
#' 
#' @return Returns a list with the currently set paths.
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Check current paths
#' paths()
#'
#' @export

paths = function(surfsuite=NULL, procorr=NULL, temporary=NULL) {
  
  if (!is.null(surfsuite)) {
    surfsuite = gsub('//','/',paste0(surfsuite,'/'))
    if (file.exists(surfsuite)) {
      if (file.exists(paste0(surfsuite,'surfsuite'))) {
        .simstar.env$path.surfsuite = surfsuite
      } else {
        stop('executable "surfsuite" does not exist in surfsuite directory')
      }
    } else {
      stop(paste0('directory does not exist: ',surfsuite))
    }
  }
  
  if (!is.null(procorr)) {
    procorr = gsub('//','/',paste0(procorr,'/'))
    if (file.exists(procorr)) {
      if (file.exists(paste0(procorr,'procorr'))) {
        .simstar.env$path.procorr = procorr
      } else {
        stop('executable "procorr" does not exist in procorr directory')
      }
    } else {
      stop(paste0('directory does not exist: ',procorr))
    }
  }
  
  if (!is.null(temporary)) {
    temporary = gsub('//','/',paste0(temporary,'/'))
    .simstar.env$path.temporary = temporary
    if (!file.exists(temporary)) {
      call = sprintf('mkdir -p %s',temporary)
      system(call)
    }
  }
  
  return(list(surfsuite = ifelse(!is.null(.simstar.env$path.surfsuite),.simstar.env$path.surfsuite,NA),
              procorr = ifelse(!is.null(.simstar.env$path.procorr),.simstar.env$path.procorr,NA),
              temporary = ifelse(!is.null(.simstar.env$path.temporary),.simstar.env$path.temporary,NA)))
  
}
