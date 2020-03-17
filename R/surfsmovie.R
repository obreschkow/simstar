#' Generate movie for SURFS particles
#' 
#' @importFrom rhdf5 h5read
#' @importFrom celestial cosdistTravelTime
#' @importFrom magick image_read image_annotate
#' @importFrom png writePNG
#'
#' @description Turns evolving particle data produced by surfsuite into an mp4 movie. Surfsuite is a fortran tool available at https://github.com/obreschkow/surfsuite. Function requires ffmpeg to be available via the terminal.
#'
#' @param file.track filename of the particle tracking file produced by calling "surfsuite trackhalo ..." (see README file of surfsuite).
#' @param radius optional value specifying the radius (in simulation units) of the image. For square-movies, this radius is equal to half the horizontal and vertical diameter of the image. Otherwise, it is the geometric mean between these two half-diameters.
#' @param aspect aspect ratio of the movie
#' @param mp4file filename of mp4-movie
#' @param fps frames per second
#' @param H0 Hubble constant of the simulation in units of km/s/Mpc
#' @param OmegaM matter density of the universe at z=0
#' @param OmegaL dark energy density of the universe at z=0
#' @param velocity.conversion conversion factor from velocity units in the simulation to length units/Gyr (at z=0). The default corresponds to the standard of Gadget-2.
#' @param dt time interval between frames in Gyr
#' @param keep.frames logical flag specifying whether the individual frame images are to be kept in a temporary folder within the directory specified by \code{\link{paths}}.
#' @param rotation either an integer (1-6) or a 3-vector specifying a rotationn of the 3D particle positions. In case of an integer: 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,y)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specifies the rotation axis and its norm the rotation angle in radians in the positive geometric sense.
#' @param show.time logical flag specifying whether the lookback time is displayed
#' @param text.size scaling factor for the text size used for the lookback time
#' @param ... additional arguments to passed to \code{\link{sphview}}.
#' 
#' @seealso \code{\link{sphview}} and \code{\link{surfsview}}
#'
#' @author Danail Obreschkow
#'
#' @export

surfsmovie = function(file.track, radius = NULL, aspect = 16/9,
                      mp4file, fps = 60,
                      H0 = 70, OmegaM = 0.3, OmegaL = 0.7,
                      velocity.conversion = 0.00102269032*(H0/100), # [-] (velocity unit)/(length unit/Gyr) at z=0
                      dt = 0.05, # [Gyr] if set to NULL, no extrapolation is chosen
                      keep.frames = TRUE,
                      rotation = 1,
                      show.time = TRUE,
                      text.size = 1,
                      ...) {
  
  # load track
  track = rhdf5::h5read(file.track,'/')
  
  # make snapshot dataframe
  snapshot_min = track$tracking$snapshot_min
  snapshot_max = track$tracking$snapshot_max
  snapshots = data.frame(index = seq(snapshot_min, snapshot_max))
  snstr = function(sn) sprintf('snapshot_%d',sn)
  for (i in seq_along(snapshots$index)) {
    snapshots$scalefactor[i] = track$particles[[snstr(snapshots$index[i])]]$scalefactor
  }
  snapshots$z = 1/snapshots$scalefactor-1
  snapshots$t = celestial::cosdistTravelTime(z=snapshots$z, OmegaM = OmegaM, OmegaL = OmegaL, H0 = H0) # [Gyr] look-back time
  snapshots$vfactor = velocity.conversion/sqrt(snapshots$scalefactor) # sqrt needed because of Gadget's internal velocity convention in cosmological runs
  
  # make temporary directory for frames
  dir = sprintf('%sframes_halo_%d/',paths()$temporary,track$halo$id)
  call = sprintf('rm -rf %s; mkdir %s',dir,dir)
  system(call)
  
  # determine scale
  if (is.null(radius)) {
    str = snstr(snapshot_max)
    x = cbind(track$particles[[str]]$rx,track$particles[[str]]$ry,track$particles[[str]]$rz)
    x0 = apply(x,2,mean) # geometric centre
    radius = sqrt(max(apply(t(t(x)-x0)^2,1,sum)))
  }
  xlim = c(-1,1)*radius*sqrt(aspect)
  ylim = c(-1,1)*radius/sqrt(aspect)
  
  # determine scale
  if (length(rotation)==3) {
    rot = t(cooltools::rotation3(rotation))
  } else {
    if (length(rotation)!=1) stop('rotation must be an integer 1,...,6 or a real 3-vector.')
    if (rotation>3) {
      e = eigen(cooltools::quadrupole(x))$vectors
    }
    if (rotation==1) {
      rot = diag(3)
    } else if (rotation==2) {
      rot = rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    } else if (rotation==3) {
      rot = rbind(c(0,1,0),c(0,0,1),c(1,0,0))
    } else if (rotation==4) {
      rot = e[,c(1,3,2)]
    } else if (rotation==5) {
      rot = e[,c(1,2,3)]
    } else if (rotation==6) {
      rot = e[,c(2,3,1)]
    } else {
      stop('rotation must be an integer 1,...,6 or a real 3-vector.')
    }
  }
  
  # function to interpolate/extrapolate positions
  .interpolate.positions = function(t) {
    
    # t = lookback time
    
    # determine interval
    n.snapshots = dim(snapshots)[1]
    if (t>=snapshots$t[1]) {
      extrapolate = T
      i0 = 1
    } else if (t<=snapshots$t[n.snapshots]) {
      extrapolate = T
      i0 = n.snapshots
    } else {
      extrapolate = F
      i1 = which(snapshots$t<=t)[1]
      i0 = i1-1
    }
    
    # interpolate/extrapolate
    if (extrapolate) {
      
      t0 = -snapshots$t[i0]
      dt = (-t)-t0
      
      str = snstr(snapshots$index[i0])
      x0 = cbind(track$particles[[str]]$rx,track$particles[[str]]$ry,track$particles[[str]]$rz)
      v0 = cbind(track$particles[[str]]$vx,track$particles[[str]]$vy,track$particles[[str]]$vz)*snapshots$vfactor[i0]
      x = x0+v0*dt
      
    } else {
      
      t0 = -snapshots$t[i0]
      t1 = -snapshots$t[i1]
      ht = t1-t0
      dt = ((-t)-t0)/ht
      
      str = snstr(snapshots$index[i0])
      x0 = cbind(track$particles[[str]]$rx,track$particles[[str]]$ry,track$particles[[str]]$rz)
      v0 = cbind(track$particles[[str]]$vx,track$particles[[str]]$vy,track$particles[[str]]$vz)*snapshots$vfactor[i0]*ht
      
      str = snstr(snapshots$index[i1])
      x1 = cbind(track$particles[[str]]$rx,track$particles[[str]]$ry,track$particles[[str]]$rz)
      v1 = cbind(track$particles[[str]]$vx,track$particles[[str]]$vy,track$particles[[str]]$vz)*snapshots$vfactor[i1]*ht
      
      q = 2*(x0-x1)+v0+v1
      p = 3*(x1-x0)-2*v0-v1
      x = x0+v0*dt+p*dt^2+q*dt^3
      
    }
    
    return(x)
    
  }
  
  # make vector of lookback times
  if (is.null(dt)) {
    t.plot = snapshots$t
  } else {
    t.plot = rev(seq(min(snapshots$t),max(snapshots$t),dt))
  }
  
  # produce frames
  for (frame in c(54,55)){#seq_along(t.plot)) {
    
    # progress update
    cat(sprintf('Make frame %d/%d\n',frame,length(t.plot)))
    
    # interpolate positions
    if (is.null(dt)) {
      str = snstr(snapshots$index[frame])
      x = cbind(track$particles[[str]]$rx,track$particles[[str]]$ry,track$particles[[str]]$rz)
    } else {
      x = .interpolate.positions(t.plot[frame])
    }
    
    # make frame
    rgb = sphview(x, track$particles$species, screen=FALSE, rotation=rot, xlim=xlim, ylim=ylim, ...)$rgb
    
    # add text to frame
    if (show.time) {
      diagonal = sqrt(prod(dim(rgb)[1:2]))
      s = 0.03*diagonal*text.size
      rgb = magick::image_read(rgb)
      rgb = magick::image_annotate(rgb, sprintf('Lookback time = %.2f Gyr',t.plot[frame]),
                                   size = s, location = sprintf('%+d%+d',round(1.8*s),round(s)), color = 'white', font='sans', degrees=90)
      rgb = as.numeric(rgb[[1]])[,,1:3]
    }
    
    # save frame
    fn = sprintf('%sframe_%0.6d.png',dir,frame)
    writePNG(rasterflip(rgb),fn)
    
  }
  
  # convert frames into movie
  cat(sprintf('Write mp4 movie file %s\n',mp4file))
  linuxspaces = function(txt) gsub(' ','\\\\ ',txt)
  call = sprintf('rm -rf %s',mp4file)
  system(call)
  call = sprintf('ffmpeg -r %d -f image2 -s %dx%d -i %sframe_%%06d.png -vcodec libx264 -crf 18 -pix_fmt yuv420p %s -loglevel quiet',
                 fps,dim(rgb)[1],dim(rgb)[2],linuxspaces(dir),linuxspaces(mp4file))
  system(call)
  
  # delete frames
  if (!keep.frames) {
    call = sprintf('rm -rf %s',dir)
    system(call)
  }
  
}
