#########################################################################
#     GcDuo - R package for GCxGC processing and PARAFAC analysis
#     Copyright (C) 2023 Maria Llambrich & Lluc Semente
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

##########################
#' readFolderCDF
#'
#' Read all CDF files in a folder and transform the data into a list of 3D arrays.
#'
#' @param folderPath path to the folder containing the CDF to process.
#' @param modulationTime modulation time used in the data.
#' @param mzRange range values used in the m/z data acquisition.
#' @return A List of arrays containing all the data in the files following the structure: |mz, time2D, time1D|
#' @export
#'
readFolderCDF <- function(folderPath, modulationTime, mzRange)
{
  # Get files path
  filenames <- list.files(folderPath,recursive = FALSE, pattern = "*.cdf")
  if (length(filenames) <= 1)
  {
    stop("Found less than 2 CDF in the path")
  }
  cat(paste("Found", length(filenames),"files in folder. \n"))

  # Setup results list
  results <- vector('list', length(filenames))

  # Processing loop
  start <- Sys.time()
  for (f in 1:length(filenames))
  {
    cat(paste0("Reading and processing file ", f,"/",length(filenames),". \n"))
    results[[f]] <- readCFD(filePath = paste0(folderPath, filenames[f]),
                            modulationTime = modulationTime,
                            mzRange = mzRange)
    cat(paste("\n"))
  }

  cat(paste("All files read. \n"))
  cat(paste("Uniforming dimentions between files.\n"))
  results <- uniformDimentions(results, mzRange)
  cat(paste("All processing complete. \n"))

  finish <- Sys.time()
  cat(paste0("Time elapsed: ", format(difftime(finish, start, units = "min"), digits = 3),". \n"))
  gc(verbose = FALSE)
  results$files <- filenames
  return(results)
}


##########################
#' readCFD
#'
#' Read a single CDF file and transform the data into a 3D array.
#'
#' @param pathFolderCDF path to the folder containing the CDF to process.
#' @param modulationTime modulation time used in the data.
#' @param mzRange range values used in the m/z data acquisition.
#' @return A 3 dimensional array containing all the data in the files following the structure: |sample, mz, time2D, time1D|
readCFD <- function(filePath, modulationTime, mzRange)
{
  suppressWarnings(library(doParallel))
  # File information Extraction
  nc <- ncdf4::nc_open(filePath)
  int <- ncdf4::ncvar_get(nc, "intensity_values", start = 1, count = -1) #variable for intensity values
  time <- ncdf4::ncvar_get(nc, "scan_acquisition_time") #variable acquisition time
  mz <- ncdf4::ncvar_get(nc, "mass_values") #variable mass to charge ratio values
  tic <- ncdf4::ncvar_get(nc, "total_intensity") #variable total intensity
  point <- ncdf4::ncvar_get(nc, "point_count")
  ncdf4::nc_close(nc) #close connection to cdf file

  N <- length(time)

  time <- time/60
  f_mostr <- abs(time[1] - time[2]) #frequency
  mod_t_ex <- round((modulationTime/60)/f_mostr, 0) * f_mostr
  cycles <- mod_t_ex/f_mostr

  dim2d <- cycles  # frequency * modulation = dimension retention time 2
  dim1d <- N/ cycles # dimension retention time 1

  mz_nodec <- round(mz,0)
  if(tail(sort(unique(mz_nodec)),1) - tail(sort(unique(mz_nodec)),2)[1] == 1){
    mz_seq <- seq(min(unique(mz_nodec)), max(unique(mz_nodec)))
  } else {
    mz_seq <- seq(min(unique(mz_nodec)), tail(sort(unique(mz_nodec)),2)[1] + 1)
  }

  period <- (max(mz_seq) + 1) - min(mz_seq) # mz period
  mz_range <- range(mz_seq) #mz range
  mzRange <- mzRange

  # clusters
  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(max = length(time), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # processing iteration
  if (foreach::getDoParRegistered())
  {
    int_list <- foreach::foreach(i = iterators::icount(length(time)), .options.snow = opts) %dopar% getSpectrumIntensity(mzs = mz, points =  point, index =  i, ints = int, mz_seq = mz_seq)
  }
  parallel::stopCluster(cl)

  # Shape data conveniently
  gc_data <- array(dim = c(length(int_list), length(int_list[[1]])))
  for (spec in 1:nrow(gc_data))
  {
    gc_data[spec, ] <- int_list[[spec]]
  }
  rm(int_list)

  # Filling time and intensity axis with zeros to ensure the folding is correct
  #needtime <- seq(0, (time[1]), by = mean(diff(time))) #correct time at the endings
  # posdelay <- length(needtime)/dim2d
  # misspos <- (trunc(posdelay) - 1) * dim2d
  # needpos <- length(needtime) - misspos
  # time_cor <- c(tail(needtime, n = needpos), time)
  time_pos <- seq(1, length(time), by = dim2d)
  time_values1d <- time[time_pos] #get time values retention time 1
  # dim1d <- length(time_values1d)
  time_values2d <- round(seq(0, (modulationTime - (modulationTime/dim2d)), length.out = dim2d), 3) #get time values retention time 2

  #create the 3D array

  array3D <- array(data = c(as.vector(t(gc_data))),
                     dim = c(period, dim2d, dim1d),
                     dimnames = list(mz_seq, time_values2d, time_values1d))

  if (length(mz_seq) > length(seq(mzRange[1], mzRange[2]))) {
    mz1 <- which(mz_seq == mzRange[1])
    mz2 <- which(mz_seq == mzRange[2])

    message(paste("mz dimension reduced from ", mz_range[1], "-", mz_range[2],
                  "to ", mzRange[1], "-", mzRange[2]))

    array3D <- array3D[mz1:mz2,,]

  } else {
    mz1 <- 1
    mz2 <- length(mz_seq)
  }

  rm(gc_data)

  results <- list(data = array3D, mz_seq = mz_seq[mz1:mz2], time_values2d = time_values2d, time_values1d = time_values1d)

  return(results)
}


##########################
#' getSpectrumIntensity
#'
#' Composes an intensity spectrum
#'
getSpectrumIntensity <- function(mzs, points, index, ints, mz_seq)
{
  point <- points[index] #number of mz
  intensity <- rep(0, times = length(mz_seq))

  if (point != 0)
  {
    position <- sum(points[1:(index - 1)])
    mz <- mzs[position + (1:point)]
    mz <- round(mz)
    mz_U <- unique(mz)
    int <- ints[position + (1:point)]
    if (length(mz_U) < length(mz))
    {
      for (j in 1:length(mz_U))
      {
        intensity[which(mz_seq == mz_U[j])] <- sum(int[which(mz == mz_U[j])])
      }
    } else
      {
        for (j in 1:length(mz))
        {
          intensity[which(mz_seq == mz_U[j])] <- int[j]
        }
      }
  }
  return(intensity)
}


##########################
#' uniformDimentions
#'
#' Uniforms the range of all dimentions between CDF files
#'
uniformDimentions <- function(data3DList, mzRange)
{
  # uniforming time1d axis
  max_time1d <- -1
  min_time1d <- 9999
  time1d_step <- c()
  for (i in 1:length(data3DList))
  {
    time1d_step <- c(time1d_step, mean(diff(data3DList[[i]]$time_values1d)))
    timeRange <- range(data3DList[[i]]$time_values1d)

    if (timeRange[1] < min_time1d)
    {
      min_time1d <- timeRange[1]
    }

    if (timeRange[2] > max_time1d)
    {
      max_time1d <- timeRange[2]
    }
  }

  # uniforming mz axis
  max_mz <- -1
  min_mz <- 9999
  mz_step <- c()
  for (i in 1:length(data3DList))
  {
    mz_step <- c(mz_step, mean(diff(data3DList[[i]]$mz_seq)))
    #mzRange <- range(data3DList[[i]]$mz_seq)
    mzRange <- mzRange

    if (mzRange[1] < min_mz)
    {
      min_mz <- mzRange[1]
    }

    if (mzRange[2] > max_mz)
    {
      max_mz <- mzRange[2]
    }
  }

  #Uniformed dimentional axis
  time1d <- seq(from = min_time1d, to = max_time1d, by = mean(time1d_step, na.rm = T))
  mz <- seq(from = min_mz, to = max_mz, by = mean(mz_step, na.rm =T))
  time2d <- data3DList[[1]]$time_values2d

  #Final data structure to fill
  GcDuoObject <- list(data4D = array(data  = 0,
                                     dim = c(length(data3DList), length(mz),
                                             length(time2d), length(time1d)),
                                     dimnames = list(1:length(data3DList), mz, time2d, time1d)),
                      time1d = time1d,
                      mz = mz,
                      time2d = time2d)

  #Loop filling the GcDuoObject
  for (i in 1:length(data3DList))
  {
    df <- array(data = 0, dim = c(length(mz), length(time2d), length(time1d)))

    if (abs(length(data3DList[[i]]$time_values1d) - length(time1d)) > 1 |
        abs(length(data3DList[[i]]$mz_seq) - length(mz)) > 100 )
    {
      stop("Test: The difference between mass axies or time1d axies in the files is to big. We need to develope new code for it!")
    }

    if (length(data3DList[[i]]$time_values1d) != length(time1d) & length(data3DList[[i]]$mz_seq) != length(mz))
    {
      timeRange <- range(data3DList[[i]]$time_values1d)
      mzRange <- range(data3DList[[i]]$mz_seq)
      if (abs(timeRange[1] - min_time1d) >  mean(time1d_step)/2  & abs(mzRange[1] - min_mz) > mean(mz_step)/2)
      {
        df[-1, , -1] <- data3DList[[i]]$data
      }
      else if (abs(timeRange[1] - min_time1d) >  mean(time1d_step)/2  & abs(mzRange[1] - min_mz) < mean(mz_step)/2)
      {
        df[-length(mz), , -1] <- data3DList[[i]]$data
      }
      else if (abs(timeRange[1] - min_time1d) <  mean(time1d_step)/2  & abs(mzRange[1] - min_mz) > mean(mz_step)/2)
      {
        df[-1, , -length(time1d)] <- data3DList[[i]]$data
      }
      else if (abs(timeRange[1] - min_time1d) < mean(time1d_step)/2  & abs(mzRange[1] - min_mz) < mean(mz_step)/2)
      {
        df[-length(mz), , -length(time1d)] <- data3DList[[i]]$data
      }
    }
    else if (length(data3DList[[i]]$time_values1d) != length(time1d))
    {
      timeRange <- range(data3DList[[i]]$time_values1d)

      if (abs(timeRange[1] - min_time1d) > 1 )
      {
        df[ , , -1] <- data3DList[[i]]$data
      }
      else if (abs(timeRange[2] - max_time1d) > 1)
      {
        df[ , , -length(time1d)] <- data3DList[[i]]$data
      }
    }
    else if (length(data3DList[[i]]$mz_seq) != length(mz))
    {
      mzRange <- range(data3DList[[i]]$mz_seq)
      mz_pos <- which(data3DList[[i]]$mz_seq %in% mz)

      df <- data3DList[[i]]$data[mz_pos,,]
    }
    else
    {
      df <- data3DList[[i]]$data
    }

    GcDuoObject$data4D[i, , , ] <- df
    rm(df)
  }

  rm(data3DList)
  gc(verbose = FALSE)
  return(GcDuoObject)
}
