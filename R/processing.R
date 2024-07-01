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
#' procesData
#'
#' Process using parafax a section of the data defined in the dataSlot
#'
#' @param data defining structure of the data processing workflow
#' @param GcDuoObject raw GCxGC-MS data uploaded with `readFolderCDF`
#' @param signalNoiseRatio1 minimum signal to noise ratio 3D
#' @param signalNoiseRatio2 minimum signal to noise ratio 2D
#' @param exclusionIon mass to charge fragments to exclude in the analysis
#' @param topAbundance percetage of mass to change fragments to consider in the model creation
#' @param win_size number of modulations to include in each window
#' @return A List
#' @export
#'
procesData <- function(GcDuoObject, signalNoiseRatio1 = 10, signalNoiseRatio2 = 500, exclusionIon = NULL, topAbundance = 0.05, win_size = 3)
{
  suppressWarnings(library(doSNOW))
  start <- Sys.time()

  # time windows
  window_T1d <- list()
  for (i in 1:(length(GcDuoObject$time1d) - (win_size - 1)))
  {
    window_T1d <- c(window_T1d, list(i:(i + (win_size - 1))))
  }

  # data slots
  dataSlots <- list()
  for (i in 1:length(window_T1d))
  {
    dataSlots[[i]] <- list(peaks_win = GcDuoObject$data4D[,,,window_T1d[[(i)]]],
                           t1 = window_T1d[[(i)]],
                           time_values1d = GcDuoObject$time1d,
                           time_values2d = GcDuoObject$time2d,
                           sampleSize = dim(GcDuoObject$data4D)[1],
                           mz_seq = GcDuoObject$mz,
                           signalNoiseRatio1 = signalNoiseRatio1,
                           signalNoiseRatio2 = signalNoiseRatio2,
                           exclusionIon = exclusionIon,
                           topAbundance = topAbundance,
                           win_size = win_size)
  }
  cat(paste0("Generated ", length(window_T1d)," data slots. \n"))
  cat(paste0("Parallel processing begins. \n"))

  # clusters
  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(max = length(window_T1d), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # processing iteration
  results <- list()
  if (foreach::getDoParRegistered())
  {
    results <- foreach::foreach(slot = dataSlots, .options.snow = opts) %dopar% procesDataSlot(slot)
  }
  parallel::stopCluster(cl)

  # results reshaping
  ind_peak <- NULL
  spec_ind <- NULL
  cat(paste0("\nMerging data slots results. \n"))
  for (slot in 1:length(results))
  {
    if (dim(results[[slot]]$ind_peak)[1] != 0)
    {
      if (is.null(ind_peak))
      {
        ind_peak <- results[[slot]]$ind_peak
        spec_ind <- results[[slot]]$spec_ind
      } else
        {
          ind_peak <- rbind(ind_peak, results[[slot]]$ind_peak)
          spec_ind <- rbind(spec_ind, results[[slot]]$spec_ind)
        }
    }
  }

  ind_peak <- ind_peak[!duplicated(ind_peak$id),]
  spec_ind <- spec_ind[!duplicated(ind_peak$id),]

  row.names(spec_ind) <- ind_peak$id


  #############################################################

  # Filter duplicates
  ind_peak_m <- ind_peak[ ,c("id", "time1D", "time2D", "mz_max", "sn", "int")]
  ind_peak_m$time1D <- as.factor(ind_peak_m$time1D)
  ind_peak_m$mz_max <- as.factor(ind_peak_m$mz_max)
  ind_peak_m <- ind_peak_m[ind_peak_m$sn > 1, ]
  ind_peak_m <- ind_peak_m |>
    dplyr::group_by(mz_max, time1D) |>
    dplyr::summarise("id" = paste(id, collapse = "/"),
                     "n" = dplyr::n(),
                     "sn" = median(sn)) |>
    dplyr::filter(n >= 3,
                  sn > 2*signalNoiseRatio2)

  # Combine peaks found in different samples
  common_mz <- NULL
  spec_mean <- data.frame()
  peaks_single <- data.frame()
  for (t in 1:nrow(ind_peak_m))
  {
    ids <- unlist(strsplit(ind_peak_m$id[t], "/")) # peak position
    dt <- ind_peak |> dplyr::filter(id %in% ids) # peaks information
    dt_spec <- spec_ind |> dplyr::filter(row.names(spec_ind) %in% dt$id) #spectra info

    #separate max mz fragments found
    if (nrow(dt) > 2)
    {
      mz_dt <- dt |> dplyr::select("mz_top") |>
        tidyr::separate(mz_top, c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")) |>
        unlist(use.names = F)
      common_mz <- rbind(common_mz, unique(mz_dt)) #detect repeated fragments among samples
      spec_tofix <- apply(dt_spec, 2, mean) #compute the mean between samples
      spec_mean <- rbind(spec_mean, spec_tofix) #save info for spectra
      peaks_single <- rbind(peaks_single, dt[1,2:4]) #save peaks info
    }
  }

  # Perform the aggregate spectra from the combined peaks
  spec_mean_cor <- data.frame() # Filtering spectra with low values
  for (m in 1:nrow(peaks_single))
  {
    dif_p <- apply(common_mz, 1, function(x) match(x, common_mz[m,]))
    del_mz <- which(dif_p[, -m] < 4)
    mz_t <- common_mz[-m, ]
    delmz <- mz_t[del_mz]
    lowmatch <- delmz[match(delmz, common_mz[m, ]) >= 5]
    spec_del <- which(GcDuoObject$mz %in% as.numeric(lowmatch))
    s1 <- spec_mean[m, ]
    s1[ , spec_del] <- 0
    spec_mean_cor <- rbind(spec_mean_cor, s1)
  }

  colnames(spec_mean_cor) <-GcDuoObject$mz
  cor_spectra <- cbind(peaks_single, spec_mean_cor)

  ##########################################################################

  # TODO At the end: results <- list(spec_mean_cor)

  DuoResults <- list(ind_peak = ind_peak, spec_ind = spec_ind, cor_spectra = cor_spectra)

  # return
  finish <- Sys.time()
  cat(paste0("Finished processing. \n"))
  cat(paste0("Time elapsed: ", format(difftime(finish, start, units = "min"), digits = 3),". \n"))
  gc(verbose = FALSE)
  return(DuoResults)
}


##########################
#' procesDataSlot
#'
#' Process using parafax a section of the data defined in the dataSlot
#'
#' @param dataSlot defining structure of the data processing workflow
#' @return A List
#'
procesDataSlot <- function(processingStructure)
{
  suppressWarnings(library(magrittr))
  ind_peak <- data.frame() # intensity vectors
  spec_ind <- data.frame() # mz vectors

  # TIC and noise computation
  databb_tic <- getTICfromMatrix(processingStructure$peaks_win, 4)
  databb_mt <- apply(databb_tic, 2, mean) # mean of all samples used as quality control
  if (databb_mt[1] > databb_mt[length(databb_mt)])
  {
    noise <- mean(tail(databb_mt, 5))
  } else
  {
    noise <- mean(databb_mt[1:5])
  }

  if (noise == 0)
  {
    noise <- min(databb_mt[databb_mt > 0])
  }

  # Watershed algorithm to detect peak regions
  rwater <- EBImage::watershed(matrix(databb_mt, ncol = processingStructure$win_size), tolerance = 20)
  blobs <- unique(c(t(rwater)))
  blobs <- blobs[blobs > 0]

  # Blob analysis
  if (length(blobs) > 1)
  {
    times <- NULL
    for (i in 1:length(processingStructure$time_values1d[processingStructure$t1]))
    {
      times <- c(times, paste(rep(round(processingStructure$time_values1d[processingStructure$t1][i], processingStructure$win_size),
                                  length(processingStructure$time_values2d)), processingStructure$time_values2d, sep = "-"))
    }

    b2t2 <- NULL # Get the retention time id for the window
    for (i in 1:(length(blobs)))
    {
      blob <- blobs[c(i, i + 1)] # we select two blobs
      blob_p <- which((rwater %in% blob) == T, arr.ind = T)
      dt <- databb_mt[blob_p] # get position of blobs
      t2t2 <- rep(processingStructure$time_values2d, processingStructure$win_size)[blob_p] # convert to retention time
      sn <- 2*max(dt)/noise # check peak noise level
      apex <- which.max(dt) # check peak shape
      left_part <- which.min(abs(dt[1:apex] - dt[apex]/2))
      right_part <- (apex - 1) + which.min(abs(dt[apex:length(dt)] - dt[apex]/2))
      sigma <- (t2t2[right_part]  - t2t2[left_part])/2
      min_width <- sigma*3*2 # minimal expected width for a Gaussian peak with that intensity

      if (abs(t2t2[1] - t2t2[length(dt)]) < min_width)
      {
        blob <- blobs[i + 2]
        blob_p <- c(blob_p, which((rwater %in% blob) == T, arr.ind = T))
        blob_p <- sort(blob_p)
      }

      if (sn > processingStructure$signalNoiseRatio1 & length(blob_p) > 4) # filter small peaks
      {
        t_blob <- times[blob_p]
        st <- strsplit(t_blob, "-")
        st2 <- sapply(st, '[[',2)
        new_b <- list(which(processingStructure$time_values2d %in% unique(st2)))

        if (!all(unlist(new_b) %in% unlist(b2t2[[length(b2t2)]])))
        {
          b2t2 <- c(b2t2, new_b)
        }
      }
    }

    b2t2 <- b2t2[!duplicated(b2t2)] #remove duplicated
    if (length(b2t2) > 1)
    {
      for (bl in 1:length(b2t2))
      {
        t2 <- b2t2[[bl]]
        for (s in 1:processingStructure$sampleSize) #for each sample
        {
          wind_check <- "yes"
          data <- processingStructure$peaks_win[s, , t2, 1:processingStructure$win_size] #data of only one sample

          if (!is.null(processingStructure$exclusionIon))
          {
            # TODO
            # Search the mz od the exclusionIon in the mass axis
            # which.min(abs(mass-exclusionIon[i]))
            # mzToExclude <- function(processingStructure$exclusionIon)
            #data[mzToExclude,,] <- 0
          }

          data[data < 0] <- 0
          data_raw <- processingStructure$peaks_win[s, ,t2,1:processingStructure$win_size] #get the raw data
          data_raw_t <- array(data_raw, dim = c(dim(data)[1], prod(dim(data)[2:3]))) #array to matrix
          var_mz <- apply(data_raw_t,1, var) #variation of each m/z
          var_5 <- order(var_mz, decreasing = T)[1:round(length(var_mz)*processingStructure$topAbundance)] #select the m/z with higher variation
          data_reduced <- data[var_5,,] #reduce data
          data_reduced_t <- array(data_reduced, dim = c(dim(data_reduced)[1], prod(dim(data_reduced)[2:3]))) #get the tic
          data_reduced_t_tic <- apply(data_reduced_t, 2, sum)
          min_start <- data_reduced_t_tic[1]
          min_end <- data_reduced_t_tic[length(data_reduced_t_tic)]

          if (min_start > min_end)
          {
            noise <- mean(tail(data_reduced_t_tic, 5))
          } else
          {
            noise <- mean(data_reduced_t_tic[1:5])
          }

          if (noise == 0)
          {
            noise <- min(data_reduced_t_tic[data_reduced_t_tic > 0])
          }

          p <- findPeaks(data_reduced_t_tic) #find peaks (2D)

          sn <- (2*data_reduced_t_tic[p])/noise #signal-to-noise
          t_p <- which(sn > processingStructure$signalNoiseRatio2) #only big peaks
          if (all(data == 0))
          {
            next
          }
          else if (length(t_p) >= 1)
          {
            parf_par <- data.frame()

            n <- 1 #iterative parafac model
            pf1 <- multiway::parafac(data_reduced, n,
                                     const = c("ortnon", "nonneg", "nonneg"),
                                     nstart = 100, verbose = FALSE)

            parf_p <- data.frame("n" = n,
                                 "Rsq" = pf1$Rsq,
                                 "SSE" = pf1$SSE,
                                 "iter" = pf1$iter)

            parf_par <- rbind(parf_par, parf_p)
            tcc <- rrcov3way::congruence(pf1$B)
            tcc[col(tcc) == row(tcc)] <- 0

            while (any(tcc > 0.90) == F)
            {
              n <- n + 1
              pf1 <- multiway::parafac(data_reduced, n,
                                       const = c("ortnon", "nonneg", "nonneg"),
                                       nstart = 100, verbose = FALSE)


              parf_p <- data.frame("n" = n,
                                   "Rsq" = pf1$Rsq,
                                   "SSE" = pf1$SSE,
                                   "iter" = pf1$iter)

              parf_par <- rbind(parf_par, parf_p)

              tcc <- rrcov3way::congruence(pf1$B)
              tcc[col(tcc) == row(tcc)] <- 0
            }

            n <- ifelse(n != 1, n - 1, 1) # perform parafac to all data with the components obtained
            pf1 <- multiway::parafac(data, n,
                                     const = c("ortnon", "nonneg", "nonneg"),
                                     nstart = 100, verbose = FALSE)
            rem <- fitted(pf1) #reconstructed matrix after parafac
            data1_t <- array(rem, dim = c(dim(rem)[1], prod(dim(rem)[2:3])))
            data1_tic <- apply(data1_t, 2, sum) # compute tic
            times <- NULL
            for (i in 1:length(processingStructure$time_values1d[processingStructure$t1]))
            {
              times <- c(times, paste(rep(round(processingStructure$time_values1d[processingStructure$t1][i],processingStructure$win_size),
                                          length(processingStructure$time_values2d[t2])),
                                      processingStructure$time_values2d[t2], sep = "-"))
            }

            p <- findPeaks(data1_tic) # detect peak points

            sn <- (2*data1_tic[p])/noise # signal to noise
            t_p <- which(sn > processingStructure$signalNoiseRatio2) # peaks over threshold
            p_sel <- p[t_p] #points selected

            # Get peak information and spectra
            peak_infs <- NULL
            spec_inds <- NULL
            if (length(p_sel) > 0)
            {
              if (n == 1)
              { # if only one compound detected
                mzs <- data.frame("mz" = processingStructure$mz_seq, "int" = pf1$A)
                peak_inf <- data.frame("id" = paste(s,processingStructure$t1[[1]], bl,1, sep = "-"),
                                       "time1D" = processingStructure$time_values1d[processingStructure$t1][which.max((pf1$C))],
                                       "time2D" = processingStructure$time_values2d[t2][which.max((pf1$B))],
                                       "mz_max" = mzs$mz[which.max(mzs$int)],
                                       "sn" = (2*sum(rem[,which.max((pf1$B)),which.max((pf1$C))]))/noise,
                                       "int" = mean(rem[,which.max((pf1$B)), which.max((pf1$C))]),
                                       "area_q1" = bayestestR::area_under_curve(seq(1, length(t2)),rem[which.max(mzs$int),,which.max((pf1$C))]),
                                       "area_q3" = bayestestR::area_under_curve(seq(1, length(t2)),rem[order(mzs$int, decreasing = T)[3],,which.max((pf1$C))]),
                                       "mz_top" = paste(mzs[order(mzs$int, decreasing = T),1][1:10], collapse = "-"),
                                       "sample" = s)
                peak_infs <- rbind(peak_infs, peak_inf)
                s_spec <- data.frame("mz" = mzs$mz[which(pf1$A != 0)],
                                     "int" = mzs[which(pf1$A != 0), 2])
                if (s_spec$mz[which.max(s_spec$int)] != peak_inf$mz_max)
                {
                  break
                }
                if(length(-which(processingStructure$mz_seq %in% s_spec$mz)) != length(processingStructure$mz_seq)) {
                  need <- data.frame("mz" = processingStructure$mz_seq[-which(processingStructure$mz_seq %in% s_spec$mz)], "int" = 0)
                  complete_mz <- rbind(s_spec, need)
                } else {
                  complete_mz <- s_spec
                }
                spec_inds <- rbind(spec_inds, t(complete_mz[order(complete_mz$mz), "int"]))
                colnames(spec_inds) <- processingStructure$mz_seq
              }
              else
              {
                for (j in 1:n)
                { # for more than 1 compound detected
                  samp_f <- which(pf1$A[,j] != 0)
                  if (length(samp_f) < 2)
                  {
                    next
                  }
                  mzs <- data.frame("mz" = processingStructure$mz_seq, "int" = pf1$A[,j])
                  peak_inf <- data.frame("id" = paste(s,processingStructure$t1[[1]], bl,j, sep = "-"),
                                         "time1D" = processingStructure$time_values1d[processingStructure$t1][which.max((pf1$C[,j]))],
                                         "time2D" = processingStructure$time_values2d[t2][which.max((pf1$B[,j]))],
                                         "mz_max" = mzs$mz[which.max(mzs$int)],
                                         "sn" = (2*sum(rem[,which.max((pf1$B[,j])), which.max((pf1$C[,j]))]))/noise,
                                         "int" = sum(rem[,which.max((pf1$B[,j])), which.max((pf1$C[,j]))]),
                                         "area_q1" = bayestestR::area_under_curve(seq(1,length(t2)),rem[which.max(mzs$int), ,which.max((pf1$C[,j]))]),
                                         "area_q3" = bayestestR::area_under_curve(seq(1,length(t2)),rem[order(mzs$int, decreasing = T)[3], ,which.max((pf1$C[,j]))]),
                                         "mz_top" =  paste(mzs[order(mzs$int, decreasing = T), 1][1:10], collapse = "-"),
                                         "sample" = s)

                  if (peak_inf$time2D != processingStructure$time_values2d[t2][1])
                  {
                    peak_infs <- rbind(peak_infs, peak_inf)
                    s_spec <- data.frame("mz" = mzs$mz[which(pf1$A[,j] != 0)], "int" = mzs[which(pf1$A[,j] != 0), 2])

                    if(length(-which(processingStructure$mz_seq %in% s_spec$mz)) != length(processingStructure$mz_seq)) {
                      need <- data.frame("mz" = processingStructure$mz_seq[-which(processingStructure$mz_seq %in% s_spec$mz)], "int" = 0)
                      complete_mz <- rbind(s_spec, need)
                    } else {
                      complete_mz <- s_spec
                    }

                    spec_inds <- rbind(spec_inds, t(complete_mz[order(complete_mz$mz), "int"]))
                    colnames(spec_inds) <- processingStructure$mz_seq
                  } else
                  {
                    wind_check <- "no"
                  }
                }
              }

              if(is.null(peak_infs) & is.null(spec_inds)){
                wind_check <- "no"
              }

              # as it is a rolling window remove duplicates
              #peaks_infs <- NULL
              if (wind_check == "no")
              {
                peaks_infs <- NULL
                spec_inds <- data.frame()
              } else if (!is.null(peak_infs)){
                if (nrow(peak_infs) > 1)
                {
                  peaks_infs <- peak_infs |> tibble::rownames_to_column("rowid") |>
                    dplyr::group_by(time1D, time2D, int, mz_max) |>
                    dplyr::slice(1L)
                  spec_inds <- spec_inds[as.numeric(peaks_infs$rowid),]
                } else if (nrow(peak_infs) == 1)
                {
                  peaks_infs <- peak_infs |> tibble::rownames_to_column("rowid")
                  spec_inds <- t(spec_inds[as.numeric(peaks_infs$rowid),])
                }

              }

              if (!is.null(peaks_infs))
              {
                if (nrow(peaks_infs) == 1)
                {
                  ind_peak <- rbind(ind_peak, peaks_infs)
                  spec_ind <- rbind(spec_ind, spec_inds)
                  colnames(spec_ind) <- processingStructure$mz_seq
                }
                else if (nrow(peaks_infs) == dim(spec_inds)[1])
                {
                  ind_peak <- rbind(ind_peak, peaks_infs)
                  spec_ind <- rbind(spec_ind, spec_inds)
                  colnames(spec_ind) <- processingStructure$mz_seq
                } else {
                }
              }
            }
          } else
            {
          }
        }
      }
    }
  }
  return(list(ind_peak = ind_peak, spec_ind = spec_ind))
}


##########################
#' getTICfromMatrix
#'
#' Computes TIC from a data array
#'
getTICfromMatrix <- function(datamatrix, dimension)
{
  if (dimension == 4)
  {
    reshape_matrix <- array(as.vector(datamatrix), dim = c(dim(datamatrix)[1], dim(datamatrix)[2], prod(dim(datamatrix)[3:4])))
    tic <- apply(reshape_matrix, c(1,3), sum)
    return(tic)
  } else if (dimension == 3)
  {
    reshape_matrix <- array(datamatrix, dim = c(dim(datamatrix)[1], prod(dim(datamatrix)[2:3])))
    tic <- apply(reshape_matrix, 2, sum)
    return(tic)
  }
}


##########################
#' findPeaks
#'
#' Finds peaks in a data array
#'
findPeaks <- function(x, m = 3)
{
  shape <- diff(sign(diff(x, na.pad = FALSE))) #Detect if the next values are increasing
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  return(pks)
}

