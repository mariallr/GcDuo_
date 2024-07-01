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
#' reprocessing
#'
#' Process using parafac with the identified compounds obtained
#'
#' @param data original dataset
#' @param DuoResults object obtained by procesData function
#' @param lib library used for identification
#' @param win_width window width for peaks
#' @return A List
#' @export
#'
ConsProcessing <- function(GcDuoObject, DuoResults, lib, win_width = 10)
  {
  require(dplyr)
  suppressWarnings(library(doSNOW))

  start <- Sys.time()
  lib_matrix <- lib$spectra_matrix

  lib_matrix <- lib_matrix[,colnames(lib_matrix) %in% GcDuoObject$mz]

  #if no ID has been done
  if(length(DuoResults) == 3){
    DuoResults$id_peaks <- data.frame("Name" = paste0("Unknown", seq(1, nrow(results$cor_spectra))),
                                      "id" = DuoResults$cor_spectra$id,
                                      "RT1" = DuoResults$cor_spectra$time1D,
                                      "RT2" = DuoResults$cor_spectra$time2D)
  }

  # Peak coordenates

  DuoPeaks <- list()

  for(i in 1:nrow(DuoResults$id_peaks)){
    DuoPeaks[[i]] <- list(peak_coord = i,
                         RT1D = DuoResults$id_peaks$RT1[i],
                         RT2D = DuoResults$id_peaks$RT2[i],
                         sim = DuoResults$id_peaks$sim[i],
                         rsim = DuoResults$id_peaks$rsim[i],
                         time2d = GcDuoObject$time2d,
                         Name = DuoResults$id_peaks$Name[i],
                         cor_spectra = DuoResults$cor_spectra[i, ],
                         win_width = win_width)
  }



  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = length(DuoPeaks), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  FinalResults <- list()

  if(foreach::getDoParRegistered()){
    FinalResults <- foreach::foreach(DuoPeak = DuoPeaks, .options.snow = opts) %dopar%
      ConstrainedParafac(DuoPeak, GcDuoObject, DuoResults, lib_matrix)
  }

  finish <- Sys.time()
  cat(paste0("Finished processing. \n"))
  cat(paste0("Time elapsed: ", format(difftime(finish, start, units = "min"), digits = 3),". \n"))

  temp_peak <- data.frame(t(sapply(FinalResults,"[[",1)))

  #temp_peak <- sapply(FinalResults,"[[",1)
  #temp_peak <- do.call(rbind.data.frame, temp_peak)
  temp_ar <- data.frame(t(sapply(FinalResults,"[[",3)))
  temp_int <- data.frame(t(sapply(FinalResults,"[[",4)))
  #temp_spek <- sapply(FinalResults,"[[",2)
  #temp_spek <- do.call(rbind.data.frame, temp_spek)

  temp_spek <- data.frame(t(sapply(FinalResults,"[[",2)))

  return(finalGCDuo = list(peaks = temp_peak,
                            area = temp_ar,
                            intensity = temp_int,
                            spectra = temp_spek))

}


#######

ConstrainedParafac <- function(DuoPeak, GcDuoObject, DuoResults, lib_matrix){

  temp <- list(temp_peak = NULL,
               temp_spek = NULL,
               temp_ar = matrix(0, ncol = length(GcDuoObject$files)),
               temp_int = matrix(0, ncol = length(GcDuoObject$files)))

    rt1_pos <- which(GcDuoObject$time1d == DuoPeak$RT1D)
    rt2_pos <- which(GcDuoObject$time2d == DuoPeak$RT2D)

    #select the window to extract the peak
    if(rt2_pos > (length(DuoPeak$time2d)) - DuoPeak$win_width) {

      rt2_pos <- length(DuoPeak$time2d) - DuoPeak$win_width

    } else if (rt2_pos < DuoPeak$win_width) {

      rt2_pos <- DuoPeak$win_width
    }

    # get the spectra from the library
    if(rlang::is_empty(grep("Unknown", DuoPeak$Name))) {

      pos_comp <- sapply(DuoPeak$Name, function(x) which(x == row.names(lib_matrix)))
      pos_comp1 <- lapply(pos_comp, '[', 1)

      spec_std <- lib_matrix[unlist(pos_comp1)[1],]
      if(all(spec_std == 0)){
        spec_std <- lib_matrix[unlist(pos_comp1)[2],]
      }

      spec_std <- rbind(matrix(spec_std, ncol= length(GcDuoObject$mz), dimnames = list(NULL, GcDuoObject$mz)))

      row.names(spec_std) <- c(row.names(lib_matrix)[unlist(pos_comp1)[1]])

    } else { # if no ID we use the cor_spectra

      spec_std <- DuoPeak$cor_spectra[-c(1:3)]
      row.names(spec_std) <- DuoPeak$Name
    }

    if(all(spec_std == 0)) {
      next
    }

    #get the number of peaks to extract
    n = ifelse(is.null(nrow(spec_std)) == F, nrow(spec_std), 1) #num of components

    # Parafac fixed
    pf1 <- multiway::parafac2(GcDuoObject$data4D[,,(rt2_pos-DuoPeak$win_width):(rt2_pos+DuoPeak$win_width),rt1_pos], n, Bfixed = matrix(t(spec_std), nrow = length(GcDuoObject$mz)),
                              const = c("orthog", "nonneg", "nonneg"), nstart = 100) #ortnon

    rem <- fitted(pf1) #reconstructed matrix

    rem_tic <- apply(rem, c(1,3), sum)

    # For each component extracted of parafac
    for (i in 1:ncol(pf1$B)){ #get the information for each component
      plot(pf1$C[,1], type = "l")

      # Get the retention time 2 of the apex
      if(max(pf1$C[,i]) > 1.4) {
        if(all(sign(pf1$C[,i]) == 1)){
          rt2 <- which.max(pf1$C[,i])
        } else if(all(sign(pf1$C[,i]) == -1)) {
          rt2 <- which.min(pf1$C[,i])
        } else {
          rt2 <- which.max(abs(pf1$C[,i]))
        }

        #Get the spectra
        if(any(sign(pf1$B[,i]) == 1) && !any(sign(pf1$B[,i]) == -1)){
          load_spec <- pf1$B[,i]
          load_spec[load_spec < 0] <- 0
        } else if(any(sign(pf1$B[,i]) == -1) && !any(sign(pf1$B[,i]) == 1)) {
          load_spec <- -1*pf1$B[,i]
          load_spec[load_spec < 0] <- 0
        } else {
          load_spec <- pf1$B[,i]
          if(load_spec[which.max(abs(load_spec))] < 0){
            load_spec <- -1*load_spec
            load_spec[load_spec < 0] <- 0
          } else {
            load_spec[load_spec < 0] <- 0
          }
        }

        #summary table
        peak_inf <- data.frame("id" = paste(DuoPeak$peak_coord, i, sep = "-"),
                               "compound" = DuoPeak$Name,
                               "time1d" = GcDuoObject$time1d[rt1_pos],
                               "time2D"= GcDuoObject$time2d[(rt2_pos - DuoPeak$win_width):(rt2_pos + DuoPeak$win_width)][rt2],
                               "max_mz" = GcDuoObject$mz[which.max(load_spec)],
                               "int" = sum(rem[which.max(abs(pf1$A[[1]][,i])), ,rt2]),
                               "sn" = (2*sum(rem[which.max(abs(pf1$A[[1]][,i])), ,rt2])),
                               "sim" = DuoPeak$sim,
                               "rsim" = DuoPeak$rsim
        )

        #print(peak_inf)
        temp$temp_peak <- peak_inf
        temp$temp_spek <- load_spec

        ## Compute area and intensity
        for(f in 1:nrow(rem_tic)) { #for each file

          cut_data <- rem[f, which.max(load_spec), ]

          if(all(cut_data == 0)) {
            next
          }

          #area of the peak
          temp$temp_ar[,f] <- (bayestestR::area_under_curve(seq(1,length(cut_data)),
                                                          cut_data))

          #intestity of the peak
          temp$temp_int[,f] <- sum(rem[f, ,which.max(abs(pf1$C[,i]))])
        }
      }
    }
    return(temp)
}


# #Gaussian fitting
# apex <- which.max(cut_data)
# if(apex >= 10 & apex <= 30) {
#   left_part <- which.min(abs(cut_data[1:apex] - cut_data[apex]/2))
#   right_part <- (apex - 1) + which.min(abs(cut_data[apex:length(cut_data)] - cut_data[apex]/2))
#   sigma <- (cut_data[right_part]  - cut_data[left_part])/2
#   min_width <- sigma*3*2 # minimal expected width for a Gaussian peak with that intensity
#
#   if(abs(cut_data[1] - cut_data[length(cut_data)]) < min_width){
#
#   }
# }



