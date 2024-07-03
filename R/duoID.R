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
#' duoID
#'
#' Compares spectras with library
#'
#' @param DuoResults object obtained from processData
#' @param lib dataset of spectras to compare.
#' @param match_cutoff value to consider a correct match (range 0 to 1)
#' @param RI retention index of reference compounds (alkanes or FAMEs) to use, time must be in seconds
#' @param RIrange tolerance of RI calculation
#' @return data table with the matched peaks
#' @export
#'
duoID <- function(DuoResults, lib, match_cutoff = 0.7, RI = NULL, RIrange = 30)
{
  require(dplyr)
  suppressWarnings(library(doSNOW))

  RI = RI; RIrange = RIrange;

  if (!is.null(RI) & length(lib) <= 1)
  {
    stop("Not valid library included. Please use...")
  }

  lib_matrix <- lib$spectra_matrix

  spec_mean_mat <- as.matrix(DuoResults$cor_spectra[,-c(1:3)])

  start <- Sys.time()

  ## check m/z ranges of both matrices

  if(!identical(range(as.numeric(colnames(lib_matrix))), range(as.numeric(colnames(spec_mean_mat))))){
    if(min(as.numeric(colnames(lib_matrix))) < min(as.numeric(colnames(spec_mean_mat))) &
       max(as.numeric(colnames(lib_matrix))) > max(as.numeric(colnames(spec_mean_mat)))){

      lib_matrix <- lib_matrix[,colnames(lib_matrix) %in% colnames(spec_mean_mat)]

    } else if (min(as.numeric(colnames(lib_matrix))) > min(as.numeric(colnames(spec_mean_mat))) &
                 max(as.numeric(colnames(lib_matrix))) < max(as.numeric(colnames(spec_mean_mat)))){

      spec_mean_mat <- spec_mean_mat[,colnames(spec_mean_mat) %in% colnames(lib_matrix)]
    } else {
      min_mz <- max(c(min(as.numeric(colnames(lib_matrix))), min(as.numeric(colnames(spec_mean_mat)))))

      max_mz <- min(c(max(as.numeric(colnames(lib_matrix))), max(as.numeric(colnames(spec_mean_mat)))))

      lib_matrix <- lib_matrix[,colnames(lib_matrix) %in% seq(min_mz, max_mz)]
      spec_mean_mat <- spec_mean_mat[,colnames(spec_mean_mat) %in% seq(min_mz, max_mz)]
      }
  }

  # Get the samples RI based in RI excel info

  if(!is.null(RI)){

    samples_ri <- getRI(DuoResults, RI)

    libs_ri <- lib$data

  } else {
    #samples_ri <- NULL

    libs_ri <- NULL
  }

  Npeak <- list()
  for(i in 1:nrow(DuoResults$cor_spectra)){
    Npeak[[i]] <- list(spectra = spec_mean_mat[i,],
                       pos = i,
                       match_cutoff = match_cutoff,
                       id = DuoResults$cor_spectra$id[i],
                       time1D = DuoResults$cor_spectra$time1D[i],
                       time2D = DuoResults$cor_spectra$time2D[i],
                       RI = RI,
                       RIrange = RIrange)
  }

  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = nrow(spec_mean_mat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  lib_peaks <- list()

  if(!is.null(RI)){

    if(foreach::getDoParRegistered()){
      lib_peaks <- foreach::foreach(Peak2ID = Npeak, .options.snow = opts) %dopar% libID(Peak2ID,
                                                                                         lib_matrix,
                                                                                         libs_ri,
                                                                                         samples_ri)
    }
  } else {

    if(foreach::getDoParRegistered()){
      lib_peaks <- foreach::foreach(Peak2ID = Npeak, .options.snow = opts) %dopar% libID(Peak2ID, lib_matrix)
    }
  }

  parallel::stopCluster(cl)


    lib_peaks <- do.call(rbind.data.frame, lib_peaks)

    # Sort results by RT
    DuoResults$id_peaks <- lib_peaks |> arrange(as.numeric(RT1), as.numeric(RT2))
    DuoResults$cor_spectra <- DuoResults$cor_spectra |> arrange(as.numeric(time1D), as.numeric(time2D))

    finish <- Sys.time()
    cat(paste0("Finished processing. \n"))
    cat(paste0("Time elapsed: ", format(difftime(finish, start, units = "min"), digits = 3),". \n"))
    gc(verbose = FALSE)

    return(DuoResults)
}

##########################
#' libID
#'
#' Annotation of spectras
#'
libID <- function(Peak2ID, lib_matrix) {

    spec_one <- Peak2ID$spectra

    cor_cut <- NULL
    for(i in 1:nrow(lib_matrix)){
      cor_cut[i] <- lsa::cosine(as.numeric(lib_matrix[i,]), as.numeric(spec_one))
    }

    names(cor_cut) <- row.names(lib_matrix)

  #Get the similarity values
  #Only matches over the cutoff
  if(any(cor_cut >= Peak2ID$match_cutoff, na.rm = T)) {

    #compute reverse similarity
    id_table <- id_rev(cor_cut, spec_one, Peak2ID$match_cutoff, lib_matrix)

    dif_sims <- id_table$rsim - id_table$sim

    id_table <- id_table[which(dif_sims >= 0),]

    id_table.s <- id_table[order(id_table$sim, id_table$rsim, decreasing = T),]

    # select ID by RI
    if(!is.null(Peak2ID$RI)){

      # Adding RI to the identifications from the library
      ris <- sapply(id_table.s$compound, function(x){
        a <- which(libs_ri$comp == x)
        return(libs_ri$ri[a[1]])})

      id_table.s$t_ri <- as.numeric(ris)

      # Get the library ID with the closest RI to the sample
      e_ri <- samples_ri[Peak2ID$pos]

      id_table.ri <- id_table.s[id_table.s$t_ri < e_ri + Peak2ID$RIrange &
                                  id_table.s$t_ri > e_ri - Peak2ID$RIrange,]

      id_table.ri <- id_table.ri[complete.cases(id_table.ri),]

      if(nrow(id_table.ri) >= 1){
        pos <- which.max(id_table.ri$sim)
        id_table.s <- id_table.ri
      } else {
        pos <- which.max(id_table.s$sim)
      }
    } else {
      pos <- which.max(id_table.s$sim)
    }

    # Report the results
    if(!rlang::is_empty(pos)) {
      process <- data.frame("Name" = id_table.s$compound[pos],
                            "id" =  Peak2ID$id,
                            "RT1" = Peak2ID$time1D,
                            "RT2" = Peak2ID$time2D,
                            "sim" = id_table.s$sim[pos],
                            "rsim" = id_table.s$rsim[pos]
      )
    } else {
      process <- data.frame("Name" = paste("Unknown_", Peak2ID$pos),
                            "id" =  Peak2ID$id,
                            "RT1" = Peak2ID$time1D,
                            "RT2" = Peak2ID$time2D,
                            "sim" = 0,
                            "rsim" = 0
      )
    }

  } else {
    process <- data.frame("Name" = paste("Unknown", Peak2ID$pos),
                          "id" =  Peak2ID$id,
                          "RT1" = Peak2ID$time1D,
                          "RT2" = Peak2ID$time2D,
                          "sim" = 0,
                          "rsim" = 0
    )
  }

  if(!is.null(Peak2ID$RI)){

    process$ri <- samples_ri[Peak2ID$pos]

  }

  return(process)

}

##########################
#' ID inverse
#'
#' Performs reverse similarity
#'
id_rev <- function(cor_cut, spec_sample, match_cutoff, lib_matrix) {

  id_table <- data.frame()
  #only compouds above the cut-off
  high <- which(cor_cut >= match_cutoff)

  #for each compound
  for (n in 1:length(high)){

    pos <- high[n]

    comp <- names(cor_cut)[pos]

    # reverse

    libmat <- lib_matrix[pos, ]

    libmat.f <- libmat[libmat > 0]

    sampspe.f <- c(t(spec_sample[names(spec_sample) %in% names(libmat.f)]))

    revmat <- (sampspe.f %*% libmat.f) / outer(sqrt(sum(sampspe.f^2)),
                                               sqrt(sum(libmat.f^2)))

    posid <- data.frame("compound" = comp, "sim" = cor_cut[pos],
                        "rsim" = revmat)

    id_table <- rbind(id_table, posid)
  }

  return(id_table)
  finish <- Sys.time()
  cat(paste0("Finished processing. \n"))
  cat(paste0("Time elapsed: ", format(difftime(finish, start, units = "min"), digits = 3),". \n"))
}

##########################
#' RI calculation
#'
#' Calculates Kovats Retention Index based in Alkanes
#'
rt2ri <- function(rtTime, lib) {

  vecrt <- sort(c(rtTime, lib$RT))

  pos <- which(vecrt == rtTime)[1]

  ri_e <- OrgMassSpecR::RetentionIndex(lib$Carbon[pos-1], rtTime, lib$RT[pos-1], lib$RT[pos])

  return(ri_e)
}

##########################
#' Get Retention Index
#'
#' Calculates Kovats Retention Index based in Alkanes
#'
getRI <- function(DuoResults, RI){

  ris <- NULL
  #for each peak in the list compute the retention index based in alkanes
  for(i in 1:nrow(DuoResults$cor_spectra)){
    ri_res <- rt2ri(DuoResults$cor_spectra$time1D[i], RI)
    if(rlang::is_empty(ri_res)){
      ris <- c(ris, NA)
    } else {
      ris <- c(ris, ri_res)
    }
  }

  DuoResults$cor_spectra$ri <- ris
  DuoResults$cor_spectra$ids <- seq(1, nrow(DuoResults$cor_spectra))

  # for peaks before first alkane and after last alkane, linear model
  feat_ri <- DuoResults$cor_spectra |>
    filter(!is.na(ri))

  feat_na <- DuoResults$cor_spectra |>
    filter(is.na(ri))

  model <- lm(ri ~ time1D, data = feat_ri)

  feat_na$ri <- predict(model, feat_na)

  samples_ri <- tibble(rbind(feat_ri, feat_na)) |>
    arrange(ids) |>
    pull(ri)

  return(samples_ri)
}

########





