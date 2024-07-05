
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
#' library_import
#'
#' Read libraries in msp format and transform to GcDuo
#'
#' @param lib_list path for the msp file/s
#' @import doParallel
#' @return a S3 object
#' @export
#'
#'
library_import <- function(msp_file){

  li <- NULL

  message("Reading")

  #Get library name
  l_id <- unlist(stringr::str_split(msp_file, "/"))
  library_id <- stringr::str_split(l_id[length(l_id)], "\\.")[[1]][1]

  msp <- readLines(msp_file)
  # remove empty lines
  msp <- msp[msp != '']
  #number of compounds
  ncomp <- grep('^NAME:', msp, ignore.case = TRUE)

  # Divide by entry
  splitFactorTmp <- rep(1:length(ncomp), diff(c(ncomp, length(msp) + 1)))

  li <- split(msp,f = splitFactorTmp)

  suppressWarnings(library(doSNOW))

  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(max = length(li), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  libs <- lapply(li, getMSP)

  # processing iteration
  if (foreach::getDoParRegistered())
  {
    libs <- foreach::foreach(lib = li, .options.snow = opts) %dopar% getMSP(lib)
  }

  libs_dic <- NULL
  for(i in 1:length(libs)){
    libs_dic <- rbind(libs_dic, data.frame("comp" = libs[[i]]$name,
                                           "exactmass" = ifelse(
                                             length(libs[[i]]$exactmass) == 0,
                                             NA,
                                             libs[[i]]$exactmass),
                                           "formula" = ifelse(
                                             length(libs[[i]]$formula) == 0,
                                             NA,
                                             libs[[i]]$formula),
                                           "inchikey" = ifelse(
                                             length(libs[[i]]$inchikey) == 0,
                                             NA,
                                             libs[[i]]$inchikey),
                                           "ri" = ifelse(
                                             length(libs[[i]]$ri) == 0,
                                             NA,
                                             libs[[i]]$ri)))
  }

  parallel::stopCluster(cl)

  message("Converting to matrix")

  no_cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(no_cores - 1)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(max = length(libs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  if (foreach::getDoParRegistered())
  {
    results <- foreach::foreach(lis = libs, .options.snow = opts) %dopar% spectramatrix(lis)
  }
  parallel::stopCluster(cl)

  spec_df <- data.frame((sapply(results,"[[",1)))

  spectra_matrix <- do.call(rbind, spec_df)

  names_sp <- data.frame((sapply(results,"[[",2)))

  if(identical(names_sp[,1], libs_dic$comp)) {
    row.names(spectra_matrix) <- names_sp[,1]
  } else {
    break
  }

  lib <- list(data = libs_dic, spectra_matrix = spectra_matrix)

  return(lib)
}

#' getMSP
#'
#' Read lines of msp file
#'
#' @param file msp file to read
#' @import doParallel
#' @return a list of list
#'
getMSP <- function(x){

    namet <- x[grep('^NAME:',x, ignore.case=TRUE)]
    name <- gsub('^NAME: ','',namet, ignore.case=TRUE)
    exactmasst <- x[grep('^EXACT.MASS:|^EXACTMASS:',x, ignore.case=TRUE)]
    exactmass <- gsub('EXACT.MASS: |^EXACTMASS: ','', exactmasst, ignore.case=TRUE)
    casnum <- x[grep('^CASNO:',x, ignore.case=TRUE)]
    cas <- gsub('CASNO: ','', casnum, ignore.case=TRUE)
    formt <- x[grep('^FORMULA: ',x, ignore.case=TRUE)]
    formula <- gsub('^FORMULA: ','',formt,ignore.case = TRUE)
    inchit <- x[grep('^INCHIKEY: |^INCHI: ',x, ignore.case=TRUE)]
    inchikey <- gsub('^INCHIKEY: |^INCHI: ','',inchit,ignore.case = TRUE)
    ids <- x[grep('ID: ',x,ignore.case = TRUE)]
    id <- gsub('ID: ','',ids,ignore.case=TRUE)
    comm <- x[grep('^COMMENT: ',x,ignore.case = TRUE)]
    rts <- x[grep('ri. |retention.index. |retentionindex. ',tolower(x),ignore.case = TRUE)][1]

    # Nist
    rt_ <- unlist(regmatches(tolower(rts),
                             gregexpr(
                               'ri.*[0-9]*|retention.index.*[0-9]*|retentionindex.*[0-9]*',
                               tolower(rts))))
    rtt <- gsub('ri.|retention.index. |retentionindex. ','',rt_, ignore.case=TRUE)

    # Number of m/z fragments
    npt <- x[grep('^Num Peaks: ',x, ignore.case=TRUE)]
    np <- gsub('^Num Peaks: ','',npt,ignore.case = TRUE)

    if(as.numeric(np) > 0){
      # matrix of masses and intensities
      massIntIndx <- which(grepl('\\s[0-9]', x) & !grepl(': ', x))
      massesInts <- unlist(strsplit(x[massIntIndx], '; '))
      massesInts <- unlist(stringr::str_squish(massesInts))
      massesInts <- strsplit(massesInts, ' ')
      mz <- unlist(lapply(massesInts, '[[', 1))
      ins <- unlist(lapply(massesInts, '[[', 2))
      # if any NAs remove from indx
      spectra <- cbind.data.frame(mz=mz,ins=ins)
    }
    return(list(name=name,exactmass=exactmass, cas = cas, formula=formula, inchikey=inchikey,
                db.id = id, lib.id = library_id, ri = rtt, ri2 = rts, np = np,spectra=spectra))
}

#' spectramatrix
#'
#' Convert the spectra list to matrix
#'
#' @param libs list obtained with `getMSP` function
#' @import doParallel
#' @import tibble
#' @import dplyr
#' @return a matrix
#'
spectramatrix <- function(lis){

  seq_mz <- seq(30,600)
    spec_ref <- tibble(lis$spectra)
    spec_ref_ok <- spec_ref %>% mutate(
      mz = round(as.numeric(as.character(mz)),0),
      ins = as.numeric(as.character(ins))) %>%
      group_by(mz) %>%
      summarise(mz = mean(mz),
                ins = sum(ins)) %>%
      filter(mz >= min(seq_mz) & mz <= max(seq_mz)) %>%
      arrange(mz)
    if (nrow(spec_ref_ok) < length(seq_mz) & nrow(spec_ref_ok) != 0) {
      spec_ref_ok <- spec_ref_ok %>%
        add_row(mz = seq_mz[!(seq_mz %in% .$mz)], ins = 0) %>%
        arrange(mz)
    } else if (nrow(spec_ref_ok) == 0) {
      ref_matrix <- rep(NA, length(seq_mz))
    }

    ref_names <- lis$name
    if (nrow(spec_ref_ok) == length(seq_mz)) {

      spec_ref_ok$rel <- spec_ref_ok$ins/max(spec_ref_ok$ins)

      ref_matrix <- spec_ref_ok$rel
    }

    res <- list(ref_matrix = ref_matrix,
                ref_names = ref_names)

  return(res)
}




