
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
#' visualizeChrom3D
#'
#' Get the graphically representation of the chromatogram in 2D. Only possible to visualize one file per time
#'
#' @param data defining structure of the data processing workflow
#' @param GcDuoObject asdasijdhjkajh
#' @param sampleNum the index of the file to visualize
#' @import patchwork
#' @return A ggplot contour
#' @export

visualizeChrom3D <- function(GcDuoObject, sampleNum = 1)
{
  dataunf <- array(GcDuoObject$data4D[sampleNum,,,], dim = c(
                                   dim(GcDuoObject$data4D)[2],
                                   prod(dim(GcDuoObject$data4D)[3:4])))

  ticred <- apply(dataunf, 2, sum)

  plotly::plot_ly(z = matrix(ticred, nrow = dim(GcDuoObject$data4D)[3],
                             dimnames = list(GcDuoObject$time2d, GcDuoObject$time1d)),
                  type = "contour", contours = list(showlabels = T, coloring = "OrRd"),
                  colorscale = "OrRd")
}

##########################
#' visualizeChrom2D
#'
#' Get the graphically representation of the chromatogram in 2D
#'
#' @param GcDuoObject asdasijdhjkajh
#' @param sampleNum the index of the file/s to visualize
#' @return An plotly graph
#' @export

visualizeChrom2D <- function(GcDuoObject, samples = 1:12)
{
  dataunf <- array(GcDuoObject$data4D[samples,,,], dim = c(length(samples), dim(GcDuoObject$data4D)[2],prod(dim(GcDuoObject$data4D)[3:4])))

  ticred <- apply(dataunf, c(1,3), sum)
  row.names(ticred) <- GcDuoObject$files[samples]

  times <- paste(round(sort(rep(GcDuoObject$time1d, length(GcDuoObject$time2d))),2),
                 GcDuoObject$time2d, sep = "-")

  plot_data2 <- data.frame(times, t(ticred)) |>
    pivot_longer(!times, names_to = "var", values_to = "val")

  pl <- plot_ly(plot_data2,
                x = ~times,
                y = ~val,
                color = ~var,
                split = ~var,
                type = "scatter",
                mode = "lines",
                colors = "Spectral")
  pl <- layout(pl,
               xaxis = list(categoryarray = names, categoryorder = "array", title = "RT"),
               yaxis = list(title = "TIC intensity"),
               legend = list(title = "Sample", y = -0.5, orientation = "h"))
  pl
}

##########################
#' PlotPeak
#'
#' Get the graphically representation of a detected peak. Only possible to visualize one peak per time
#'
#' @param GcDuoObject asdasijdhjkajh
#' @param finalGCDuo final table with the peaks detected
#' @param peakid index of the selected peak to visualize
#' @param type select if the visualization is for tic or for eic
#' @return A ggplot multiline
#' @export

PlotPeak <- function(GcDuoObject, finalGCDuo, peakid = "68-68-1", type = "eic", win_width = 20){
  require(patchwork)

  peak_id <- which(finalGCDuo$peaks$id == peakid)
  rt1_pos <- which(GcDuoObject$time1d == finalGCDuo$peaks$time1d[peak_id])
  rt2_pos <- which(GcDuoObject$time2d == finalGCDuo$peaks$time2D[peak_id])

  nb.cols <- length(GcDuoObject$files)
  mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "PuBuGn"))(nb.cols)


  if(type == "eic") {
    var1 = readline(prompt = 'Enter an m/z : ')
    mz_pos <- which(GcDuoObject$mz == var1)

    raw_eic <- GcDuoObject$data4D[,mz_pos,(rt2_pos-win_width):(rt2_pos+win_width),rt1_pos]
    row.names(raw_eic) <- GcDuoObject$files

    spec_std <- finalGCDuo$spectra[peak_id,]

    pf1 <- multiway::parafac(GcDuoObject$data4D[,,(rt2_pos-win_width):(rt2_pos+win_width),rt1_pos], 1, Bfixed = matrix(t(spec_std), nrow = length(GcDuoObject$mz)),
                             const = c("ortnon", "nonneg", "nonneg"), nstart = 100)

    rem <- fitted(pf1) #reconstructed matrix

    rem_eic <- rem[,mz_pos,]
    row.names(rem_eic) <- GcDuoObject$files

    plot_data1 <- data.frame(colnames(raw_eic), t(raw_eic)) |>
      pivot_longer(!colnames.raw_eic., names_to = "var", values_to = "val")

    pl1 <- ggplot(plot_data1, aes(x = colnames.raw_eic., y = val, group = var, color = var)) +
      geom_line() + theme_minimal() +
      xlab("") + ylab("Intensity") + ggtitle("EIC Raw") +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_color_manual(values = mycolors)

    plot_data2 <- data.frame(colnames(raw_eic), t(rem_eic)) |>
      pivot_longer(!colnames.raw_eic., names_to = "var", values_to = "val")

    pl2 <- ggplot(plot_data2, aes(x = colnames.raw_eic., y = val, group = var, color = var)) +
      geom_line() + theme_minimal() +
      xlab("Retention time 2 (s)") + ylab("Intensity") + ggtitle("EIC GcDuo") +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_color_manual(values = mycolors)

    fig <- pl1 / pl2

    fig <- fig + plot_annotation(title = finalGCDuo$peaks$compound[peak_id],
                          subtitle = paste('RT1: ',round(unlist(finalGCDuo$peaks$time1d[peak_id]),2),"s",
                                           ' - m/z frag.: ',GcDuoObject$mz[mz_pos],  sep = ""))
  } else if (type == "tic"){

    raw_eic <- apply(GcDuoObject$data4D[,,(rt2_pos-win_width):(rt2_pos+win_width),rt1_pos], c(1,3), sum)
    row.names(raw_eic) <- GcDuoObject$files

    spec_std <- finalGCDuo$spectra[peak_id,]

    pf1 <- multiway::parafac(GcDuoObject$data4D[,,(rt2_pos-win_width):(rt2_pos+win_width),rt1_pos], 1, Bfixed = matrix(t(spec_std), nrow = length(GcDuoObject$mz)),
                             const = c("ortnon", "nonneg", "nonneg"), nstart = 100)

    rem <- fitted(pf1) #reconstructed matrix

    rem_eic <- apply(rem, c(1,3), sum)
    row.names(rem_eic) <- GcDuoObject$files

    plot_data1 <- data.frame(colnames(raw_eic), t(raw_eic)) |>
      pivot_longer(!colnames.raw_eic., names_to = "var", values_to = "val")

    pl1 <- ggplot(plot_data1, aes(x = colnames.raw_eic., y = val, group = var, color = var)) +
      geom_line() + theme_minimal() +
      xlab("") + ylab("Intensity") + ggtitle("TIC Raw") +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_color_manual(values = mycolors)

    plot_data2 <- data.frame(colnames(raw_eic), t(rem_eic)) |>
      pivot_longer(!colnames.raw_eic., names_to = "var", values_to = "val")

    pl2 <- ggplot(plot_data2, aes(x = colnames.raw_eic., y = val, group = var, color = var)) +
      geom_line() + theme_minimal() +
      xlab("Retention time 2 (s)") + ylab("Intensity") + ggtitle("TIC GcDuo") +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_color_manual(values = mycolors)

    fig <- pl1 / pl2

    fig <- fig + plot_annotation(title = finalGCDuo$peaks$compound[peak_id],
                                 subtitle = paste('RT1: ',round(unlist(finalGCDuo$peaks$time1d[peak_id]),2),"s",  sep = ""))
  }


  return(fig)
}
