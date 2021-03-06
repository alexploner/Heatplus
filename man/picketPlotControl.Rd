% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annHeatmap.R
\name{picketPlotControl}
\alias{picketPlotControl}
\title{Default parameter settings for picketPlot}
\usage{
picketPlotControl()
}
\value{
A named list
}
\description{
This function returns a named list of parameters that affect how a
picketPlot is generated. This list can be used as a template for 
overriding the defaults partially or completely.
}
\details{
The following parameter affects the overall appearance of the plot:
\itemize{ 
  \item\code{cex.label} is the expansion factor for the size of the cluster labels at
    the bottom of the plot; default is 1.5.
}

The following parameters directly affect how binary indicator variables are 
displayed:
\itemize{ 
  \item\code{boxw} is the relative length of the short side of a box marking 
     (width for a horizontal plot); default is 1.
  \item\code{boxh} is the relative length of the long side of a box marking
    (default: 4) 
  \item\code{hbuff} is the relative distance between two box markings for the 
    same variable (horizontal buffer for a horizontal plot); default is 0.1 
  \item\code{vbuff} is the relative distance between two box
     markings for the same subject, but different variables (default: 0.1)
  \item\code{nacol} is the color for box markings indicating missing values 
    (default: \code{gray(0.85)}) 
 }
 Note that \code{boxh} and \code{vbuff} also affect the display of numerical
 variables as a scatter plot: as the amount of vertical space allowed for a 
 single numerical variable (see also \code{numfac} below) and the vertical 
 space between two neighboring variable panels (binary or ornumerical), 
 respectively.

The following parameters only affect the display of a numerical variable:
\itemize{
  \item\code{numfac} is the expansion factor indicating how much higher (for a
    horizontal plot) or wider (for a vertical plot) panels with numerical
    variables are than a panels for a single binary indicator
  \item\code{span} is the span argument for the loess smoother. Default is 1/3; 
    setting this to zero switches off smoothing.
  \item\code{degree} is the degree of loess smoothing. Default is 1; setting
    this to zero switches off smoothing 
  \item\code{pch} is the plotting character for numerical variables; uses the 
    device default. 
  \item\code{cex.pch} is the size of the plotting character for numerical 
    variables; uses the device default. 
  \item\code{col.pch} is the color of the plotting character for numerical 
    variables; uses the device default.
  \item\code{label_axis_shrink} controls the range of the axis for which axis 
    ticks are labeled: by default, labels covering the whole observed range 
    are defined (via a call to \code{pretty}); if set to a number between zero
    and one, the range covered by labels is shortened by that fraction and 
    centered within the observed range (and the fed to \code{pretty}); this can 
    be used to avoid overlapping labels for multiple adjacent panes with 
    numerical variables. 
  \item\code{plot_baseline} is a logical value indicating whether to draw 
    a baseline for panes showing numerical variables: FALSE by default, this
    can be useful to visually separate multiple adjacent panes with numerical
    variables. 
}
}
\seealso{
\code{\link{picketPlot}}, \code{\link{par}}, \code{\link{pretty}}
}
