% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sits_label_classification.R
\name{sits_label_classification}
\alias{sits_label_classification}
\title{Post-process a classified data raster probs to obtain a labelled image}
\usage{
sits_label_classification(
  cube,
  smoothing = NULL,
  multicores = 1,
  memsize = 1,
  output_dir = "./",
  version = "v1"
)
}
\arguments{
\item{cube}{Classified image data cube.}

\item{smoothing}{(deprecated)}

\item{multicores}{Number of process to label the classification in
snow subprocess.}

\item{memsize}{Maximum overall memory (in GB) to label the
classification.}

\item{output_dir}{Output directory where to out the file}

\item{version}{Version of resulting image
(in the case of multiple tests)}
}
\value{
A data cube
}
\description{
Takes a set of classified raster layers with probabilities,
             and label them based on the maximum probability for each pixel.
}
\examples{
\dontrun{
# Retrieve the samples for Mato Grosso
# select band "ndvi"

samples_ndvi <- sits_select(samples_mt_4bands, bands = "NDVI")

# select a random forest model
rfor_model <- sits_train(samples_ndvi, sits_rfor(num_trees = 500))

# Classify a raster file with 23 instances for one year
files <- c(system.file("extdata/raster/mod13q1/sinop-crop-ndvi.tif",
    package = "sits"
))

# create a data cube based on the information about the files
sinop <- sits_cube(
    type = "BRICK", satellite = "TERRA",
    sensor = "MODIS", name = "Sinop-crop",
    timeline = timeline_modis_392,
    output_dir = tempdir(),
    bands = c("NDVI"), files = files
)

# classify the raster image
sinop_probs <- sits_classify(sinop,
    ml_model = rfor_model,
    output_dir = tempdir(),
    memsize = 4, multicores = 2
)

# label the classification and smooth the result with a bayesian filter
sinop_label <- sits_label_classification(sinop_probs, output_dir = tempdir())
}

}
\author{
Rolf Simoes, \email{rolf.simoes@inpe.br}
}
