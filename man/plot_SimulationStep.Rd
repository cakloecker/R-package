% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_SimulationStep}
\alias{plot_SimulationStep}
\title{This function can be called within the observer function to save a multiplot for each timestep, 
containing the landscape rasters and species abundance & richness, incl. indicated observation area
NOTE: requires plot_raster()}
\usage{
plot_SimulationStep(data, pa_matrix, abundance_matrix)
}
\arguments{
\item{data}{data list as available in Gen3sis simulation}

\item{pa_matrix}{presence-absence matrix with x and y coordinates in first cols and p/a entries per species in subsequent cols (as generated in Gen3sis observer)}

\item{abundance_matrix}{abundance matrix with x and y coordinates in first cols and abundance entries per species in subsequent cols (as generated in Gen3sis observer)}

\item{c}{list with the experiment configuration from  yaml file}
}
\value{
no return value - plot saved for each timestep
}
\description{
This function can be called within the observer function to save a multiplot for each timestep, 
containing the landscape rasters and species abundance & richness, incl. indicated observation area
NOTE: requires plot_raster()
}
