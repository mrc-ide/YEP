% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_setup.R
\name{template_region_xref}
\alias{template_region_xref}
\title{template_region_xref}
\usage{
template_region_xref(template = list(), regions = c())
}
\arguments{
\item{template}{List containing one or more sets of template data (serological data, case data or burden data)}

\item{regions}{Vector of individual regions}
}
\description{
Cross-reference template data with individual regions
}
\details{
Examines template data (serological, case or burden) and compares with vector of region names in order
 to check which lines in each set of template data require model data from which region(s). For example, if a
 line in a set of serological data has its region given as "AGO.1_1,AGO.2_1,AGO.3_1" and the compared vector of
 regions is c("AGO.1_1","AGO.2_1","AGO.3_1",...), then that line requires data from regions 1, 2 and 3.

 This function is used when generating a dataset from one or more templates; it is normally used by the functions
 Generate_Dataset, Generate_VIMC_Burden_Dataset and Generate_Multiple_Datasets. It returns a list containing [TBA].
}
