\name{riskCpGpattern}
\alias{riskCpGpattern}

\title{
	Plot the detected association pattern with one phenotype using heatmap
}
\description{
    The detected association pattern is a way to visualize the p-values providied by the HIRE function.
}
\usage{
riskCpGpattern(pval_matr, main_title = "Detected association pattern", hc_row_ind = FALSE)
}

\arguments{
  \item{pval_matr}{
	the p-value matrix for one phenotype, where one row represents a CpG site and one column indicates one cell type. 
}
  \item{main_title}{
	the title name. The default is "Detected association pattern".
}
  \item{hc_row_ind}{
	whether we conduct hierarchical clustering in the row. The default is FALSE. 
}
}
\details{
	This function depends on the heatmap.2 function in the gplots R package.  
}

\value{
	return a heatmap 
}

\author{
 Xiangyu Luo
}

\examples{
#a p-value matrix from the uniform distribution
pvalues <- matrix(runif(600), 100, 6)

#Visualize this p-value matrix
riskCpGpattern(pvalues, 
		main_title="An example", hc_row_ind = FALSE)
}

\keyword{ HIRE }
\keyword{ riskCpGpattern }
