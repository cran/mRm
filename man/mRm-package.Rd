\name{mRm-package}
\alias{mRm-package}
\alias{mRm}
\docType{package}
\title{Parameter estimation and model selection in mixed Rasch models.}
\description{
This package provides routines for cML estimation and model selection in binary mixed Rasch models (Rost 1990). For a detailed discussion of the accuracy of parameter estimates and reliability of AIC- and BIC- based model selection techniques see Preinerstorfer and Formann (2011). 

The core part of the algorithm has been implemented in C++, using parts of the \emph{Scythe Statistical Library} (2007) for matrix manipulations. Rows with missing values and constant rows are excluded. 
}

\author{
David Preinerstorfer

\email{david.preinerstorfer@univie.ac.at}

\url{http://homepage.univie.ac.at/david.preinerstorfer}.
}

\references{
  Pemstein, D., Quinn, K. M. and Martin, A. D. (2007). 
  Scythe Statistical Library: Application programmer's interface. 
  \url{http://scythe.berkeley.edu/}. 

  Preinerstorfer, D. and Formann, A. K. (2012)
  Parameter recovery and model selection in mixed Rasch models.
  \emph{British Journal of Mathematical and Statistical Psychology, 65}, 251-262. 

  Rost (1990).
  Rasch models in latent classes: An integration of two approaches to item analysis. 
  \emph{Applied Psychological Measurement, 14}, 271-282. 
}

