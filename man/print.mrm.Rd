\name{print.mrm}
\alias{print.mrm}
\title{Print method for Objects of Class mrm.}
\description{
 Prints arguments of an object of class mrm.
}
\usage{
\method{print}{mrm}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{mrm}.}
  \item{\dots}{Additional parameters to \code{print}.}
 }

\author{
David Preinerstorfer; 

\email{david.preinerstorfer@univie.ac.at};

\url{http://homepage.univie.ac.at/david.preinerstorfer}.
}

\examples{
#Simulate data matrix conforming to a mixed Rasch model with two classes
data <- sim.mrm(1000, 20, c(.5, .5))

#Parameter estimation
fit <- mrm(data$data.matrix, 2)

print(fit)

}


