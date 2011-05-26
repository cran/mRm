\name{plot.mrm}
\alias{plot.mrm}
\title{Plot method for Objects of Class mrm.}
\description{
 Figures of item parameters and conditional score probabilities are generated.
}
\usage{
\method{plot}{mrm}(x, \dots)
}
\arguments{
  \item{x}{An object of class 'mrm'.}
  \item{\dots}{Additional parameters to plot.}
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

plot(fit)

}


