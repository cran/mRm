\name{sim.mrm}
\alias{sim.mrm}
\title{Generating data matrices conforming to a mixed Rasch model}
\description{
This function generates data matrices conforming to a mixed Rasch model (Rost 1990). Both, person and item parameters may be provided by the user. Otherwise, person parameters are randomly drawn from a standard normal distribution; random equidistant partitions of the interval [-2, 2] are used as item parameters. Class membership of each object is based on a realization of a multinomial random variable with sample size and class proportions as parameters (see Preinerstorfer and Formann 2011 for details).
}

\usage{
sim.mrm(N.sample, N.items, cl.prob, item.para = NULL, 
pers.para = NULL, seed = NULL) 
}

\arguments{
  \item{N.sample}{Sample size.}
  \item{N.items}{Number of items.}
  \item{cl.prob}{Vector of relative class sizes.}
  \item{item.para}{Matrix of item (easiness) parameters. Rows indicate items, columns indicate classes. If no parameters are provided by the user, random permutations of an equidistant partition of the interval [-2, 2] are used in each class.}
  \item{pers.para}{Vector of person parameters. If no parameters are provided by the user, person parameters are drawn from a standard normal distribution.}
  \item{seed}{Seed value.}
}

\value{
  \item{data.matrix}{0/1 data matrix of item responses.}
  \item{beta}{Generated/Provided easiness parameters.}
  \item{emp.probs}{Observed class sizes}
  \item{xi}{Generated/Provided person parameters.}
}

\author{
David Preinerstorfer

\email{david.preinerstorfer@univie.ac.at}

\url{http://homepage.univie.ac.at/david.preinerstorfer}.
}


\references{
  Preinerstorfer, D. and Formann, A. K. (2012)
  Parameter recovery and model selection in mixed Rasch models.
  \emph{British Journal of Mathematical and Statistical Psychology, 65}, 251-262. 

  Rost (1990).
  Rasch models in latent classes: An integration of two approaches to item analysis. 
  \emph{Applied Psychological Measurement, 14}, 271-282. 
}


\examples{

#Simulate a data matrix conforming to a 2-class 
#mixed Rasch model with sample size 1000 and 20 items.

data <- sim.mrm(1000, 20, c(.5, .5))

}

