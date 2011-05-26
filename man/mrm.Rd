\name{mrm}
\alias{mrm}
\title{cMLE and model selection in binary mixed Rasch models.}
\description{
Fits a binary mixed Rasch model via cML estimation (Rost 1990) and provides information criteria for model selection purposes. The core part of the routine has been written in C++, incorporating parts of the Scythe Statistical Library (2007) for matrix manipulations. Starting values for item parameters are uniformly drawn from the interval [-2, 2] and standardized appropriately (sum = 0). 
Starting values for latent score probabilities are uniformly drawn from the interval [0, 1] and standardized to sum up to one in each class. Rows with missing values, as well as constant rows and columns are excluded. The function returns an object of class 'mrm'.
}

\usage{
mrm(data.matrix, cl, max.it = 1000, conv.crit = .0001)
}

\arguments{
  \item{data.matrix}{Input 0/1 data matrix or data frame, rows representing individuals and columns representing items. Rows with missing values as well as constant rows are excluded from the analysis.}
  \item{cl}{The number of classes to be fitted.}
  \item{max.it}{Maximum number of iterations.}
  \item{conv.crit}{If the absolute difference between two successive log-likelihoods falls below this value, the iteration procedure is terminated.}
}

\value{
  \item{beta}{Item easiness parameters.}
  \item{pi.r.c}{Latent score probabilities.}
  \item{class.size}{Estimated class sizes.}
  \item{logLik}{Conditional log-likelihood.}
  \item{AIC}{AIC.}
  \item{BIC}{BIC.}
  \item{number.of.iterations}{Total number of iterations required.}
  \item{number.of.parameters}{Number of parameters.}
  \item{conv.to.bound}{Either 0 or 1, where 1 indicates termination due to divergence to the boundary of the parameter space, i.e. the modulus of an item parameter exceeds 20 (see Preinerstorfer and Formann 2011 for details).}
}

\author{
David Preinerstorfer

\email{david.preinerstorfer@univie.ac.at}

\url{http://homepage.univie.ac.at/david.preinerstorfer}
}

\references{
  Rost (1990).
  Rasch models in latent classes: An integration of two approaches to item analysis. 
  \emph{Applied Psychological Measurement, 14}, 271-282. 

  Pemstein, D., Quinn, K. M. and Martin, A. D. (2007). 
  Scythe Statistical Library: Application programmer's interface. 
  \url{http://scythe.wustl.edu/api/index.hgml}. 

  Preinerstorfer, D. and Formann, A. K. (2011)
  Parameter recovery and model selection in mixed Rasch models.
  \emph{British Journal of Mathematical and Statistical Psychology,} in press.
}


\examples{

#Simulate data matrix conforming to a mixed Rasch model with two classes
data <- sim.mrm(1000, 20, c(.5, .5))

#Parameter estimation
fit <- mrm(data$data.matrix, 2)

}

