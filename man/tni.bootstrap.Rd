\name{tni.bootstrap}
\alias{tni.bootstrap}

\title{
Inference of consensus transcriptional networks.
}

\description{
This function takes a TNI object and returns the consensus transcriptional network.
}

\usage{
tni.bootstrap(object, estimator="pearson", nBootstraps=100, consensus=95, parChunks=10, verbose=TRUE)
}

\arguments{

\item{object}{
a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the method 
\code{\link[RTN:tni.permutation]{tni.permutation}}.
}
\item{estimator}{
a character string indicating which estimator to be used for mutual information computation. 
One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
}
\item{nBootstraps}{
a single integer or numeric value specifying the number of bootstraps 
for deriving a consensus between every TF-target association inferred 
in the mutual information analysis. If running in parallel, nBootstraps 
should be greater and multiple of parChunks.
}
\item{consensus}{
a single integer or numeric value specifying the consensus fraction (in percentage) 
under which a TF-target association is accepted.
}
\item{parChunks}{
an optional single integer value specifying the number of bootstrap chunks to 
be used in the parallel analysis.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}

}

\value{
a matrix in the slot "results" containing a reference transcriptional network, 
see 'tn.ref' option in \code{\link[RTN:tni.get]{tni.get}}.
}

\author{Mauro Castro}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(dt4rtn)

# select 5 regulatoryElements for a quick demonstration!
tfs4test <- dt4rtn$tfs[c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")]

\dontrun{

# preprocessing
rtni <- tni.constructor(expData=dt4rtn$gexp, regulatoryElements=tfs4test, 
        rowAnnotation=dt4rtn$gexpIDs)

# linear version!
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)

# parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(3, "SOCK"))
rtni<-tni.permutation(rtni)
rtni<-tni.bootstrap(rtni)
stopCluster(getOption("cluster"))
}
}


\keyword{methods}

