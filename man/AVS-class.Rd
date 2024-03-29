\name{AVS-class}
\Rdversion{1.1}
\docType{class}
\alias{AVS-class}
\alias{avs.vse,AVS-method}
\alias{avs.evse,AVS-method}
\alias{avs.pevse,AVS-method}
\alias{avs.get,AVS-method}

\title{Class \code{"AVS"}: an S4 class for variant set enrichment analysis.}

\description{
This S4 class includes a series of methods to do enrichment analyses in Associated Variant Sets (AVSs).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AVS", markers)}.
}

\section{Slots}{
  \describe{
    \item{\code{markers}:}{Object of class \code{"character"}, a data frame,
    a 'BED file' format with rs# markers mapped to the same genome build of the LD source in the RTNdata package.
    }
    \item{\code{validatedMarkers}:}{Object of class \code{"data.frame"}, 
    a data frame with genome positions of the validated markers.
    }  
    \item{\code{variantSet}:}{Object of class \code{"list"}, 
    an associated variant set.
    }
    \item{\code{randomSet}:}{Object of class \code{"list"}, 
    a random associated variant set. 
    }
    \item{\code{para}:}{Object of class \code{"list"},
    a list of parameters for variant set enrichment analysis.
    }
    \item{\code{results}:}{Object of class \code{"list"},
    a list of results (see return values in the AVS methods).
    }
    \item{\code{summary}:}{Object of class \code{"list"},
    a list of summary information for \code{markers}, \code{para}, 
    and \code{results}.    
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the AVS object 
    based on the available methods.
    }    
  }
}

\section{Methods}{
  \describe{
    \item{avs.vse}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.vse]{avs.vse}} } 
    \item{avs.evse}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.evse]{avs.evse}} }
    \item{avs.pevse}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.pevse]{avs.pevse}} }    
    \item{avs.get}{\code{signature(object = "AVS")}: see \code{\link[RTN:avs.get]{avs.get}} }    
	 }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}
}

\examples{

\dontrun{
#This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDHapMapRel27.hg18)
data(bcarisk, package = "RTNdata.LDHapMapRel27.hg18")
avs <- avs.preprocess.LDHapMapRel27.hg18(bcarisk, nrand=100)
}
}

\keyword{classes}
