\name{importData}
\alias{importData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import Data from a Folder
}
\description{
Function to import data (Nifti images, text data) from a folder for processing and analysis in R.
}
\usage{
DATA <- importData(mainFolder, input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mainFolder}{
name of the folder inside the current working directory where the data is located.
}
  \item{input}{
a list containing at least two named input variables (see details below).
}
}
\details{
\code{importData} takes a list of different data types and imports them into R for analysis and processing (for example using machine leaning algorithms). \code{input} should contain at least two named variables (usually x an y). Each one is a list with details about the data being imported.

For \strong{images}, \code{type} should be set to "image", \code{path} should be the name of the folder with the Nifti image files (located inside the \code{mainFolder}) and \code{mask} is the mask file (Analyze or Nifti) used to import the images (could be a whole brain mask or a mask for grey matter or certain brain region - also located inside the \code{mainFolder}).

For \strong{text}, \code{type} should be set to "text" and \code{path} to the name of the text file. If a csv extension is found, text will be imported accordingly (comma as separator, dot as decimal). Otherwise, text will be assumed tab/space delimited. Missing values can be written as a space, lowercase x or NA. Important: the first collumn of the text file should always contain the IDs of subjects.
}
\value{
The function returns a list where items have the same name as the input list and contain the imported data. Images are imported as a n*m matrix with n subjects and m features (voxels). Text is imported as a data.frame. The first collumn contains the subjects IDs and each row represents a subject. The imported data can be accessed using the $name operator.

The function also produces to global variables: \code{maskArray} is the imported mask file (as a x*y*z*1 array) and \code{maskIndex} is a n*m matrix where n=features and m=coordinates of features. This can be used when printing the results back to an image (for example a weight map).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Manuel Schutze
}
\note{
This function needs the AnalyzeFMRI package to run.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create input list
input <- list(
  x=list(type="image", path="img", mask="mask.nii"),
  y=list(type="text", path="sample.tsv")
)

## Import data
DATA <- importData("sample",input)

## Results
DATA$x #prints matrix with the image data
DATA$y #prints table with the text data
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
