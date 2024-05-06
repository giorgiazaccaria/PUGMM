#' Penguins
#'
#' The data set contains five measurements made on 342 penguins which are classified into three species.
#'
#' @docType data
#'
#' @usage data(penguins)
#'
#' @format
#' A data frame with 342 observations and 5 variables, which are described as follows.
#' \describe{
#'   \item{species}{Penguin species (Chinstrap, Adélie, or Gentoo)}
#'   \item{culmen_length_mm}{Culmen length (mm)}
#'   \item{culmen_depth_mm}{Culmen depth (mm)}
#'   \item{flipper_length_mm}{Flipper length (mm)}
#'   \item{body_mass_g}{Body mass (g)}
#' }
#' @details
#' Data were collected and made available by Dr. Kristen Gorman and the Palmer Station, Antarctica LTER, a member of the Long Term Ecological Research Network.
#' The categorical variables 'island' and 'sex' have been removed from the original dataset, as well as the incomplete observations on the five variables reported herein.
#' @keywords Datasets
#'
#' @references Gorman, K.B., Williams T.D., Fraser W.R. (2014). Ecological sexual dimorphism and environmental variability within a community of Antarctic penguins (genus Pygoscelis). \emph{PLoS ONE}, 9(3), e90081.
#'
#' @source Dataset downloaded from Kaggle \href{https://www.kaggle.com/code/parulpandey/penguin-dataset-the-new-iris}{https://www.kaggle.com/code/parulpandey/penguin-dataset-the-new-iris}.

#' @examples
#' data(penguins)
#'
"penguins"
