#' FarmSelect: Factor Adjusted Robust Model Selection
#'
#' This R package implements a consistent model selection strategy for high dimensional sparse regression when the covariate dependence can be reduced through factor models. By separating the latent factors from idiosyncratic components, the problem is transformed from model selection with highly correlated covariates to that with weakly correlated variables. It is appropriate for cases where we have many variables compared to the number of samples. Moreover, it implements a robust procedure to estimate distribution parameters wherever possible, hence being suitable for cases when the underlying distribution deviates from Gaussianity, which is commonly assumed in the literature. See the paper on this method, Fan et al.(2017) <https://arxiv.org/abs/1612.08490>, for detailed description of methods and further references. For detailed information on how to use and install see \url{https://kbose28.github.io/FarmSelect}.
#'
#' @name FarmSelect
#' @docType package
NULL
