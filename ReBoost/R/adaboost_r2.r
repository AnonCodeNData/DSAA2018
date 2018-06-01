#' AdaBoost.R2
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iters
#' @param power power
#'
#' @export
AdaBoost.R2.train <- function(form,data,niter=100,power=2) {

  #require(spatstat)

  models <- list()
  betas <- c()

  ind.y <- which(colnames(data)==as.character(form[[2]]))

  n <- nrow(data)
  weights <- rep(1/n,n)
  err_t <- 0

  for (i in 1:niter) {
    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart::rpart(form,data[train.ind,])

    models[[i]] <- m

    f <- predict(m,data)
    ar <- abs(f-data[,ind.y])
    ar <- (ar/max(ar))^power

    err_t <- sum(weights*ar)

    if(err_t>=0.5) break

    beta_t <- err_t / (1-err_t)
    betas[[i]] <- beta_t

    weights <- weights*(beta_t^(1-err_t))
    weights <- weights/sum(weights)
  }
  names(models) <- paste0("M",seq_along(models))

  list(models=models,beta=betas)
}

#' AdaBoost.R2 predict method
#'
#' @param models list of models
#' @param newdata newdata
#'
#' @export
AdaBoost.R2.predict <-
  function(models, newdata) {
    preds <- lapply(models$models, rpart::predict, newdata)
    preds <- as.data.frame(preds)

    apply(preds, 1,
          function(o) {
            stats::weighted.mean(o, models$beta)
          })
  }

#' SMOTEd AdaBoost.R2
#'
#' @param form formula
#' @param data training data
#' @param niter no iterations
#' @param power power?
#' @param perc.O per
#'
#' @export
SMOTEBoost.R2.train <-
  function(form,
           data,
           niter = 100,
           power = 2,
           perc.O = 1.5) {

    models <- list()
    betas <- c()

    ind.y <- which(colnames(data) == as.character(form[[2]]))

    n <- nrow(data)
    weights <- rep(1 / n, n)
    err_t <- 0

    pc <-
      UBL::phi.control(y = data[, ind.y],
                       method = "extremes",
                       coef = 1.5)

    for (t in 1:niter) {
      train.ind <- sample(1:n, n, replace = TRUE, prob = weights)
      new.train <- data[train.ind, ]

      new.train <- adaSMOTE(form, new.train, perc.O, 0.9, k = 3, pc)

      m <- rpart::rpart(form, new.train)

      models[[t]] <- m

      f <- predict(m, data)

      ar <- abs(f - data[, ind.y])
      ar <- (ar / max(ar)) ^ power

      err_t <- sum(weights * ar)

      if (err_t >= 0.5)
        break

      beta_t <- err_t / (1 - err_t)
      betas[[t]] <- beta_t

      weights <- weights * (beta_t ^ (1 - err_t))
      weights <- weights / sum(weights)
    }

    names(models) <- paste0("M",seq_along(models))

    list(models=models,beta=betas)
  }
