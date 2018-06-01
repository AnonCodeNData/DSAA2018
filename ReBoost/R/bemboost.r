#' BEMBoost
#'
#' @param form formual
#' @param data trianing data
#' @param niter boosting iters
#' @param BEM ??
#'
#' @export
BEMBoost.train <-
  function(form, data, niter = 100, BEM=0.5) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(data)==as.character(form[[2]]))

  n <- nrow(data)
  weights <- rep(1/n,n) #initialize weights

  err_t <- 0
  for (t in 1:niter) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart::rpart(form, data[train.ind,])

    models[[t]] <- m

    f <- predict(m,data)

    ae <- abs(f-data[,ind.y])

    grtBEM <- ae>BEM
    errCnt <- sum(grtBEM)

    if(errCnt==0) break

    upfactor <- n/errCnt
    downfactor <- 1/upfactor

    lwrBEM <- !grtBEM

    weights[grtBEM] <- weights[grtBEM] * upfactor
    weights[lwrBEM] <- weights[lwrBEM] * downfactor

    weights <- weights/sum(weights)
  }
  names(models) <- paste0("M",seq_along(models))

  models
}

#' predict method
#'
#' @param model list of models
#' @param newdata new data
#'
#' @export
BEMBoost.predict <-
  function(model, newdata) {
    preds <- lapply(model, rpart::predict, ts)
    preds <- as.data.frame(preds)

    rowMeans(preds)
  }

#' SMOTEd BEMBoost
#'
#' @param form form
#' @param data train data
#' @param niter boost iters
#' @param BEM ??
#' @param perc.O ??
#'
#' @export
SMOTEBoostBEM.train <-
  function(form,
           data,
           niter = 100,
           BEM = 0.5,
           perc.O = 1.5) {
    models <- list()
    betas <- c()
    pred.mat <- c()

    ind.y <- which(colnames(data) == as.character(form[[2]]))

    n <- nrow(data) #size of train

    weights <- rep(1 / n, n) #initialize weights

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

      f <- predict(m, data) ## isto ta bem?

      ae <- abs(f - data[, ind.y])

      grtBEM <- ae > BEM
      errCnt <- sum(grtBEM)

      if (errCnt == 0)
        break

      upfactor <- n / errCnt
      downfactor <- 1 / upfactor

      lwrBEM <- !grtBEM

      weights[grtBEM] <- weights[grtBEM] * upfactor
      weights[lwrBEM] <- weights[lwrBEM] * downfactor

      weights <- weights / sum(weights)
    }
    names(models) <- paste0("M",seq_along(models))

    models
  }
