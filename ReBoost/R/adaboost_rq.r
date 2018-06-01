#' AdaBoost.RQ
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iters
#' @param power power??
#'
#' @export
AdaBoost.RQ.train <-
  function(form,
           data,
           niter = 100,
           power = 2) {
    models <- list()
    betas <- c()
    pred.mat <- c()

    ind.y <- which(colnames(data) == as.character(form[[2]]))

    n <- nrow(data)
    weights <- rep(1 / n, n)
    err_t <- 0

    for (t in 1:niter) {
      train.ind <- sample(1:n, n, replace = TRUE, prob = weights)
      m <- rpart::rpart(form, data[train.ind,])
      models[[t]] <- m

      f <- predict(m, data)

      abs.err <- abs(as.numeric((f - data[, ind.y])))

      large.err.ind <-
        which(abs.err > grDevices::boxplot.stats(abs.err)$stats[3])

      err_t <- sum(weights[large.err.ind])

      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[] <- beta_t
      weights[large.err.ind] <- 1

      weights <- weights / sum(weights)
    }
    names(models) <- paste0("M",seq_along(models))

    list(models=models,beta=betas)
  }

#' AdaBoost.RQ predict method
#'
#' @param models list of models
#' @param betas betas
#' @param newdata newdata
#'
#' @export
AdaBoost.RQ.predict <-
  function(models, betas, newdata) {
    pred.mat<-c()
    num <- 0
    for (i in 1:length(models)) {
      preds <- rpart::predict(models[[i]], newdata)
      pred.mat <- cbind(pred.mat, preds)
      num <- num + (log(1 / betas[i]) * preds)
    }
    finalpreds <- num / sum(log(1 / betas))
    finalpreds
  }

#' SMOTEd AdaBoost.RQ
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iters
#' @param power pw
#' @param perc.O po
#'
#' @export
SMOTEBoost.RQ.train <-
  function(form,
           data,
           niter = 100,
           power = 2,
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

    for (t in 1:t_final) {
      train.ind <- sample(1:n, n, replace = TRUE, prob = weights)
      new.train <- data[train.ind, ]

      new.train <- adaSMOTE(form, new.train, perc.O, 0.9, k = 3, pc)

      m <- rpart(form, new.train)

      models[[t]] <- m

      f <- predict(m, data)

      abs.err <- abs(as.numeric((f - data[, ind.y])))

      large.err.ind <- which(abs.err > grDevices::boxplot.stats(abs.err)$stats[3])

      err_t <- sum(weights[large.err.ind])
      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[] <- beta_t
      weights[large.err.ind] <- 1

      weights <- weights / sum(weights)

    }
    names(models) <- paste0("M",seq_along(models))

    models

    list(models=models,beta=betas)
}
