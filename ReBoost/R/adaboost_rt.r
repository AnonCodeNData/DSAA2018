#' AdaBoost.RT
#'
#' @param form formaula
#' @param data trianing data
#' @param niter boosting ites
#' @param thr dd
#' @param power pw
#'
#' @export
AdaBoost.RT.train <-
  function(form,
           data,
           niter = 100,
           thr = 0.1,
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
      m <- rpart::rpart(form, data[train.ind, ])

      models[[t]] <- m
      f <- predict(m, data)

      are <- abs((f - data[, ind.y]) / data[, ind.y])

      err_t <- sum(weights[are > thr])

      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[are <= thr] <- beta_t
      weights[are > thr] <- 1
      weights <- weights / sum(weights)
    }

    names(models) <- paste0("M",seq_along(models))

    list(models=models,beta=betas)
  }

#' AdaBoost.RT predict method
#'
#' @param models list of models
#' @param newdata newdata
#' @param betas betas
#'
#' @export
AdaBoost.RT.predict <-
  function(models, newdata, betas) {
    num <- 0
    pred.mat <- c()
    for (t in 1:length(models)) {
      preds <- rpart::predict(models[[t]], newdata)
      pred.mat <- cbind(pred.mat, preds)
      num <- num + (log(1 / betas[t]) * preds)
    }

    num / sum(log(1 / betas))
  }

#' SMOTEd AdaBoost.RT
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iters
#' @param thr thr
#' @param power pw
#' @param perc.O pwe
#'
#' @export
SMOTEBoost.RT.train <-
  function(form,
           data,
           niter = 100,
           thr = 0.1,
           power = 2,
           perc.O = 1.5) {
    models <- list()
    betas <- c()
    pred.mat <- c()

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

      are <- abs(as.numeric((f - data[, ind.y]) / data[, ind.y]))
      err_t <- sum(weights[are > thr])

      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[are <= thr] <- beta_t
      weights[are > thr] <- 1
      weights <- weights / sum(weights)
    }

    names(models) <- paste0("M",seq_along(models))

    list(models=models,beta=betas)
  }
