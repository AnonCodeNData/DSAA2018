#' SMOTEd AdaBoost.RT+
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iterations
#' @param thr thr
#' @param power pw
#' @param sigma sad
#'
#' @export
AdaBoost.RTPlus.train <-
  function(form,
           data,
           niter = 100,
           thr = 0.01,
           power = 2,
           sigma = 0.5) {

    models <- list()
    betas <- c()
    pred.mat <- c()
    train_pred.mat <- c()

    ind.y <- which(colnames(data) == as.character(form[[2]]))

    n <- nrow(data)
    weights <- rep(1 / n, n) #initialize weights
    err_t <- 0

    for (t in 1:niter) {
      train.ind <- sample(1:n, n, replace = TRUE, prob = weights)
      m <- rpart::rpart(form, data[train.ind, ])

      models[[t]] <- m

      f <- predict(m, train)

      are <- abs((f - data[, ind.y]) / data[, ind.y])

      err_t <- sum(weights[are > thr])

      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[are <= thr] <- beta_t
      weights[are > thr] <- 1
      weights <- weights / sum(weights)
    }

    for (t in 1:niter) {
      preds <- predict(models[[t]], data)
      train_pred.mat <- cbind(train_pred.mat, preds)
    }

    delta <- 0

    if (is.null(sigma)) {
      # no regularization
      delta <- t(MASS::ginv(t(train_pred.mat))) %*% data[, ind.y]
    } else {
      delta <-
        t(ginv(train_pred.mat %*% t(train_pred.mat) +
                 sigma * diag(nrow(data))) %*% train_pred.mat) %*% data[, ind.y]
    }

    names(models) <- paste0("M",seq_along(models))

    list(models = models, delta=delta)
  }

#' predict method for AdaBoost.RTPlus
#'
#' @param models lsit of models
#' @param newdata new data
#' @param delta delta
#'
#' @export
AdaBoost.RTPlus.predict <-
  function(models, newdata, delta) {
    pred.mat<-c()
    for (t in 1:length(models)) {
      preds <- rpart::predict(models[[t]], test)
      pred.mat <- cbind(pred.mat, preds)
    }
    pred.mat %*% delta
  }

#' SMOTEd AdaBoost.RT+
#'
#' @param form formula
#' @param data training data
#' @param niter boosting iterations
#' @param thr thr
#' @param power pw
#' @param sigma sad
#' @param perc.O perco
#'
#' @export
SMOTEBoost.RTPlus.train <-
  function(form,
           data,
           niter = 100,
           thr = 0.01,
           power = 2,
           sigma = 0.5,
           perc.O = 1.5) {

    models <- list()
    betas <- c()
    pred.mat <- c()
    train_pred.mat <- c()

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
      are <- abs((f - data[, ind.y]) / data[, ind.y])

      err_t <- sum(weights[are > thr])

      beta_t <- err_t ^ power
      betas[[t]] <- beta_t

      weights[are <= thr] <- beta_t
      weights[are > thr] <- 1
      weights <- weights / sum(weights)
    }

    for (t in 1:niter) {
      preds <- predict(models[[t]], data)
      train_pred.mat <- cbind(train_pred.mat, preds)
    }

    delta <- 0

    if (is.null(sigma)) {
      # no regularization
      delta <- t(MASS::ginv(t(train_pred.mat))) %*% data[, ind.y]
    } else {
      delta <-
        t(MASS::ginv(train_pred.mat %*% t(train_pred.mat) +
                 sigma * diag(nrow(data))) %*% train_pred.mat) %*% data[, ind.y]
    }

    names(models) <- paste0("M",seq_along(models))

    list(models = models, delta=delta)
  }
