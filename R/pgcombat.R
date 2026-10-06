library(R6)

pgCombat = R6Class(
  "pgCombat",
  public = list(
    batches = NULL,
    L = NULL,
    S = NULL,
    gamma.star = NULL,
    delta.star = NULL,
    Xc = NULL,
    X0 = NULL,
    
    fit = function (dat,
                    batch,
                    mod = NULL,
                    par.prior = TRUE,
                    ref.batch = NULL,
                    mean.only = FALSE) {
      if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
      }
      batch <- as.factor(batch)
      batchmod <- model.matrix( ~ -1 + batch)
      if (!is.null(ref.batch)) {
        if (!(ref.batch %in% levels(batch))) {
          stop("reference level ref.batch is not one of the levels of the batch variable")
        }
        #cat("Using batch =", ref.batch, "as a reference batch (this batch won't change)\n")
        ref <-
          which(levels(as.factor(batch)) == ref.batch)
        batchmod[, ref] <- 1
      }
      else {
        ref <- NULL
      }
      message("Found", nlevels(batch), "batches")
      n.batch <- nlevels(batch) # number of batches
      batches <-
        list() # a list with length 1..n.batch with each element the column indexes of the members of corresponding batch in dat
      for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
      }
      n.batches <-
        sapply(batches, length) # number of dat columns (samples) for each batch
      if (!mean.only & any(n.batches == 1)) {
        stop("At least one batch has only 1 observation, consider using L-model")
      }
      n.array <- sum(n.batches)
      design <- cbind(batchmod, mod)
      check <-
        apply(design, 2, function(x)
          all(x == 1))
      if (!is.null(ref)) {
        check[ref] <- FALSE
      }
      
      # check if any IDs have sd = 0
      zero_sd_indices <- which(apply(dat, 1, function(row) sd(row) == 0))
      if (length(zero_sd_indices > 0)){
        stop(paste0("Error: variables (peptides) with 0 standard deviation! Row indices: ", paste(zero_sd_indices, collapse = ", "), 
                    ". Please remove them prior to ComBat."))
      }
      
      
      design <- as.matrix(design[,!check])
      # message("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)")
      if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n.batch + 1)) {
          stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if (ncol(design) > (n.batch + 1)) {
          if ((qr(design[,-c(1:n.batch)])$rank < ncol(design[,-c(1:n.batch)]))) {
            stop(
              "The covariates are confounded! Please remove one or more of the covariates so the design is not confounded"
            )
          }
          else {
            stop(
              "At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat"
            )
          }
        }
      }
      NAs <- any(is.na(dat))
      
      if (NAs) {
        stop("Missing values are not allowed")
      }
      #cat("Standardizing Data across genes\n")
      
      B.hat <-
        solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
      
      
      if (!is.null(ref.batch)) {
        grand.mean <- t(B.hat[ref,])
      }
      else {
        grand.mean <- crossprod(n.batches / n.array, B.hat[1:n.batch,])
      }
      
      if (!is.null(ref.batch)) {
        ref.dat <- dat[, batches[[ref]]]
        var.pooled <-
          ((ref.dat - t(design[batches[[ref]],] %*% B.hat)) ^ 2) %*% rep(1 / n.batches[ref], n.batches[ref])
      }
      else {
        var.pooled <- ((dat - t(design %*% B.hat)) ^ 2) %*%
          rep(1 / n.array, n.array)
      }
      
      stand.mean <-
        t(grand.mean) %*% t(rep(1, n.array))
      if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
      }
      s.data <-
        (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1,
                                                         n.array)))
      #message("Fitting L/S model and finding priors")
      batch.design <- design[, 1:n.batch]
      
      gamma.hat <-
        solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                  as.matrix(s.data)))
      
      
      delta.hat <- NULL
      for (i in batches) {
        delta.hat <-
          rbind(delta.hat, private$rowVars(s.data[, i, drop = FALSE], na.rm = TRUE))
      }
      gamma.bar <- rowMeans(gamma.hat)
      t2 <- private$rowVars(gamma.hat)
      a.prior <-
        apply(delta.hat, 1, private$aprior)
      b.prior <-
        apply(delta.hat, 1, private$bprior)
      
      gamma.star <-
        delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
      if (par.prior) {
        results <- lapply(1:n.batch, function(i) {
          if (mean.only) {
            gamma.star <- private$postmean(gamma.hat[i,], gamma.bar[i],
                                           1, 1, t2[i])
            delta.star <- rep(1, nrow(s.data))
          }
          else {
            temp <- private$it.sol(s.data[, batches[[i]]],
                                   gamma.hat[i,],
                                   delta.hat[i,],
                                   gamma.bar[i],
                                   t2[i],
                                   a.prior[i],
                                   b.prior[i])
            gamma.star <- temp[1,]
            delta.star <- temp[2,]
          }
          list(gamma.star = gamma.star, delta.star = delta.star)
        })
        for (i in 1:n.batch) {
          gamma.star[i,] <- results[[i]]$gamma.star
          delta.star[i,] <- results[[i]]$delta.star
        }
      }
      else {
        #message("Finding nonparametric adjustments")
        results <- lapply(1:n.batch, function(i) {
          temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                             gamma.hat[i,], delta.hat[i,])
          list(gamma.star = temp[1,], delta.star = temp[2,])
        })
        for (i in 1:n.batch) {
          gamma.star[i,] <- results[[i]]$gamma.star
          delta.star[i,] <- results[[i]]$delta.star
        }
      }
      if (!is.null(ref.batch)) {
        gamma.star[ref,] <- 0
        delta.star[ref,] <- 1
      }
      # message("Adjusting the Data\n")
      bayesdata <- s.data
      j <- 1
      for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,] %*% gamma.star)) 
          (sqrt(delta.star[j,]) %*% t(rep(1,
                                          n.batches[j])))
        j <- j + 1
      }
      bayesdata <-
        (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
                                                 n.array)))) + stand.mean
      if (!is.null(ref.batch)) {
        bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
      }
      self$X0 = dat
      self$Xc = bayesdata
      self$batches = batch
      self$L = stand.mean[, 1]
      self$S = var.pooled
      self$gamma.star = gamma.star
      self$delta.star = delta.star
      return(self)
    },
    
    apply = function(dat, batch) {
      if (is.null(self$batches))
        stop("Empty combat model, use pgCombat:fit before using apply")
      batch = as.factor(batch)
      if (length(batch) != dim(dat)[2])
        stop("Data matrix and batch variable don't match.")
      lbx = levels(self$batches)
      Xc = matrix(nrow = dim(dat)[1], ncol = dim(dat)[2])
      for (i in 1:dim(dat)[2]) {
        Xc[, i] = (dat[, i] - self$L) / sqrt(self$S)
        bIdx = (1:length(lbx))[batch[i] == lbx]
        Xc[, i] = (Xc[, i] - self$gamma.star[bIdx, ]) /
          sqrt(self$delta.star[bIdx, ])
        Xc[, i] = Xc[, i] * sqrt(self$S) + self$L
      }
      return(Xc)
    }
  ),
  private = list(
    aprior = function(gamma.hat) {
      m <- mean(gamma.hat)
      s2 <- var(gamma.hat)
      (2 * s2 + m ^ 2) / s2
    },
    
    bprior =  function(gamma.hat) {
      m <- mean(gamma.hat)
      s2 <- var(gamma.hat)
      (m * s2 + m ^ 3) / s2
    },
    
    postmean =  function(g.hat, g.bar, n, d.star, t2) {
      (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
    },
    
    rowVars = function(X, ...) {
      v = apply(X, 1, var, ...)
    },
    
    postvar = function(sum2, n, a, b) {
      (.5 * sum2 + b) / (n / 2 + a - 1)
    },
    
    it.sol  = function(sdat, g.hat, d.hat, g.bar, t2, a, b, conv =
                         .0001) {
      n <- rowSums(!is.na(sdat))
      g.old <- g.hat
      d.old <- d.hat
      change <- 1
      count <- 0
      while (change > conv) {
        g.new <- private$postmean(g.hat, g.bar, n, d.old, t2)
        sum2 <-
          rowSums((sdat - g.new %*% t(rep(1, ncol(
            sdat
          )))) ^ 2, na.rm = TRUE)
        d.new <- private$postvar(sum2, n, a, b)
        change <-
          max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)
        g.old <- g.new
        d.old <- d.new
        count <- count + 1
      }
      ## cat("This batch took", count, "iterations until convergence\n")
      adjust <- rbind(g.new, d.new)
      rownames(adjust) <- c("g.star", "d.star")
      adjust
    }
  )
)