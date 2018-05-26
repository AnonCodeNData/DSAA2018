adaSMOTE <- function(form,dat,perc.o,rel.thr,k,pc=NULL) {
  
  require(UBL)
  
  y <- dat[,as.character(form[[2]])]
  
  if(length(pc)!=3) {
    pc <- UBL::phi.control(y = y,method = "extremes",coef=1.5)
  }
  
  pc.out <- pc.m <- matrix(pc$control.pts,nrow = (length(pc$control.pts)/3),ncol = 3,byrow = TRUE)
  
  new.dat <- c()
  
  if(any(pc.m[,2]==1)) {
    
    percs <- list()
    
    if(nrow(pc.out)==3) {
      if(pc.out[1,2]==1) {
        
        if(sum(y<=pc.out[1,1]) >= 1) {
          percs <- c(percs,perc.o)
        } else {
          pc.m <- pc.m[-1,]
        }
      } 
      
      percs <- c(percs,1)
      
      if(pc.out[3,2]==1) {
        if(sum(y>=pc.out[3,1]) >= 1) {
          percs <- c(percs,perc.o)
        } else {
          pc.m <- pc.m[-3,]
        }
      } 
    } else {
      
      if(pc.m[1,2]==1) {
        
        if(sum(y<=pc.m[1,1]) >= 1) {
          percs <- c(percs,perc.o)
        } else {
          percs <- c(percs,1)
        }
        
      } else {
        percs <- c(percs,1)
      }
      
      if(pc.m[2,2]==1) {
        
        if(sum(y>=pc.m[2,1]) >= 1) {
          percs <- c(percs,perc.o)
        } else {
          percs <- c(percs,1)
        }
        
      } else {
        percs <- c(percs,1)
      }
      
    }
    
    if(length(percs)>1) {
      
      #' Snippet: SMOTE requires distinct cases. If there's only one case and its repetitions, 
      #' We randomly select a numerical column (except for the target) and add Gaussian noise (sd=0.001)
      
      if(nrow(pc.out)==3) {
        
        dat.high <- dat[which(y>=pc$control.pts[7]),]
        dat.low <- dat[which(y<=pc$control.pts[1]),]
        tgt <- which(colnames(dat)==as.character(form[[2]]))
        
        num.colname <- names(which(sapply(dat[,-tgt],is.numeric)))
        rnd.col <- as.numeric(which(colnames(dat)==sample(num.colname,1)))
        
        if(nrow(unique(dat.high))==1 & nrow(dat.high)>1) {
          dat[which(y>=pc$control.pts[7]),rnd.col] <- rnorm(n=nrow(dat.high),mean = dat.high[1,rnd.col],sd = 0.001)
          dat[which(y>=pc$control.pts[7]),tgt] <- rnorm(n=nrow(dat.high),mean = dat.high[1,tgt],sd = 0.001)
        }
        
        if(nrow(unique(dat.low))==1 & nrow(dat.low)>1) {
          
          dat[which(y<=pc$control.pts[1]),rnd.col] <- rnorm(n=nrow(dat.low),mean = dat.low[1,rnd.col],sd = 0.001)
          dat[which(y<=pc$control.pts[1]),tgt] <- rnorm(n=nrow(dat.low),mean = dat.low[1,tgt],sd = 0.001)
        }
        
        rm(dat.high,dat.low,rnd.col)
        
      } else {
        
        
        
      }
      
      #' End of snippet
      
      if(any(sapply(dat,is.numeric)==FALSE)){ #If there's any nominal predictor, use HEOM distance
        
        new.dat <- UBL::SmoteRegress(form,dat,rel=pc.m,thr.rel=rel.thr,C.perc=percs,k=k,dist="HEOM")
        
      } else { #If all predictors are numerical, use Euclidean distance
        
        new.dat <- UBL::SmoteRegress(form,dat,rel=pc.m,thr.rel=rel.thr,C.perc=percs,k=k,dist="Euclidean")
        
      }
      
    } else {
      
      warning("Did not found any extreme cases. Returning the original train set.")
      new.dat <- dat
      
    }
    
    
  } else {
    
    new.dat <- dat
    warning("Did not found any extreme cases. Returning the original train set.")
  }
  
  new.dat
  
}

eval.stats <- function(form,train,test,preds,cf=1.5,thr=0.99,beta=1) {
  
  require(UBL)
  require(uba)
  
  trues.y <- resp(form,train)
  trues <- test[,as.character(form[[2]])]
  
  ph <- uba::phi.control(trues.y, method="extremes",coef=cf)
  ls <- uba::loss.control(trues.y)
  
  mu <- util(preds,trues,ph,ls,util.control(umetric="MU",event.thr=thr))
  u_new <- util(preds,trues,ph,ls,util.control(umetric="P",event.thr=thr),return.uv=TRUE)
  
  phi.trues <- UBL::phi(trues,control.parms = ph)
  phi.preds <- UBL::phi(preds,control.parms = ph)
  
  pr_frm <- data.frame(Utility=u_new)
  pr_frm["phiTrues"] <- phi.trues
  pr_frm["phiPreds"] <- phi.preds
  
  rmse= sqrt(mean((trues-preds)^2))
  
  # rmse_phi= sqrt(mean(phi.trues[phi.trues>thr]*(trues[phi.trues>thr]-preds[phi.trues>thr])^2))
  
  prec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>thr,]$phiPreds)
  rec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>thr,]$phiTrues)
  F1.u <- (1+beta) * prec.u * rec.u / ( beta^2 * prec.u + rec.u)
  
  prec0.9 <- sum(1+pr_frm[pr_frm$phiTrues>=0.9 & pr_frm$phiPreds>=0.9,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>=0.9,]$phiPreds)
  rec0.9 <- sum(1+pr_frm[pr_frm$phiTrues>=0.9 & pr_frm$phiPreds>=0.9,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>=0.9,]$phiTrues)
  F1.0.9 <- (1+beta) * prec0.9 * rec0.9 / ( beta^2 * prec0.9 + rec0.9)
  
  prec0.8 <- sum(1+pr_frm[pr_frm$phiTrues>=0.8 & pr_frm$phiPreds>=0.8,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>=0.8,]$phiPreds)
  rec0.8 <- sum(1+pr_frm[pr_frm$phiTrues>=0.8 & pr_frm$phiPreds>=0.8,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>=0.8,]$phiTrues)
  F1.0.8 <- (1+beta) * prec0.8 * rec0.8 / ( beta^2 * prec0.8 + rec0.8)
  
  prec0.7 <- sum(1+pr_frm[pr_frm$phiTrues>=0.7 & pr_frm$phiPreds>=0.7,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>=0.7,]$phiPreds)
  rec0.7 <- sum(1+pr_frm[pr_frm$phiTrues>=0.7 & pr_frm$phiPreds>=0.7,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>=0.7,]$phiTrues)
  F1.0.7 <- (1+beta) * prec0.7 * rec0.7 / ( beta^2 * prec0.7 + rec0.7)
  
  prec0.6 <- sum(1+pr_frm[pr_frm$phiTrues>=0.6 & pr_frm$phiPreds>=0.6,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>=0.6,]$phiPreds)
  rec0.6 <- sum(1+pr_frm[pr_frm$phiTrues>=0.6 & pr_frm$phiPreds>=0.6,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>=0.6,]$phiTrues)
  F1.0.6 <- (1+beta) * prec0.6 * rec0.6 / ( beta^2 * prec0.6 + rec0.6)
  
  prec0.5 <- sum(1+pr_frm[pr_frm$phiTrues>=0.5 & pr_frm$phiPreds>=0.5,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>=0.5,]$phiPreds)
  rec0.5 <- sum(1+pr_frm[pr_frm$phiTrues>=0.5 & pr_frm$phiPreds>=0.5,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>=0.5,]$phiTrues)
  F1.0.5 <- (1+beta) * prec0.5 * rec0.5 / ( beta^2 * prec0.5 + rec0.5)
  
  c(rmse=rmse, mu=mu, prec.u=prec.u, rec.u=rec.u, F1.u=F1.u, F1.0.9=F1.0.9, F1.0.8=F1.0.8, F1.0.7=F1.0.7, F1.0.6=F1.0.6, F1.0.5=F1.0.5)
  
  
}

# bias <- function(mat,trues) {
#   
#   bias.vec <- c()
#   
#   for(i in 1:ncol(mat)) {
#     bias.vec <- c(bias.vec,mean(mat[,i])-mean(trues))
#   }
#   
#   mean(bias.vec)^2
#   
# }

bias <- function(mat,trues) {
  
  mean(apply(mat,2,function(ch) {
    mean(ch) - mean(trues)
  })) ^ 2
  
}


# variance <- function(mat) {
#   
#   var.vec <- c()
#   
#   for(i in 1:ncol(mat)) {
#     var.vec <- c(var.vec,var(mat[,i]))
#   }
#   
#   mean(var.vec) / ncol(mat)
#   
# }

variance <- function(mat) {
  
  mean(apply(mat,2, function(x) {
    mean((x - mean(x))^2)
  })) / ncol(mat)
  
}
  


# covariance <- function(mat) {
#   
#   n <- ncol(mat)
#   ens_cov <- cov(mat)
#   diag(ens_cov) <- 0
#   ens_cov <- sum(ens_cov)
#   ens_cov <- ens_cov / (n * (n-1))
#   ens_cov / ((ncol(mat)-1)/ncol(mat))
#   
# }

covariance <- function(mat) {
  
  M <- ncol(mat)
  ens_cov <- cov(mat)
  diag(ens_cov) <- 0
  ens_cov <- sum(ens_cov) 
  ens_cov <- ens_cov / (M * (M-1))
  ens_cov * (1 - 1/M)
  
}

ambiguity <- function(mat,preds,betas) {
  
  w <- log(1/betas) / sum(log(1/betas))
  
  amb.vec <- c()
  
  for(i in 1:ncol(mat)) {
    amb.vec <- c(amb.vec,w[i]*(mat[,i]-preds)^2)
  }
  
  sum(amb.vec)
  
}
