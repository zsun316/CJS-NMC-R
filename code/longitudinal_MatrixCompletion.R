library(softImpute)
library(missMDA)
library(jomo)
library(mice)
library(denoiseR)
library(h2o)
library(matrixcalc)
library(MASS)
library(lmmlasso)
library(rje)
library(CVTuningCov)
library(Matrix)  
library(stats)
library(Deriv)
library(pracma)
library(matrixStats)
library(fBasics)
library(doParallel)
library(npmr)
library(mice)

#Dim_Nn <- 1000
#Dim_Rr <- 2;Dim_Pp <- 5;Dim_Qq <- 5
#Dim_Tt <- 100

#Mat_LL <- matrix(rnorm(Dim_Nn*Dim_Rr),Dim_Nn,Dim_Rr) %*% 
#   matrix(rnorm(Dim_Rr*Dim_Tt),Dim_Rr,Dim_Tt)

#Mat_HH <- matrix(rnorm(Dim_Pp*Dim_Qq,1,0.1),Dim_Pp,Dim_Qq)

#zero_sampleH <- sample(seq(Dim_Pp*Dim_Qq),
#                       round(Dim_Pp*Dim_Qq * 0.7),replace=FALSE)
#Mat_HHsparse <- Mat_HH
#Mat_HHsparse[zero_sampleH] <- 0

mean_impute <- function(TrueData){
   n <- dim(TrueData$YYmiss)[1]
   p <- dim(TrueData$YYmiss)[2]
   
   omega <- TrueData$OmegaInd
   ImpValue <- sum(TrueData$YYmiss,na.rm = T)/sum(omega)
   
   ImpMatY <- TrueData$YYmiss
   ImpMatY[which(omega == 0)] <- ImpValue
   
   rmse.mean = norm(ImpMatY - TrueData$YYcomp,"F")/sqrt(n*p)
   
   test.error.mean <- norm((1-omega)* (ImpMatY - TrueData$YYcomp),"F")^2/norm((1-omega)* TrueData$YYcomp,"F")^2
   
   return(list(rmse = rmse.mean, test.error = test.error.mean, ImputedY = ImpMatY))
}



cv.softimpute <- function(TrueData,
                          N = 10,
                          thresh = 1e-5,
                          maxit = 100,
                          parallel = T,
                          len = 10,
                          trace.it=T){
   y <- TrueData$YYmiss
   Y2 <- y
   Y2[is.na(Y2)] <- 0
   d <- dim(y)
   n <- d[1]
   p <- d[2]
   m <- sum(!is.na(y))
   na_func <- function(x, prob = 0.1) {
      x <- as.matrix(x)
      omega <- !is.na(x)
      obs.idx <- which(omega)
      yp <- x
      yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
      return(yp)
   }
   
   lambda1.max <- max(svd(Y2)$d)
   lambda1.min <- 1e-4*lambda1.max
   grid.lambda1 <-
      exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
   # Parallel Computing
   NumofCore <- detectCores()
   Cluster_make <- makeCluster(NumofCore)
   registerDoParallel(Cluster_make)
   res.cv <-
      foreach(k = 1:N,
              .packages = c("softImpute", "parallel")) %dopar% (
                 sapply(1:len,
                        function(i) {
                           #yy <- na_func(as.matrix(y), prob = 0.1)
                           yy <- as.matrix(y)
                           if (trace.it){print(paste("lambda", i))}
                           res <-
                              softImpute(as.matrix(yy),
                                         lambda = grid.lambda1[i],
                                         #rank.max = 2,
                                         maxit = 1000)
                           u <- res$u
                           d <- res$d
                           v <- res$v
                           if (is.null(dim(u))) {
                              res <- d * u %*% t(v)
                           } else {
                              res <- u %*% diag(d) %*% t(v)
                           }
                           imp <- as.matrix(yy)
                           imp[is.na(yy)] <- res[is.na(yy)]
                           return(sqrt(sum((res - y) ^ 2, na.rm = T))  )
                        })
              )
   
   res.cv <- colMeans(do.call(rbind, res.cv))
   min_res_cv <- which.min(res.cv)
   lambda <- grid.lambda1[min_res_cv]
   
   BestImp <- softImpute(as.matrix(y),
                         lambda = grid.lambda1[min_res_cv],
                         #rank.max = 2,
                         maxit = 500)
   
   uB <- BestImp$u; dB <- BestImp$d; vB <- BestImp$v 
   if (is.null(dim(uB))) {
      res.best <- as.matrix(dB * uB %*% t(vB))
   }else {
      res.best <- as.matrix(uB %*% diag(dB) %*% t(vB))
   }
   
   dat <-
      data.frame(training.errors = res.cv,
                 lambda1 = grid.lambda1)
   
   rmse.SoftImp = norm(res.best - TrueData$YYcomp,"F")/sqrt(n*p)
   
   test.error.SoftImp <- norm((1-TrueData$OmegaInd)* (res.best - TrueData$YYcomp),"F")^2/norm((1-TrueData$OmegaInd)* TrueData$YYcomp,"F")^2
   
   return(list(
      lambda = lambda,
      errors = dat,
      ImputedMat = res.best,
      rmse = rmse.SoftImp,
      test.error = test.error.SoftImp 
   ))
}




Prob_R1_given_XZ <- function(X_data,Z_data,gamma_r){
   nn <- dim(X_data)[1]; pp <- dim(X_data)[2]
   tt <- dim(Z_data)[1]; qq <- dim(Z_data)[2]
   
   logit_pi <- gamma_r[1] + 
      as.numeric(kronecker(as.numeric(X_data %*% gamma_r[2:(pp+1)]),as.numeric(rep(1,tt))  ) )+
      as.numeric(kronecker(as.numeric(rep(1,nn)), as.numeric(Z_data %*% gamma_r[(pp+2):(pp+qq+1)])  ) ) 
   
   prob_pi <- exp(logit_pi)/(1+exp(logit_pi))
   return(prob_pi)
}



Data_generate <- function(Dim_N = 100, Dim_T = 30, 
                          Dim_P = 5, Dim_Q = 5, Dim_R = 2,
                          missfracY = 0.3,sparsefracH = 0.7,
                          meanX = 0.5,varX = 0.1, 
                          meanZ = 1.5, varZ = 0.1,
                          meanH = 1, varH = 0.1,
                          #Mat_H,Mat_L,
                          SNRatio = NULL,TF_noise = 1,missrate = 70){
   
   Total_NT <- Dim_N * Dim_T 
   
   vec_alpha <- rnorm(Dim_P,0,1)
   vec_beta <- rnorm(Dim_Q,0,1)
   
   #gamma_r <- round(  ( (-(Dim_P + Dim_Q +1)/2):  (Dim_P + Dim_Q+1) )/(Dim_P + Dim_Q + 1), 3) [1:(Dim_P + Dim_Q + 1)]
   #gamma_r <- seq(-0.70,0.5,0.1)[1:(Dim_P + Dim_Q + 1)] ###72%missing
   if(missrate == 70){
      gamma_r <- seq(-0.70,0.5,0.1)[1:(Dim_P + Dim_Q + 1)] ###72%missing
   }else if(missrate == 50){
      gamma_r <- seq(-0.61,0.5,0.1)[1:(Dim_P + Dim_Q + 1)]  ###50%
   }else{
      gamma_r <- seq(-0.54,0.5,0.1)[1:(Dim_P + Dim_Q + 1)]  ###30%
   }
   
   #miss_sampleY <- sample(seq(Total_NT),round(Total_NT*missfracY),replace=FALSE)
   #zero_sampleH <- sample(seq(Dim_P*Dim_Q),round(Dim_P*Dim_Q*sparsefracH),replace=FALSE)
   
   Mat_L <- matrix(rnorm(Dim_N*Dim_R),Dim_N,Dim_R) %*% matrix(rnorm(Dim_R*Dim_T),Dim_R,Dim_T)
   Mat_H <- matrix(rnorm(Dim_P*Dim_Q,meanH,varH),Dim_P,Dim_Q)
   Mat_Hsparse <- Mat_H
   #Mat_Hsparse[zero_sampleH] <- 0
   
   Mat_X <- matrix(rnorm(Dim_N*Dim_P,meanX,varX),Dim_N,Dim_P)
   
   Mat_Z <- matrix(rnorm(Dim_T*Dim_Q,meanZ,varZ),Dim_T,Dim_Q)
   
   miss_prob <- Prob_R1_given_XZ(Mat_X,Mat_Z,gamma_r)
   
   #miss_prob <- missfracY * rep(1,length(miss_prob))
   
   miss_sampleY <- rbinom(Total_NT,1,miss_prob)
   Omega_Y <- matrix(miss_sampleY,Dim_N,Dim_T,byrow = TRUE)
   miss_prob_mat <- matrix(miss_prob,Dim_N,Dim_T,byrow = TRUE)
   
   #print(summary(miss_prob_mat))
   
   
   Mat_X_tilde <- cbind(Mat_X,diag(1,Dim_N))
   Mat_Z_tilde <- cbind(Mat_Z,diag(1,Dim_T))
   Mat_HHsparse_tilde_up <- cbind(Mat_Hsparse,as.vector(vec_alpha))
   Mat_HHsparse_tilde_down <- cbind(t(as.vector(vec_beta)),0)
   
   AA <- rbind(cbind(diag(1,Dim_P),matrix(0,Dim_P,1)),
               cbind(matrix(0,Dim_N,Dim_P),matrix(1,Dim_N,1) ) )
   
   BB <- rbind( cbind(diag(1,Dim_Q),matrix(0,Dim_Q,Dim_T)  ),
                cbind(matrix(0,1,Dim_Q),matrix(1,1,Dim_T) )   )
   
   XAA  <- Mat_X_tilde %*% AA
   ZBB <- Mat_Z_tilde %*% t(BB)
   
   Mat_HHsparse_tilde <- rbind(Mat_HHsparse_tilde_up,Mat_HHsparse_tilde_down)
   
   zero_sampleHtilde <- sample(seq((Dim_P + 1)*(Dim_Q+1)-1),
                               round(((Dim_P+1)*(Dim_Q+1)-1)*sparsefracH),
                               replace=FALSE)
   
   Mat_HHsparse_tilde[zero_sampleHtilde] <- 0
   
   #Mat_HHsparse_tilde <- sim_Mat_HHsparse_tilde2
   
   Mat_Hsparse <- Mat_HHsparse_tilde[1:Dim_P,1:Dim_Q]
   
   #print(dim(XAA))
   #print(dim(ZBB))
   #print(dim(Mat_HHsparse_tilde))
   
   Mat_Yf <- XAA %*% Mat_HHsparse_tilde %*% t(ZBB)
   
   Mat_Yfmiss <- Mat_Yf +  Mat_L
   #Mat_Yfmiss[miss_sampleY] <- NA
   
   Mat_Yfmiss[which(Omega_Y == 0)] <- NA 
   
   if(is.null(SNRatio)){SNRatio <- 1}
   SNR_Mat_Y <- var(as.numeric(Mat_Yf))
   
   if(TF_noise == 1){
      noise <- matrix(rnorm(Dim_N*Dim_T,0,SNR_Mat_Y/SNRatio),Dim_N,Dim_T)
         mvrnorm(n = Dim_N, rep(0,Dim_T), 0.05 * AR1(Dim_T,rho = 0.5))
   }else{
      noise <- 0
   }
   
   
   Mat_Y <- Mat_Yf +  Mat_L + noise
   Mat_Ycomplete <- Mat_Y
   Mat_Yna <- Mat_Ycomplete
   Mat_Yna[which(Omega_Y == 0)] <- NA 
   #Mat_Yna[miss_sampleY] <- NA
   
   
   #Omega_Y <- matrix(1,Dim_N,Dim_T)
   #Omega_Y[miss_sampleY] <- 0
   
   
   return(list(YYcomp = Mat_Ycomplete, YYmiss = Mat_Yna, XX = Mat_X,ZZ = Mat_Z, 
               HH = Mat_Hsparse, HHtilde = Mat_HHsparse_tilde, LL = Mat_L, OmegaInd = Omega_Y, 
               noise = noise, XHZLmiss = Mat_Yfmiss,Mat_XAA = XAA, Mat_ZBB = ZBB,
               alpha = vec_alpha,beta = vec_beta,gamma = gamma_r,
               Miss_Prob_Mat = miss_prob_mat))
}


mean_impute <- function(TrueData){
   n <- dim(TrueData$YYmiss)[1]
   p <- dim(TrueData$YYmiss)[2]
   
   omega <- TrueData$OmegaInd
   ImpValue <- sum(TrueData$YYmiss,na.rm = T)/sum(omega)
   
   ImpMatY <- TrueData$YYmiss
   ImpMatY[which(omega == 0)] <- ImpValue
   
   rmse.mean = norm(ImpMatY - TrueData$YYcomp,"F")/sqrt(n*p)
   
   test.error.mean <- norm((1-omega)* (ImpMatY - TrueData$YYcomp),"F")^2/norm((1-omega)* TrueData$YYcomp,"F")^2
   
   return(list(rmse = rmse.mean, test.error = test.error.mean, ImputedY = ImpMatY))
}


cv.softimpute <- function(TrueData,
                          N = 10,
                          thresh = 1e-5,
                          maxit = 100,
                          parallel = T,
                          len = 20,
                          trace.it=T){
   y <- TrueData$YYmiss
   Y2 <- y
   Y2[is.na(Y2)] <- 0
   d <- dim(y)
   n <- d[1]
   p <- d[2]
   m <- sum(!is.na(y))
   na_func <- function(x, prob = 0.1) {
      x <- as.matrix(x)
      omega <- !is.na(x)
      obs.idx <- which(omega)
      yp <- x
      yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
      return(yp)
   }
   
   lambda1.max <- max(svd(Y2)$d)
   lambda1.min <- 1e-4*lambda1.max
   grid.lambda1 <-
      exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
   # Parallel Computing
   NumofCore <- detectCores()
   Cluster_make <- makeCluster(NumofCore)
   registerDoParallel(Cluster_make)
   res.cv <-
      foreach(k = 1:N,
              .packages = c("softImpute", "parallel")) %dopar% (
                 sapply(1:len,
                        function(i) {
                           #yy <- na_func(as.matrix(y), prob = 0.1)
                           yy <- as.matrix(y)
                           if (trace.it){print(paste("lambda", i))}
                           res <-
                              softImpute(as.matrix(yy),
                                         lambda = grid.lambda1[i],
                                         rank.max = 2,
                                         maxit = 500)
                           u <- res$u
                           d <- res$d
                           v <- res$v
                           if (is.null(dim(u))) {
                              res <- d * u %*% t(v)
                           } else {
                              res <- u %*% diag(d) %*% t(v)
                           }
                           imp <- as.matrix(yy)
                           imp[is.na(yy)] <- res[is.na(yy)]
                           return(sqrt(sum((res - y) ^ 2, na.rm = T))  )
                        })
              )
   
   res.cv <- colMeans(do.call(rbind, res.cv))
   min_res_cv <- which.min(res.cv)
   lambda <- grid.lambda1[min_res_cv]
   
   BestImp <- softImpute(as.matrix(y),
                         lambda = grid.lambda1[min_res_cv],
                         rank.max = 2,
                         maxit = 500)
   
   uB <- BestImp$u; dB <- BestImp$d; vB <- BestImp$v 
   if (is.null(dim(uB))) {
      res.best <- as.matrix(dB * uB %*% t(vB))
   }else {
      res.best <- as.matrix(uB %*% diag(dB) %*% t(vB))
   }
   
   dat <-
      data.frame(training.errors = res.cv,
                 lambda1 = grid.lambda1)
   
   rmse.SoftImp = norm(res.best - TrueData$YYcomp,"F")/sqrt(n*p)
   
   test.error.SoftImp <- norm((1-TrueData$OmegaInd)* (res.best - TrueData$YYcomp),"F")^2/norm((1-TrueData$OmegaInd)* TrueData$YYcomp,"F")^2
   
   return(list(
      lambda = lambda,
      errors = dat,
      ImputedMat = res.best,
      rmse = rmse.SoftImp,
      test.error = test.error.SoftImp 
   ))
}


quad_approx <- function(yy,xx,zz,hh,ll,dh,dl,lambda_L,lambda_H){
   
   #st_gamma <- 0.5
   
   dim_y <- dim(yy);dim_xx <- dim(xx); dim_zz <- dim(zz)
   nn <- dim_y[1]; tt <- dim_y[2]
   #pp <- dim_xx[2]; qq <- dim_zz[2]
   omega <- matrix(0,nn,tt)
   
   omega[which(!is.na(yy) == TRUE)] <- 1
   
   yy[which(is.na(yy) == TRUE)] <- 0
   
   quad <- -2/(nn*tt) * sum(omega * (yy - (xx%*%dh%*%t(zz)) ) * (xx%*%dh%*%t(zz) + dl)  ) + 1/(nn*tt) * sum(omega * ( (xx%*%dh%*%t(zz)) + dl )^2)
   return(quad)
}

quad_approx_h <- function(yy,xx,zz,hh,ll,dh,dl,lambda_L,lambda_H,omega){
   
   gamma_h <- 0.5
   
   dim_y <- dim(yy);dim_xx <- dim(xx); dim_zz <- dim(zz)
   nn <- dim_y[1]; tt <- dim_y[2]
   #pp <- dim_xx[2]; qq <- dim_zz[2]
   
   yy[which(is.na(yy) == TRUE)] <- 0
   
   quad <- -2/(nn*tt) * sum(omega * (yy - (xx%*%hh%*%t(zz) + ll) ) * (xx%*%dh%*%t(zz))  ) 
   + gamma_h * sum(omega * ( xx%*%dh%*%t(zz) )^2) 
   + lambda_H * (sum(hh + dh) - sum(hh) )
   return(quad)
}

quad_approx_l <- function(yy,xx,zz,hh,ll,dh,dl,lambda_L,lambda_H,omega){
   
   gamma_l <- 0.5
   
   yy[which(is.na(yy) == TRUE)] <- 0
   
   quad <- -2 * sum(omega * (yy - (xx%*%hh%*%t(zz) + ll) ) * (dl)  ) +   gamma_l * sum(omega * dl^2) + lambda_L * (nuclear(ll+dl) - nuclear(ll))
   return(quad)
}



Armijo_func <- function(yy,xx,zz,hh,ll,dh,dl,lambda_H,lambda_L,is_alpha_H){
   ###
   # hh is old H
   ###
   
   beta <- 1
   zeta <- 0.5
   nok = TRUE
   
   dim_y <- dim(yy);dim_xx <- dim(xx); dim_zz <- dim(zz)
   nn <- dim_y[1]; tt <- dim_y[2]
   omega <- matrix(0,nn,tt)
   
   omega[which(!is.na(yy) == TRUE)] <- 1
   
   if(is_alpha_H == TRUE){
      dl <- matrix(0,dim(yy)[1],dim(yy)[2])
      quad <- quad_approx_h(yy,xx,zz,hh,ll,dh,dl,lambda_L,lambda_H,omega)
      diff <- sum(omega * (yy - xx%*%(hh + beta *dh)%*%t(zz) - ll )^2 ) - sum(  omega * (yy - xx%*%hh%*%t(zz) - ll )^2 )
      while(nok == TRUE){
         if(diff < beta * zeta * quad){
            nok = FALSE
         }else{
            beta <- beta/2
         }
         if(beta < 1e-3){nok = FALSE}
      }
      #print(paste("The step size for H chosen in the current step is ", beta))
      return(hh + dh*beta)
      
   }else{
      
      dh <- matrix(0,dim(xx)[2],dim(zz)[2]) 
      quad <- quad_approx_l(yy,xx,zz,hh,ll,dh,dl,lambda_L,lambda_H,omega)
      diff <-  sum(omega * (yy - xx%*%(hh)%*%t(zz) - (ll+beta * dl))^2 ) - sum(  omega * (yy - xx%*%hh%*%t(zz) - ll )^2 )
      while(nok == TRUE){
         if(diff < beta * zeta * quad){
            nok = FALSE
         }else{
            beta <- beta/2
         }
         if(beta < 1e-8){nok = FALSE}
      }
      #print(paste("The step size for L chosen in the current step is ", beta))
      return(ll+dl*beta)
   }
}


Target_function_value <- function(Y_mat,XI_mat,ZI_mat,H_tilde_mat,
                                  L_mat,lambda_h,lambda_l,
                                  dh,dl,lr,is_update_h = TRUE,prob_mat){
   dim_y <- dim(Y_mat);dim_xx <- dim(Y_mat); dim_zz <- dim(Y_mat)
   nn <- dim_y[1]; tt <- dim_y[2]
   
   omega <- matrix(0,nn,tt)
   omega[which(!is.na(Y_mat) == TRUE)] <- 1
   
   if(is_update_h == TRUE){
      no_penalty <- sum( (( omega/prob_mat) * (Y_mat - XI_mat %*% (H_tilde_mat + lr * dh)%*% t(ZI_mat) - L_mat) )^2 )
      no_penalty_original <- sum( (( omega/prob_mat) * (Y_mat - XI_mat %*% (H_tilde_mat)%*% t(ZI_mat) - L_mat) )^2 )
      
      penalty_term1 <- lambda_h * sum(abs(H_tilde_mat + lr * dh))
      penalty_term1_original <- lambda_h * sum(abs(H_tilde_mat))
      
      penalty_term2 <- lambda_l * nuclear(L_mat)
      penalty_term2_original <- penalty_term2 
      
   }else{
      no_penalty <- sum(( omega/prob_mat)* (Y_mat - XI_mat %*% (H_tilde_mat)%*% t(ZI_mat) - (L_mat + lr * dl))^2 )
      no_penalty_original <- sum(( omega/prob_mat)* (Y_mat - XI_mat %*% (H_tilde_mat)%*% t(ZI_mat) - L_mat)^2 )
      
      penalty_term1 <- lambda_h * sum(abs(H_tilde_mat))
      penalty_term1_original <- penalty_term1 
      
      penalty_term2 <- lambda_l * nuclear(L_mat + lr * dl)
      penalty_term2_original <- lambda_l * nuclear(L_mat)
   }
   
   value_3_term <- no_penalty + penalty_term1 + penalty_term2
   value_3_term_orginal <- no_penalty_original + penalty_term1_original  + penalty_term2_original
   descent_value <- value_3_term - value_3_term_orginal
   
   
   return(c(value_3_term,value_3_term_orginal,descent_value))
}

Check_LR_GRID <- function(Y_mat,XI_mat,ZI_mat,H_tilde_mat,
                          L_mat,lambda_h,lambda_l,dh,dl,is_update_h = TRUE,
                          len = 40,lr.min = 1e-6,lr.max = 0.1,prob_mat){
   
   lr_vec <- seq(lr.min, lr.max,length.out = len)
   target_value_list  <- c()
   
   for(i in 1: len){
      target_value_list[i] <- Target_function_value(Y_mat,XI_mat,ZI_mat,H_tilde_mat,
                                                    L_mat,lambda_h,lambda_l,
                                                    dh,dl,lr_vec[i],is_update_h,prob_mat)[3]
   }
   #print(target_value_list)
   best_lr <- lr_vec[which(target_value_list == min(target_value_list))]
   
   #print(target_value_list)
   #print(best_lr)
   
   return(best_lr)
   
}

iter.main <- function(TrueData,thresh_H = 1e-5,thresh_L=1e-5,
                      max.iter = 100,len = 20,warm_start_H = NULL, warm_start_L = NULL){
   yy <- TrueData$YYmiss; xx <- TrueData$Mat_XAA; zz <- TrueData$Mat_ZBB 
   probmat <- TrueData$Miss_Prob_Mat
   
   dim_y <- dim(yy);dim_xx <- dim(xx); dim_zz <- dim(zz)
   nn <- dim_y[1]; tt <- dim_y[2]
   pp <- dim_xx[2] - 1; qq <- dim_zz[2] - 1
   
   #H0 <- matrix(rnorm(pp*qq),pp,qq)
   #L0 <- matrix(rnorm(nn*tt),nn,tt)
   fitx_m <- kronecker(zz,xx)
   
   omega <- matrix(0,nn,tt)
   omega[which(!is.na(yy) == TRUE)] <- 1
   
   
   if(is.null(warm_start_H)){
      
      fitx_m2 <- as.matrix(kronecker(zz,xx))
      yymiss <- matrix(c(yy),dim_y[1]*dim_y[2],1)
      
      y_initial_comp <- cv.softimpute(TrueData,
                                      N = 5,
                                      thresh = 1e-5,
                                      maxit = 100,
                                      parallel = T,
                                      len = 20,
                                      trace.it=T)
      
      y_initial <- y_initial_comp$ImputedMat
      yymeanimp <- y_initial
      
      #H0 <- matrix(rnorm((pp+1)*(qq+1)),pp+1,qq+1)
      #yymeanimp <- mean_impute(TrueData)$ImputedY
      
      #data_int <- data.frame(yymiss,fitx_m2)
      
      #print(c(row(data_int1),col(data_int1)))
      #imp_int <- mice(data_int,m = 1,maxit = 3)
      #data_int_comp <- complete(imp_int)
      
      #yy_int_impt <- data_int_comp$yymiss
      #yymeanimp <- as.numeric(yy_int_impt)
      #print(yymeanimp[1:100])
      
      #IntY.fit <- lm(c(yymeanimp) ~ fitx_m2 - 1)
      
      IntY.fit <- glmnet::cv.glmnet(fitx_m2[,-ncol(fitx_m2)],
                                    c(yymeanimp),
                                    intercept = F,
                                    lambda = seq(1e-12,1e-4,length.out = 100),
                                    type.measure = "mse",
                                    nfold = 10,
                                    weights = c(omega)/(c(probmat) * nn*tt),
                                    alpha = 1)
      
      
      #vec_H0 <- as.numeric(IntY.fit$coef)
      vec_H0 <- as.numeric(coef(IntY.fit,"lambda.min"))[-1]
      
      vec_H0[36] <- 0
      H0 <- matrix(vec_H0,(pp+1),(qq+1),byrow = F)
      L0_temp <- yymeanimp - xx%*% H0 %*% t(zz)
      
      
      ### 
      
   }else{
      H0 <- warm_start_H
   }
   
   if(is.null(warm_start_L)){
      #L0 <- matrix(rnorm(nn*tt),nn,tt)
      L0 <- L0_temp 
   }else{
      L0 <- warm_start_L
   }
   
   
   yy[which(!is.na(yy) == FALSE)] <- 0
   
   ratio_H <- 1
   ratio_L <- 1
   iter = 0
   Hn <- H0; Ln <- L0; dl <- 0
   
   
   while (  (ratio_H > thresh_H) || (ratio_L > thresh_L) ||(iter > 20)){
      
      #update H
      
      yy_n <- yy
      yy_n[which(omega == 0)] <- (xx%*%Hn%*%t(zz) + Ln)[which(omega == 0)]
      
      fity_m <- yy_n - Ln
      
      ridge1_cv <- glmnet::cv.glmnet(fitx_m[,-ncol(fitx_m)],
                                     c(fity_m),
                                     intercept = F,
                                     #lambda = seq(1e-12,1,length.out = 100),
                                     type.measure = "mse",
                                     weights = c(omega)/(c(probmat) * nn * tt),
                                     nfold = 10,
                                     alpha = 0)
      
      best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
      
      fit.res <- glmnet::cv.glmnet(fitx_m[,-ncol(fitx_m)],
                                   c(fity_m), 
                                   intercept = F,
                                   #lambda = seq(1e-12,1,length.out = 100),
                                   family = "gaussian", 
                                   type.measure = "mse",
                                   weights = c(omega)/(c(probmat) * nn*tt),
                                   alpha = 1,
                                   penalty.factor = 1 / abs(best_ridge_coef)
      )
      
      #fity_m_nomis <- c(fity_m * omega)[which(c(omega) !=  0)]
      #fitx_m_nomis <- c(fitx_m * omega)[which(c(omega) !=  0)]
      #fit.res <- lm(c(fity_m_nomis) ~ fitx_m_nomis[,-ncol(fitx_m_nomis)] - 1 )
      
      
      #fit.res$coef
      
      
      lambda_H <- fit.res$lambda.min
      #print(lambda_H)
      #print(coef(fit.res,"lambda.min"))
      #print(fit.res$lambda)
      #print(paste("The best tuning parameter for H is ", fit.res$lambda.min))
      #print(fit.res)
      
      vec_H_Transpose <- coef(fit.res,"lambda.min")[-1]
      #print(dim(fitx_m[,-1]))
      #print(vec_H_Transpose)
      
      H_new <- matrix(c(vec_H_Transpose,0),(pp+1),(qq+1),byrow = F)
      #print(Hn)
      #print(H_new)
      
      #coef_H <- seq(10);dh_seq <- seq(10)
      
      #for (i in 1:10){
      #    coef_H <- coef(fit.res,s = 0.0100000 + (i-1)* (4-0.01)/9)[-1]
      #    assign(paste("H_new", i, sep=""), 
      #           matrix(coef_H[i],pp,qq,byrow = F) )
      #    
      #}
      
      ###apply armijo rule
      dh <- H_new - Hn
      
      #H_new <- Armijo_func(yy_n,xx,zz,Hn,Ln,dh,0,lambda_H,0,TRUE)
      
      #decay <- 0.9
      Learn_Rate <- Check_LR_GRID(yy_n,xx,zz,Hn,Ln,lambda_H,0,dh,0,is_update_h = TRUE,20,1e-6,0.5,probmat)
      #H_new <- Hn + 0.5* (decay^iter)* dh
      
      H_new <- Hn + Learn_Rate *  dh
      
      #print(H_new)
      
      #ratio_H <- norm(H_new - Hn, "F")/norm(Hn, "F")
      ratio_H <- norm(H_new - Hn, "F")/norm(Hn, "F")
      #ratio_H <- norm(H_new - Data_sim$HHtilde, "F")
      
      
      #print(vec_H_Transpose)
      Hn <- H_new
      
      #Update L
      fitL_m <- yy - xx%*%Hn%*%t(zz)
      fitL_m[which(omega == 0)] <- NA
      
      fitL_m2 <- fitL_m
      fitL_m2[which(omega == 0)] <- 0
      
      lambdal.max <- max(svd(fitL_m2)$d)
      lambdal.min <- 1e-4*lambdal.max
      
      grid.lambda.l <-
         exp(seq(log(lambdal.min), log(lambdal.max), length.out = len))
      
      
      NumofCore <- detectCores()
      Cluster_make <- makeCluster(NumofCore)
      registerDoParallel(Cluster_make)
      res.cv <-
         foreach(k = 1:2,
                 .packages = c("softImpute", "parallel")) %dopar% (
                    sapply(1:len,
                           function(i) {
                              
                              ll <- as.matrix(fitL_m)
                              #print(paste("lambda", i))
                              res.L <- 
                                 softImpute(as.matrix(ll),
                                            lambda = grid.lambda.l[i],
                                            #rank.max = min(nn,tt)-1,maxit = 1000)
                                            rank.max = 5,maxit = 500)
                              u <- res.L$u
                              d <- res.L$d
                              v <- res.L$v
                              if (is.null(dim(u))) {
                                 res <- d * u %*% t(v)
                              } else {
                                 res <- u %*% diag(d) %*% t(v)
                              }
                              imp <- as.matrix(ll)
                              imp[is.na(ll)] <- res[is.na(ll)]
                              return(sqrt(sum((res - ll) ^ 2, na.rm = T)))
                           }
                    )
                 )
      stopCluster(Cluster_make)
      res.cv <- colMeans(do.call(rbind, res.cv))
      min_res_cv <- which.min(res.cv)
      lambda.min <- grid.lambda.l[min_res_cv]
      #print(paste("The best tuning parameter for L is ", lambda.min))
      
      L_new_svd <- softImpute(as.matrix(fitL_m),
                              lambda = lambda.min,
                              #rank.max = min(nn,tt)-1,maxit = 1000)
                              rank.max = 3, maxit = 500)
      
      u <- L_new_svd$u
      d <- L_new_svd$d
      v <- L_new_svd$v
      if (is.null(dim(u))) {
         L_new <- d * u %*% t(v)
      } else {
         L_new <- u %*% diag(d) %*% t(v)
      }
      
      
      dl <- L_new - Ln
      #L_new <- Armijo_func(yy_n,xx,zz,Hn,Ln,0,dl,0,lambda.min,FALSE)
      
      Learn_Rate <- Check_LR_GRID(yy_n,xx,zz,Hn,Ln,0,lambda.min,0,dl,FALSE,20,1e-6,0.5,probmat)
      
      L_new <- Ln + Learn_Rate * dl
      
      #ratio_L <- norm(L_new - TrueData$LL, "F")/norm(Ln, "F")
      ratio_L <- norm(L_new - Ln, "F")/norm(Ln, "F")
      
      
      error_L <- norm(L_new - TrueData$LL, "F")/(nn*tt)
      error_H <- norm(H_new - TrueData$HHtilde, "F")/((pp+1)*(qq+1)  )
      
      test_err_Y <- norm((1-omega) * (xx%*%H_new%*%t(zz) + L_new - TrueData$YYcomp),"F")^2/norm((1-omega) * TrueData$YYcomp,"F")^2
      rmse_Y <- norm(xx%*%H_new%*%t(zz) + L_new - TrueData$YYcomp,"F")/sqrt(nn*tt)
      
      
      #print(c(ratio_H,ratio_L,error_H,error_L,test_err_Y,rmse_Y))
      
      Ln <- L_new
      
      #dat <- data.frame(errors = res.cv, lambda1 = grid.lambda.l)
      iter <- iter + 1
      #print(Ln)
      
   }
   
   return(list(Test.Error = test_err_Y, RMSE = rmse_Y, 
               H = Hn,L = Ln,LambdaH = lambda_H))
}

Data_sim <- Data_generate(Dim_N = 1000, Dim_T = 30, 
                          Dim_P = 5, Dim_Q = 5, Dim_R = 2,
                          missfracY = 0.7,sparsefracH = 0.7,
                          meanX = 0.5,varX = 0.1, meanZ = 1.5,
                          varZ = 0.1,meanH = 1, varH = 0.1,
                          SNRatio = NULL,TF_noise = 1,missrate = 70)


start.time <- Sys.time()
result_sim <- iter.main(Data_sim,1e-3,3e-2,
                        #max.iter = 200,
                        #warm_start_H = Data_sim$HHtilde + 0.5,
                        #warm_start_L = Data_sim$LL,
                        len = 20)
#result_sim222 <- iter.main(Data_sim,1e-5,1e-5,max.iter = 200,len = 20,
#                           warm_start_H = Data_sim$HH, warm_start_L = Data_sim$LL)

end.time <- Sys.time()
computing_time <- end.time - start.time


#result_softimp <- cv.softimpute(Data_sim)

#start.time <- Sys.time()
#result_sim1.1 <- iter.main(Data_sim,1e-5,1e-5,max.iter = 200,
#                           len = 20)
#end.time <- Sys.time()

#computing_time1.1 <- end.time - start.time

#print(Data_sim$HHtilde)
#print(result_sim$H)
#print(as.numeric(c(result_sim1.1$Test.Error,result_sim1.1$RMSE,
#                   computing_time1.1,result_sim1.1$LambdaH)))


print(as.numeric(c(result_sim$Test.Error,result_sim$RMSE,
                   result_softimp$test.error,result_softimp$rmse)))

as.numeric(c(result_sim$Test.Error,result_sim$RMSE,computing_time,
             result_softimp$test.error,result_softimp$rmse))

#as.numeric(c(result_sim1.1$Test.Error,result_sim1.1$RMSE,computing_time1.1))

