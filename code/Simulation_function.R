simulation.function <- function(nSims,n,K,t.max,betas,phi,ntrial,alpha,nu,gamma,D,cens){
  
  # parameters collection
  beta0<-c()
  beta1<-c()
  
  
  sigmau<-c()
  sigmav<-c()
  Phi <- c()
  
  Alpha<-c()
  gamma0<-c()
  Nu <- c()
  
  

  TS.betas.coef <- c()
  TS.betas.se <- c() 
  
  TS.rand.coef<-c()
  TS.sigmau.coef <- c()
  TS.sigmav.coef <- c()
  TS.sigmau.se <- c()
  TS.sigmav.se <- c()
  
  TS.psi.coef<-c()
  TS.psi.se<-c()
  
  TSalph.coef<-c()
  TSalph.se<-c() 
  significance.TS <- 0
  
  JMno.conv <- TSno.conv <- 0  
  
  out_data <- list()
  
  for(i in 1:nSims){
    
    ###### JM ######
    JMconv <- 1 
    
    while( JMconv == 1 ){
      
      dataset <- Simulation_data(n,K,t.max,betas,phi,ntrial,alpha,nu,gamma,D,cens)
      print("data set simulated")
      
      # GAUS LEGENDRE QUADRATRUE  
      glq <- gauss.quad(15, kind = "legendre")
      xk <- glq$nodes   # nodes
      wk <- glq$weights # weights
      k <- length(xk)   # K-points
      
      dataset_stan <- list(N= nrow(dataset$longitudinal) ,
                           n=nrow(dataset$survival), 
                           m = ntrial,
                           y=dataset$longitudinal$y,
                           times=dataset$longitudinal$time, 
                           ID=as.numeric(dataset$longitudinal$id),
                           Time=dataset$survival$Time,
                           status=dataset$survival$event, K= k, xk= xk, wk= wk)
      print("stan data done")
      
      jm_rstan  <- stan(file = "Stan/JMBB.stan",
                        data = dataset_stan,pars=c("betas","Alpha","Sigma","phi","gamma","nu"),
                        init = 0, 
                        chains = 3,
                        iter = 2000,
                        cores = getOption("mc.cores",3))
     
      
      print("stan model done") 
      
      jm <- summary(jm_rstan)$summary
      
      JMconv <- as.numeric(jm[3,10] > 1.05 || jm[3,9] < 600) #bad rhat or low n_eff
      print(paste("convergence was", JMconv))    
      print(jm)
      if (JMconv==1){
        JMno.conv <- JMno.conv + JMconv
        print(paste("convergence was",JMconv))
        next}
       }
    
    ### TSBB ###
    TSconv <- 1 
    
    while (TSconv == 1 ){
      
    #STEP 1
      
    rZ.i <- model.matrix(~id-1,data=dataset$longitudinal)
    rZ.t<-model.matrix(~id:time-1,data = dataset$longitudinal)
    rZ<-cbind(rZ.i,rZ.t)
   
    BB<- try(PROreg::BBmm(fixed.formula = y~time,
                          Z=rZ,
                          nRandComp =c(n,n),
                          m=ntrial,data=dataset$longitudinal,
                          maxiter = 100,		
                          show = TRUE))
    sum.BB<-try(summary(BB))
    
    
    if (class(BB)[1]=="try-error"){
      TSconv <- 1
      TSno.conv <- TSno.conv + TSconv
      print(paste("convergence was",TSconv))
      next} else if (BB$conv=="no"){
      TSconv <- 1
      TSno.conv <- TSno.conv + TSconv
      print(paste("convergence was",TSconv))
      next} else if (class(sum.BB)[1]=="try-error"){
      TSconv <- 1
      TSno.conv <- TSno.conv + TSconv
      print(paste("convergence was",TSconv))
      next } else{
      TSconv <-0
      TSno.conv <- TSno.conv + TSconv
      print(paste("convergence was",TSconv))}
    
    #STEP 2
    X. <- cbind(1,dataset$survival$Time)
    Z..i <- model.matrix(~id-1, data = dataset$survival)
    Z..t <- model.matrix(~id:Time-1, data = dataset$survival)
    Z.. <- cbind(Z..i,Z..t)
    etaY <- try(X.%*%BB$fixed.coef+Z..%*%BB$random.coef)
    exY <- exp(etaY)/(1+exp(etaY))
    
    dataset$survival$exY <- exY
    TS <- try(survival::coxph(survival::Surv(Time,event)~exY, data = dataset$survival))
   
    sum.TS <- summary(TS)
    
    if(is.na(TS$coefficients)){
      TSconv <- 1
      TSno.conv <- TSno.conv + TSconv
      print(paste("convergence was",TSconv))
      next}
    }
    print("TSBB done")
    print(TS)

   print(paste("iter number",i))
    

    
    ######################
    # RESULTS COLLECTION #
    ######################
    
    # JM # 
    
    beta0<-rbind(beta0,jm[1,])
    beta1<-rbind(beta1,jm[2,])
    
    sigmau<-rbind(sigmau,jm[4,])
    sigmav<-rbind(sigmav,jm[7,])
    Phi<-rbind(Phi,jm[8,])
    
    
    Alpha<-rbind(Alpha,jm[3,])
    gamma0<-rbind(gamma0,jm[9,])
    Nu<-rbind(Nu,jm[10,])
    
    

    # TS #
    
    TS.betas.coef <- cbind(TS.betas.coef,BB$fixed.coef)
    TS.betas.se <- cbind(TS.betas.se,sqrt(diag(BB$fixed.vcov))) 
    
    TS.rand.coef<-cbind(TS.rand.coef,BB$random.coef)
    TS.sigmau.coef <- cbind(TS.sigmau.coef,BB$sigma.coef[1])
    TS.sigmav.coef <- cbind(TS.sigmav.coef,BB$sigma.coef[2])
    TS.sigmau.se <- cbind(TS.sigmau.se,sum.BB$sigma[1,2])
    TS.sigmav.se <- cbind(TS.sigmav.se,sum.BB$sigma[2,2])
    
    TS.psi.coef<-cbind(TS.psi.coef, BB$psi.coef)
    TS.psi.se<-cbind(TS.psi.se,sum.BB$psi.table[1,2])
    
    TSalph.coef<-c(TSalph.coef,TS$coefficients)
    TSalph.se<-c(TSalph.se,sqrt(diag(TS$var))) 
    if (sum.TS$coefficients[,5]<0.05){
      significance.TS <- significance.TS+1
    }
    
    print("results collected for this iter")
    
    #save the generated data
    out_data[i] <- list(dataset)
  }
  
  
  long_est<-list(beta0=beta0,beta1=beta1,sigmau=sigmau,sigmav=sigmav,Phi=Phi)
  surv_est <- list(Alpha=Alpha,gamma0=gamma0,Nu=Nu)
  
  result_JM <- list(long =long_est, surv=surv_est)
  
  
  
  longitudinal <- list(beta0 = rbind(TS.betas.coef[1,],TS.betas.se[1,]),
                       beta1 = rbind(TS.betas.coef[2,],TS.betas.se[2,]),
                       sigmau = rbind(TS.sigmau.coef,TS.sigmau.se),
                       sigmav = rbind(TS.sigmav.coef,TS.sigmav.se),
                       psi = rbind(TS.psi.coef,TS.psi.se),
                       rand.coef = TS.rand.coef)
  
  survival = list(Alpha = cbind(TSalph.coef,TSalph.se,significance.TS))
  
  result_TS <- list(longitudinal = longitudinal,
                    survival = survival)
  
  
  out <- list(result_JM = result_JM, 
              result_TS = result_TS,
              parameters = list(alpha,phi,betas,D),
              Data = out_data,
              convergence = c(JMno.conv,TSno.conv))
  
}
