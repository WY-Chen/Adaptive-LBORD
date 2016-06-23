LBORD = function(P,alpha){
  # "level based on recent discoveries", Javanmard & Montanari 2015
  # at time i reject if P[i]<=alpha[i]
  # set alpha[i] = beta[i-tau(i)]
  # where tau(i) is the time of the most recent discovery
  # and beta is a sequence with sum_{i=1}^{infty} beta[i] = alpha
  # scheme: beta[i] proportional to i^{-1.5}
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  for(i in 1:n){
    alphas[i] = beta[i - tau]
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      tau = i
    }
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  output  
}

LBORD_adaptive = function(P,alpha, r, m){
  
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  ct= 0
  gap=0
  w= rep(beta[1],n)
  
  for(i in 1:n){
    alphas[i] = beta[i - tau]
    if (ct > r) {        # magic number r=3
      alphas[i] =m*alpha #here we let m=1.45
    }
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      tau = i
      ct = ct+1
    } else {
      ct= 0
      gap = max(gap, i-max(which(tau==1),0))
    }
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  output  
}

LBORD_adaptive = function(P,alpha, r, m){
  
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  ct= 0
  gap=0
  w= rep(beta[1],n)
  
  for(i in 1:n){
    alphas[i] = beta[i - tau]
    if (ct > r) {        # magic number r=3
      alphas[i] =m*alpha #here we let m=1.45
    }
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      tau = i
      ct = ct+1
    } else {
      ct= 0
      gap = max(gap, i-max(which(tau==1),0))
    }
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  output  
}

LBORD_adaptive = function(P,alpha, r, m){
  
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  ct= 0
  gap=0
  w= rep(beta[1],n)
  
  for(i in 1:n){
    alphas[i] = beta[i - tau]
    if (ct > r) {        # magic number r=3
      alphas[i] =m*alpha #here we let m=1.45
    }
    if(P[i]<=alphas[i]){
      discoveries[i] = 1
      tau = i
      ct = ct+1
    } else {
      ct= 0
      gap = max(gap, i-max(which(tau==1),0))
    }
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  output  
}

LBORD_forgive = function(P,alpha){
  n = length(P)
  alphas = rep(0,n)
  beta = (1:n)^(-1.5); beta = beta * alpha / sum(beta)
  discoveries = rep(0,n)
  tau = 0
  ct= 0
  for (i in 1:n) {
    alphas[i]=beta[i-tau]
    if (P[i]<=alphas[i]) {
      discoveries[i]=1
      tau=i
      ct=ct+1
    } else
      if (P[i]<=ct*alpha/2) {
        ct=ct-floor(1/alpha)
        discoveries[i]=1
        alphas[i]=min(ct*alpha,1)
      } else {
        ct=0
      }
  }
  output = list()
  output$discoveries = discoveries
  output$alphas = alphas
  output 
}

run_methods = function(n,pi0,mu,alpha,r,m,Pvals){
  LBORD_result = LBORD(Pvals,alpha)
  LBORD_AD_result = LBORD_adaptive(Pvals,alpha, r, m)
  LBORD_Forg_result=LBORD_forgive(Pvals,alpha)
  outLBORD = LBORD_result$discoveries
  outLBORD_AD = LBORD_AD_result$discoveries
  outLBORD_Forg=LBORD_Forg_result$discoveries
  #toggle plot option
  if (1==0) {
    plot(Pvals,pch=20,xlim = c(0,n*1.1))
    points(LBORD_AD_result$alphas,col='blue',type='l')
    points(LBORD_result$alphas,col='red',type='l')
    points(LBORD_Forg_result$alphas,col='green',type='l')
    legend('topright',legend=c('LBORD','LBORD_adaptive','LBORD_Forgive'),
           fill=c('red','blue','green'),cex = 0.5)
  }
  list(outLBORD = outLBORD, outLBORD_AD =outLBORD_AD, 
       outLBORD_Forg=outLBORD_Forg)
}

