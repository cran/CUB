#' @title Variance-covariance matrix of a CUB model without covariates
#' @description Compute the variance-covariance matrix of parameter estimates of a CUB model without covariates.
#' @aliases varcovgecub
#' @usage varcovgecub(ordinal,Y,W,X,bet,gama,omega,shelter)
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Y Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param X Matrix of selected covariates to explain the shelter component (default: no covariate is included 
#' in the model)
#' @param bet Parameter vector for the Uncertainty component
#' @param gama Parameter vector for the Feeling component
#' @param omega Parameter vector for the shelter component
#' @param shelter Cateogry corresponding to the shelter effect
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.
#' @seealso \code{\link{probgecub}}
#' @keywords internal



varcovgecub<-function(ordinal,Y,W,X,bet,gama,omega,shelter){
  
  Y<-as.matrix(Y);W<-as.matrix(W);  X<-as.matrix(X);
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(X)==1){
    X<-as.numeric(X)
  }
  
  probi<-probgecub(factor(ordinal,ordered=TRUE),Y,W,X,bet,gama,omega,shelter);
  vvi<-1/probi; dicotom<-ifelse(as.numeric(ordinal)==shelter,1,0)   

  #   Y=as.matrix(Y)
  #   W=as.matrix(W)
  #   X=as.matrix(X)
  paii<-logis(Y,bet);
  csii<-logis(W,gama); deltai<-logis(X,omega);
  m<-length(levels(factor(ordinal,ordered=TRUE)))
  
  bierrei<-bitgama(m,factor(ordinal,ordered=TRUE),W,gama);
  YY<-cbind(1,Y);
  WW<-cbind(1,W); XX<-cbind(1,X);
  np<-NCOL(YY)+NCOL(WW)+NCOL(XX);
  #  /* vettori (n,1)  utili per deriv...prime*/
  mconi<-m-ordinal-(m-1)*csii;
  
  #  /* deriv.prime/prob_i */
  vettAA<-paii*(1-paii)*(1-deltai)*(bierrei-1/m)*vvi
  AA <- Hadprod(YY,vettAA)  ###   n x p+1
  #AA=YY.*paii.*(1-paii).*(1-deltai).*(bierrei-1/m).*vvi;
  vettBB<-paii*(1-deltai)*mconi*bierrei*vvi;
  BB <- Hadprod(WW,vettBB)
  vettCC<-deltai*(dicotom-probi)*vvi
  CC <- Hadprod(XX,vettCC);
  # CC=XX.*deltai.*(dicotom-probi).*vvi;
  
  # /* utili per deriv.seconde/prob_i*/
  dconi<- paii*(1-paii)*(1-2*paii)*(1-deltai)*(bierrei-1/m)*vvi;
  #dconi=paii.*(1-paii).*(1-2*paii).*(1-deltai).*(bierrei-1/m).*vvi;
  gconi<- paii*(1-deltai)*bierrei*(mconi^2-(m-1)*csii)*(1-csii)*vvi;
  #   gconi=paii.*(1-deltai).*bierrei.*(mconi^2-(m-1)*csii.*(1-csii)).*vvi;
  lconi <- deltai*(1-2*deltai)*(dicotom-probi)*vvi;
  #   lconi=deltai.*(1-2*deltai).*(dicotom-probi).*vvi;
  econi <- paii*(1-paii)*(1-deltai)*(bierrei)*mconi*vvi;
  #   econi=paii.*(1-paii).*(1-deltai).*bierrei.*mconi.*vvi;
  fconi<- (-paii)*(1-paii)*deltai*(1-deltai)*(bierrei-1/m)*vvi;
  #   fconi=-paii.*(1-paii).*deltai.*(1-deltai).*(bierrei-1/m).*vvi;
  hconi<-(-paii)*deltai*(1-deltai)*bierrei*mconi*vvi;#   hconi=-paii.*deltai.*(1-deltai).*bierrei.*mconi.*vvi;
  #/* Observed information matrix (gi??? con il segno meno) */
  
  inf11<-t(AA)%*%AA-t(YY)%*%Hadprod(YY,dconi);   #inf11=AA'AA-YY'(YY.*dconi);          #/* (p+1,p+1) */
  inf22<-t(BB)%*%BB-t(WW)%*%Hadprod(WW,gconi);          #/* (q+1,q+1) */
  inf33<-t(CC)%*%CC-t(XX)%*%Hadprod(XX,lconi);          #/* (s+1,s+1) */
  
  inf21<-t(BB)%*%AA-t(WW)%*%Hadprod(YY,econi);          #/* (q+1,p+1) */
  inf31<-t(CC)%*%AA-t(XX)%*%Hadprod(YY,fconi);          #/* (s+1,p+1) */
  inf32<-t(CC)%*%BB-t(XX)%*%Hadprod(WW,hconi);          #/* (s+1,q+1) */  @...prima era: BB'CC @
  
  inf12<-t(inf21);  inf13<-t(inf31);  inf23<-t(inf32);
  
  
  matinf<-rbind(cbind(inf11,inf12,inf13),cbind(inf21,inf22,inf23),cbind(inf31,inf32,inf33));  # Information matrix 
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=np,ncol=np)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=np,ncol=np)
    } else {
      varmat<-solve(matinf)
    }
  }
  
  # 
  # if(det(matinf)<=0){
  #   notpd=1;
  #   cat("=======================================================================","\n")
  #   cat("Variance-covariance matrix NOT positive definite","\n")
  #   cat("=======================================================================","\n")
  # } else {
  #   notpd=0;
  #   varmat=solve(matinf);
  # }                         
  return(varmat)
  
}
