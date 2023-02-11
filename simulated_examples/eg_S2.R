library(RobustCalibration)
library(RobustGaSP)
library(nloptr)
library(lhs)

limetal02non <- function(xx)
{
 
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- 30 + 5*x1*sin(5*x1)
  fact2 <- 4 + exp(-5*x2)
  
  y <- (fact1*fact2 - 100) / 6
  return(y)
}



cm_model <- function(x, theta)
{

  x1 <- x[,1]
  x2 <- x[,2]
  y=theta[2]*sin(5*x1)+rep(theta[1],dim(x)[1])

  return(y)
}


neg_profile_lik<-function(param){
  theta=(param[1:p_theta])

  beta=exp(param[(p_theta+1):(p_theta+p_x)])
  
  eta=exp(param[p_theta+p_x+1])
  

  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(1,p_x))
  

  R_tilde=R+eta*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  

  output_tilde=output-cm_model(input,theta)
  
  S_2=t(output_tilde)%*%R_tilde_inv%*%output_tilde
  
  -(-1/2*determinant(R_tilde)$modulus[1]-length(output)/2*log(S_2))
  
}

neg_profile_lik_fixed_sigma_2<-function(param){
  theta=(param[1:p_theta])

  beta=exp(param[(p_theta+1):(p_theta+p_x)])
  
  eta=exp(param[p_theta+p_x+1])
  

  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(1,p_x))
  

  R_tilde=R+eta*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  
  output_tilde=output-cm_model(input,theta)
  
  S_2=t(output_tilde)%*%R_tilde_inv%*%output_tilde
  
  -(-1/2*determinant(R_tilde)$modulus[1]-S_2/(2*sigma_2)-length(output)/2*log(sigma_2))
  
}


neg_profile_lik_SGaSP<-function(param){
  theta=(param[1:p_theta])
  beta=exp(param[(p_theta+1):(p_theta+p_x)])
  
  eta=exp(param[p_theta+p_x+1])
  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(1,4))

  R_z=solve(solve(R)+100/sqrt(n)*diag(length(output)))  ##S-GaSP: tilde_lambda=1/2
  

  
  R_tilde=R_z+eta*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  
  output_tilde=output-cm_model(input,theta)
  
  S_2=t(output_tilde)%*%R_tilde_inv%*%output_tilde
  
  -(-1/2*determinant(R_tilde)$modulus-length(output)/2*log(S_2))
  
}


neg_profile_lik_SGaSP_fixed_sigma_2<-function(param){
  theta=(param[1:p_theta])

  beta=exp(param[(p_theta+1):(p_theta+p_x)])
  
  eta=exp(param[p_theta+p_x+1])
  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(1,4))
  
  R_z=solve(solve(R)+100/sqrt(n)*diag(length(output)))  ##S-GaSP: tilde_lambda=1/2
  
  R_tilde=R_z+eta*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  output_tilde=output-cm_model(input,theta)

  S_2=t(output_tilde)%*%R_tilde_inv%*%output_tilde
  
  -(-1/2*determinant(R_tilde)$modulus[1]-S_2/(2*sigma_2)-n/2*log(sigma_2))
  
}


neg_profile_lik_SGaSP_only_theta<-function(param,beta,eta){
  theta=(param[1:p_theta])
  R=separable_kernel(R0,beta,kernel_type='matern_5_2',alpha=rep(1,4))
  
  R_z=solve(solve(R)+100/sqrt(n)*diag(length(output)))  ##S-GaSP: tilde_lambda=1/2
  
  R_tilde=R_z+eta*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  output_tilde=output-cm_model(input,theta)

  S_2=t(output_tilde)%*%R_tilde_inv%*%output_tilde
  
  -(-1/2*determinant(R_tilde)$modulus[1]-S_2/(2*sigma_2)-n/2*log(sigma_2))
  
}

set.seed(1)

p_x=2
p_theta=2
input=maximinLHS(n=30,k=p_x)

n=dim(input)[1]
output=rep(0,n)
for(i in 1:n){
  output[i]=limetal02non(input[i,])
}
output=output+rnorm(n,mean=0,sd=0.05)
m1=rgasp(input,output,nugget.est=T)


R0=m1@R0



test_funct_gasp=neg_profile_lik

m_gasp_record=try(lbfgs( c(rep(-1,p_theta),rep(-1,p_x+1)), test_funct_gasp),silent = T)

while(is.character(m_gasp_record[1])){
  initial_start=-5+10*runif(p_theta+p_x+1)
  
  m_gasp_record=try(lbfgs(initial_start, test_funct_gasp),silent = T)
  
}

for(i_initial in 1:10){
  print(i_initial)
  set.seed(i_initial)
  
  initial_start=-5+10*runif(p_theta+p_x+1)
  
  m_gasp_try=try(lbfgs(initial_start, test_funct_gasp),silent = T)
  if(!is.character(m_gasp_try[1])){
    if( (-m_gasp_record$value)< (-m_gasp_try$value) ){ 
      m_gasp_record=m_gasp_try
    }
  }
}



m_gasp_record




testing_input=cbind(runif(1000),runif(1000))

r0=as.list(1:p_x)
for(i in 1:p_x){
  r0[[i]]=abs(outer(input[,i], testing_input[,i],'-'))
}

testing_output=rep(0,dim(testing_input)[1])
for(i in 1:dim(testing_input)[1]){
  testing_output[i]=limetal02non(testing_input[i,])
}


theta_gasp=(m_gasp_record$par[1:p_theta])
beta_gasp=exp(m_gasp_record$par[(p_theta+1):(p_theta+p_x)])
eta_gasp=exp(m_gasp_record$par[p_theta+p_x+1])
r=separable_kernel(r0,beta_gasp,kernel_type='matern_5_2',alpha=rep(1,p_x))

R=separable_kernel(R0,beta_gasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
median(R)


R_tilde=R+eta_gasp*diag(length(output))

R_tilde_inv=solve(R_tilde)



output_tilde=output-cm_model(input,theta_gasp)


sigma_2_hat_gasp=t(output_tilde)%*%R_tilde_inv%*%(output_tilde)/length(output)
sigma_2_hat_gasp



cm_plus_mean_prediction_gasp=cm_model(testing_input,theta_gasp)

mean_prediction_gasp=cm_plus_mean_prediction_gasp+t(r)%*%R_tilde_inv%*%output_tilde


mean( (cm_plus_mean_prediction_gasp- testing_output)^2 )
mean( (mean_prediction_gasp-testing_output )^2)
c(theta_gasp,sigma_2_hat_gasp,1/beta_gasp, (sigma_2_hat_gasp*eta_gasp) )


test_funct_sgasp=neg_profile_lik_SGaSP

m_sgasp_record=try(lbfgs( c(rep(-1,p_theta),rep(-1,p_x+1)), test_funct_sgasp),silent = T)

while(is.character(m_sgasp_record[1])){
  initial_start=-5+10*runif(p_theta+p_x+1)
  
  m_sgasp_record=try(lbfgs(initial_start, test_funct_sgasp),silent = T)
  
}


for(i_initial in 1:10){
  print(i_initial)
  set.seed(i_initial)
  
  initial_start=-5+10*runif(p_theta+p_x+1)
  
  m_sgasp_try=try(lbfgs(initial_start, test_funct_sgasp),silent = T)
  if(!is.character(m_sgasp_try[1])){
    if( (-m_sgasp_record$value)< (-m_sgasp_try$value) ){ 
      # print("T")
      m_sgasp_record=m_sgasp_try
    }
  }
}

m_sgasp_record


theta_sgasp=(m_sgasp_record$par[1:p_theta])
beta_sgasp=exp(m_sgasp_record$par[(p_theta+1):(p_theta+p_x)])
eta_sgasp=exp(m_sgasp_record$par[p_theta+p_x+1])
r=separable_kernel(r0,beta_sgasp,kernel_type='matern_5_2',alpha=rep(1,p_x))

R=separable_kernel(R0,beta_sgasp,kernel_type='matern_5_2',alpha=rep(1,p_x))

R_z=solve(solve(R)+100/sqrt(n)*diag(length(output)))  ##S-GaSP: tilde_lambda=1/2

R_tilde=R_z+eta_sgasp*diag(length(output))


R_tilde_inv=solve(R_tilde)


output_tilde=output-cm_model(input,theta_sgasp)

sigma_2_hat_sgasp=t(output_tilde)%*%R_tilde_inv%*%(output_tilde)/length(output)
sigma_2_hat_sgasp


R_tilde_middle=R+sqrt(n)/100*diag(length(output))

r_z_t=t(r)%*%t(diag(length(output))-R%*%solve(R_tilde_middle))

r_z=r-R%*%solve(R+sqrt(n)/100*diag(n))%*%r

###computer model plus mean: fit
cm_plus_mean_prediction_sgasp=cm_model(testing_input,theta_sgasp)
mean_prediction_sgasp=cm_plus_mean_prediction_sgasp+t(r_z)%*%R_tilde_inv%*%output_tilde


mean( (cm_plus_mean_prediction_sgasp- testing_output)^2 )
mean( (mean_prediction_sgasp-testing_output )^2)

c(theta_sgasp,sigma_2_hat_sgasp,1/beta_sgasp, (sigma_2_hat_sgasp*eta_sgasp) )

###sim
n_sim=101
mean_prediction_gasp_record=rep(0,n_sim)
mean_prediction_sgasp_record=rep(0,n_sim)

cm_plus_mean_prediction_gasp_record=rep(0,n_sim)
cm_plus_mean_prediction_sgasp_record=rep(0,n_sim)

record_param_gasp=matrix(0, n_sim,p_x+p_theta+2)
record_param_sgasp=matrix(0, n_sim,p_x+p_theta+2)

record_median_cor_gasp=rep(0,n_sim)
record_median_cor_sgasp=rep(0,n_sim)

  
#gasp
for(i_simulation in 1:n_sim){
  print(i_simulation)
  
  sigma_2=1+2*(i_simulation-1)
  

  test_funct_gasp=neg_profile_lik_fixed_sigma_2
  
  m_gasp_record=try(lbfgs( c(rep(-1,p_theta),rep(-1,p_x+1)), test_funct_gasp),silent = T)
  
  while(is.character(m_gasp_record[1])){
    initial_start=-5+10*runif(p_theta+p_x+1)
    
    m_gasp_record=try(lbfgs(initial_start, test_funct_gasp),silent = T)
    
  }
  
  for(i_initial in 1:10){
    set.seed(i_initial)
    
    initial_start=-5+10*runif(p_theta+p_x+1)
    
    m_gasp_try=try(lbfgs(initial_start, test_funct_gasp),silent = T)
    if(!is.character(m_gasp_try[1])){
      if( (-m_gasp_record$value)< (-m_gasp_try$value) ){ 
        m_gasp_record=m_gasp_try
      }
    }
  }
  

  
  
  theta_gasp=(m_gasp_record$par[1:p_theta])
  beta_gasp=exp(m_gasp_record$par[(p_theta+1):(p_theta+p_x)])
  eta_gasp=exp(m_gasp_record$par[p_theta+p_x+1])
  r=separable_kernel(r0,beta_gasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
  
  R=separable_kernel(R0,beta_gasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
  


  #
  R_tilde=R+eta_gasp*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  

  output_tilde=output-cm_model(input,theta_gasp)

  
  cm_plus_mean_prediction_gasp=cm_model(testing_input,theta_gasp)
  
  mean_prediction_gasp=cm_plus_mean_prediction_gasp+t(r)%*%R_tilde_inv%*%output_tilde
  
  
  mean_prediction_gasp_record[i_simulation]=mean( (mean_prediction_gasp- testing_output)^2 )
  cm_plus_mean_prediction_gasp_record[i_simulation]=mean( (cm_plus_mean_prediction_gasp-testing_output )^2)
  record_param_gasp[i_simulation,]=c(theta_gasp,sigma_2,1/beta_gasp, (sigma_2*eta_gasp) )
  
  
}

for(i_sim in 1:n_sim){
  beta_gasp=1/record_param_gasp[i_sim,(p_theta+2):(p_theta+1+p_x) ]
  R=separable_kernel(R0,beta_gasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
  record_median_cor_gasp[i_sim]=median(R)
}

##sgasp

for(i_simulation in 1:n_sim){
  print(i_simulation)
  
  sigma_2=1+2*(i_simulation-1)
  
  
  

   
  test_funct_sgasp=neg_profile_lik_SGaSP_only_theta
  
  m_sgasp_record=try(lbfgs( c(rep(-1,p_theta)), test_funct_sgasp,beta=1/(record_param_gasp[i_simulation,4:5]),eta=(record_param_gasp[i_simulation,6]/record_param_gasp[i_simulation,3])),silent = T)
  
  while(is.character(m_sgasp_record[1])){
    initial_start=-5+10*runif(p_theta)
    
    m_sgasp_record=try(lbfgs(initial_start, test_funct_sgasp,beta=1/(record_param_gasp[i_simulation,4:5]),eta=(record_param_gasp[i_simulation,6]/record_param_gasp[i_simulation,3])),silent = T)
    
  }
  
  
  for(i_initial in 1:10){
    set.seed(i_initial)
    
    initial_start=-5+10*runif(p_theta)
    
    m_sgasp_try=try(lbfgs(initial_start, test_funct_sgasp,beta=1/(record_param_gasp[i_simulation,4:5]),eta=(record_param_gasp[i_simulation,6]/record_param_gasp[i_simulation,3])),silent = T)
    if(!is.character(m_sgasp_try[1])){
      if( (-m_sgasp_record$value)< (-m_sgasp_try$value) ){
        m_sgasp_record=m_sgasp_try
      }
    }
  }
  
  m_sgasp_record
  
  
  theta_sgasp=(m_sgasp_record$par[1:p_theta])

  beta_sgasp=1/(record_param_gasp[i_simulation,4:5])
  eta_sgasp=(record_param_gasp[i_simulation,6]/record_param_gasp[i_simulation,3])
  
  
  r=separable_kernel(r0,beta_sgasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
  
  R=separable_kernel(R0,beta_sgasp,kernel_type='matern_5_2',alpha=rep(1,p_x))
  

  R_z=solve(solve(R)+100/sqrt(n)*diag(length(output)))  ##S-GaSP: tilde_lambda=1/2
  
  R_tilde=R_z+eta_sgasp*diag(length(output))
  
  R_tilde_inv=solve(R_tilde)
  
  
  output_tilde=output-cm_model(input,theta_sgasp)

  
  R_tilde_middle=R+sqrt(n)/100*diag(length(output))
  
  r_z_t=t(r)%*%t(diag(length(output))-R%*%solve(R_tilde_middle))

  r_z=r-R%*%solve(R+sqrt(n)/100*diag(n))%*%r

  cm_plus_mean_prediction_sgasp=cm_model(testing_input,theta_sgasp)
  mean_prediction_sgasp=cm_plus_mean_prediction_sgasp+t(r_z)%*%R_tilde_inv%*%output_tilde
  
  

  cm_plus_mean_prediction_sgasp_record[i_simulation]=  mean( (cm_plus_mean_prediction_sgasp- testing_output)^2 )
   mean_prediction_sgasp_record[i_simulation]=  mean( (mean_prediction_sgasp-testing_output )^2)
  record_param_sgasp[i_simulation,]=  c(theta_sgasp,sigma_2_hat_sgasp,1/beta_sgasp, (sigma_2_hat_sgasp*eta_sgasp) )

  
}

for(i_sim in 1:n_sim){
  beta_sgasp=1/record_param_sgasp[i_sim,(p_theta+2):(p_theta+1+p_x) ]
  R=separable_kernel(R0,beta_sgasp,kernel_type='matern_5_2',alpha=rep(1,p_x))

  record_median_cor_gasp[i_sim]=median(R)
}

 #par(mfrow=c(1,2))
 #pdf("prediction_with_discrepancy.pdf",height=4,width=6)
 plot(1+seq(0,100,1)*2,(mean_prediction_gasp_record),col='red',pch=15,xlab=expression(tau^2),ylab=expression(MSE[f^M]+delta),ylim=c(0,0.025))
 axis(3, at=c(2,50,100,150,200),labels=round(record_median_cor_gasp[c(1,25,50,75,100)],digit=2))
 mtext(expression(hat(rho)), side=3, line=3, cex.lab=1)
 lines(1+seq(0,100,1)*2,mean_prediction_sgasp_record,col='blue',type='p',pch=16)
 legend("topright", col=c('red','blue'), pch=c(15, 16),
                legend=c('GaSP', 'S-GaSP'))
 #dev.off()
 #pdf("prediction_without_discrepancy.pdf",height=4,width=6)
 plot(1+seq(0,100,1)*2,(cm_plus_mean_prediction_gasp_record),col='red',pch=15,xlab=expression(tau^2),ylab=expression(MSE[f^M]),ylim=c(0,160))
 axis(3, at=c(2,50,100,150,200),labels=round(record_median_cor_gasp[c(1,25,50,75,100)],digit=2))
 mtext(expression(hat(rho)), side=3, line=3, cex.lab=1)
 lines(1+seq(0,100,1)*2,cm_plus_mean_prediction_sgasp_record,col='blue',type='p',pch=16)
 legend("topleft", col=c('red','blue'), pch=c(15, 16),
        legend=c('GaSP', 'S-GaSP'))
 #dev.off()
 
