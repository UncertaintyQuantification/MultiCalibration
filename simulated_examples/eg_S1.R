library(RobustGaSP)
library(mvtnorm)
library(rospca)
library(pracma)


M=100000
n_max=100
n_min=5
beta_real=50
sigma_2_real=1
output=rep(0,n)
MSE_record_avg=rep(0, n_max-n_min+1)
MSE_z_record_avg=rep(0, n_max-n_min+1)
theta_hat_record=rep(0,M)
theta_z_hat_record=rep(0,M)

for( i_n in n_min:n_max){
  print(i_n)
  n=i_n
  input=seq(0,1,1/(n-1))
  
  
  gamma_real=1/beta_real
  R0_00=abs(outer(input,input,'-'))
  #R=matern_5_2_funct(R0_00,beta_real)
  R=pow_exp_funct(R0_00,beta_real,1)
  
  L=t(chol(R))
  
  
  
  X=rep(1,n)
  rho_real=exp(-beta_real*1/(n-1))
  
  set.seed(i_n)
  
  output_matrix=sigma_2_real*L%*%matrix(rnorm(n*M),n,M)
  
  theta_hat_record=(output_matrix[1,]+output_matrix[n,]+(1-rho_real)*colSums(output_matrix[2:(n-1),]) )/(n-(n-2)*rho_real)
  
  MSE_record_avg[i_n-n_min+1]=mean((theta_hat_record-0)^2)


 ###s-gasp
 tilde_lambda=1/2
 B=(R+1/(tilde_lambda)*diag(n) )
 L_B=t(chol(B))
 R_z=R- R%*%backsolve(t(L_B),forwardsolve(L_B,t(R)))
 L_z=t(chol(R_z))

 output_z_matrix=sigma_2_real*L_z%*%matrix(rnorm(n*M),n,M)
 
 R_z_inv=solve(R_z)
 
 theta_z_hat_record=solve(t(X)%*%R_z_inv%*%X)%*%(t(X)%*%R_z_inv%*%output_z_matrix)
 
 MSE_z_record_avg[i_n-n_min+1]=mean((theta_z_hat_record-0)^2)
 
}

#pdf("eg_1_beta_50.pdf",height=4,width=5)
limiting_MSE=2*sigma_2_real*gamma_real/(1+2*gamma_real)
plot(n_min:n_max, MSE_record_avg,ylim=c(0, 0.25),col='red',pch=20,xlab='n',ylab='MSE')
abline(a=limiting_MSE,b=0)
legend("topright", legend=c("GaSP"),
       col=c("red"), pch=c(20), cex=1)
#dev.off()

