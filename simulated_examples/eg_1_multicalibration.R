library(RobustGaSP)
library(RobustCalibration)


##  simulation for plots
set.seed(1)
k=10
beta_measurement=rep(10,k)

beta_model=50

n=100
input_measurement=list()
input_model=matrix(0,n,1)

R0_measurement=list()
R_measurement=list()
L_measurement=list()
delta_measurement=list()

input_real=seq(0,1,1/(n-1))
var_measurement_real=seq(0.4,0.8,(0.8-0.4)/(k-1))

var_model_real=.2^2


delta_measurement_vec=matrix(0,k,n)


cor_real=.5
cov_function_real=matrix(cor_real,k,k)
diag(cov_function_real)=1
chol_cov_function_real=t(chol(cov_function_real))

input_model=input_real
R0_model=abs(outer(input_model,input_model,'-'))

R_model=matern_5_2_funct(R0_model,beta_model)
L_model=t(chol(var_model_real*R_model))

for(i in 1:k){
  input_measurement[[i]]=input_real
  R0_measurement[[i]]=abs(outer(input_measurement[[i]],input_measurement[[i]],'-'))
  R_measurement[[i]]=matern_5_2_funct(R0_measurement[[i]],beta_measurement[i])
  L_measurement[[i]]=t(chol(var_measurement_real[i]*R_measurement[[i]]))
  
}



output=list()
error_sd=rep(0.05,k)



theta_real=c(pi/2)

for(i in 1:k){

  delta_measurement[[i]]= L_measurement[[i]]%*%rnorm(n,0,1)
}



delta_model=L_model%*%rnorm(n,0,1)


math_model_sin<-function(x,theta){
  sin(theta[1]*x)
}

reality=list()

for(i in 1:k){
  reality[[i]]=math_model_sin(input_real,theta_real)+delta_model+delta_measurement[[i]]
  output[[i]]=reality[[i]]+rnorm(n,0,sd=error_sd[i])
}



math_model=list()
for(i in 1:length(input_measurement)){
  math_model[[i]]=math_model_sin
}

##here the true discrepancy (delta) is from GaSP so the GaSP is expected to behave better. Here S-GaSP is misspecified but it behaves well

##1. s-gasp with measurement bias
system.time(
  for(i in 1:1){
    model_sgasp=rcalibration_MS(design=input_measurement,observations=output,p_theta=1,math_model=math_model,
                                simul_type=rep(1, length(input_measurement)),
                                S=20000,S_0=4000,thinning=10,measurement_bias=T,shared_design=input_model,
                                have_measurement_bias_recorded=T,
                                discrepancy_type=c(rep('GaSP',length(input_measurement)),'S-GaSP'),
                                theta_range=matrix(c(-2*pi,2*pi),1,2),sd_proposal_theta=rep(0.02,2))
  }
)
# 

##2. gasp with measurement bias
system.time(
  for(i in 1:1){
    model_gasp=rcalibration_MS(design=input_measurement,observations=output,p_theta=1,math_model=math_model,
                               simul_type=rep(1, length(input_measurement)),
                               S=20000,S_0=4000,thinning=10,measurement_bias=T,shared_design=input_model,
                               have_measurement_bias_recorded=T,
                               discrepancy_type=c(rep('GaSP',length(input_measurement)),'GaSP'),
                               theta_range=matrix(c(-2*pi,2*pi),1,2),sd_proposal_theta=rep(0.02,2))
  }
)

##predict reality
M=dim(model_gasp@post_individual_par[[1]])[1] ###number of samples saved 
pred_reality_gasp_invididual=pred_reality_sgasp_invididual=0

for(i_M in 1 :M ){
  pred_reality_gasp_invididual=pred_reality_gasp_invididual+math_model_sin(input_real,model_gasp@post_theta[i_M,])+model_gasp@post_delta[i_M,]
  pred_reality_sgasp_invididual=pred_reality_sgasp_invididual+math_model_sin(input_real,model_sgasp@post_theta[i_M,])+model_sgasp@post_delta[i_M,]
}
pred_reality_gasp_invididual=pred_reality_gasp_invididual/M
pred_reality_sgasp_invididual=pred_reality_sgasp_invididual/M


output_stack=0

for(i in 1:k){
  reality[[i]]=math_model_sin(input_real,theta_real)+delta_model+delta_measurement[[i]]
  output[[i]]=reality[[i]]+rnorm(n,0,sd=error_sd[i])
  output_stack=output_stack+ output[[i]]
}
output_stack=output_stack/k


#3. stack image, SGaSP

model_sgasp_stack=rcalibration(design=as.matrix(input_real),observations=as.matrix(output_stack),p_theta=1,math_model=math_model_sin,
                               simul_type=1,
                               S=20000,S_0=4000,thinning=10,
                               discrepancy_type=c('S-GaSP'),
                               theta_range=matrix(c(-2*pi,2*pi),1,2),sd_proposal=c(rep(0.02,2),rep(0.25,1),0.25))

Pred_sgasp_stack=predict(model_sgasp_stack,as.matrix(input_real),math_model=math_model_sin)

#4. stack image, GaSP

model_gasp_stack=rcalibration(design=as.matrix(input_real),observations=as.matrix(output_stack),p_theta=1,math_model=math_model_sin,
                              simul_type=1,
                              S=20000,S_0=4000,thinning=10,
                              discrepancy_type=c('GaSP'),
                              theta_range=matrix(c(-2*pi,2*pi),1,2),sd_proposal=c(rep(0.02,2),rep(0.25,1),0.25))

Pred_gasp_stack=predict(model_gasp_stack,as.matrix(input_real),math_model=math_model_sin)


#5. average of data 
average_data=rep(0,length(pred_reality_gasp_invididual))

for(i in 1:k){
  average_data=average_data+output[[i]]
}

average_data=average_data/k

##define truth
truth=delta_model+math_model_sin(input_real,theta_real)

###compute MSE
MSE_record_truth=rep(NA,5)
MSE_record_delta=rep(NA,4)
MSE_record_measurement_bias=matrix(NA,k,4)

MSE_record_truth[1]=(mean( (pred_reality_gasp_invididual-truth)^2))
MSE_record_truth[2]=(mean( (pred_reality_sgasp_invididual-truth)^2))
MSE_record_truth[3]=(mean( (Pred_gasp_stack@mean-truth)^2))
MSE_record_truth[4]=(mean( (Pred_sgasp_stack@mean-truth)^2))
MSE_record_truth[5]=(mean( (average_data-truth)^2))

MSE_record_delta[1]=(mean( (colMeans(model_gasp@post_delta)-delta_model)^2))
MSE_record_delta[2]=(mean( (colMeans(model_sgasp@post_delta)-delta_model)^2))
MSE_record_delta[3]=(mean( (Pred_gasp_stack@mean-Pred_gasp_stack@math_model_mean-delta_model)^2))
MSE_record_delta[4]=(mean( (Pred_sgasp_stack@mean-Pred_sgasp_stack@math_model_mean-delta_model)^2))

for(i_k in 1:k){
  MSE_record_measurement_bias[i_k,1]=(mean( (colMeans(model_gasp@post_measurement_bias[[i_k]])-delta_measurement[[i_k]])^2))
  MSE_record_measurement_bias[i_k,2]=(mean( (colMeans(model_gasp@post_measurement_bias[[i_k]])-delta_measurement[[i_k]])^2))
  MSE_record_measurement_bias[i_k,3]=(mean( (output[[i_k]]-Pred_gasp_stack@mean-delta_measurement[[i_k]])^2))
  MSE_record_measurement_bias[i_k,4]=(mean( (output[[i_k]]-Pred_sgasp_stack@mean-delta_measurement[[i_k]])^2))
}

colMeans(MSE_record_measurement_bias)

###make some plots
i=1 ##plot the first one source

model_sgasp_LB_95_post_measurement_bias_i=rep(NA,n)
model_sgasp_UB_95_post_measurement_bias_i=rep(NA,n)

for(j in 1:n){
  model_sgasp_LB_95_post_measurement_bias_i[j]=quantile(model_sgasp@post_measurement_bias[[i]][,j],c(0.025))
  model_sgasp_UB_95_post_measurement_bias_i[j]=quantile(model_sgasp@post_measurement_bias[[i]][,j],c(0.975))
  
}

input_plot=seq(0,1,1/(n-1))
#pdf("measurement_bias_1_k_10.pdf",height=4,width=6)
plot(input_plot,delta_measurement[[i]],xlab='x',ylab=expression(delta[1](x)),ylim=c(-2,2),
     main='Measurement Bias',mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)
polygon( c(input_plot,rev(input_plot)),c(model_sgasp_LB_95_post_measurement_bias_i,
                                               rev(model_sgasp_UB_95_post_measurement_bias_i)),col = "grey80", border = F)
lines(input_plot,delta_measurement[[i]],xlab='x',ylab=expression(delta[1](x)),ylim=c(-2,2),main='Measurement Bias',mgp=c(2.5,1,0),
      cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)
lines(input_plot,colMeans(model_gasp@post_measurement_bias[[i]]),xlab='x',ylim=c(-2,2),type='l',col='red',lty=2,lwd=2)
lines(input_plot,colMeans(model_sgasp@post_measurement_bias[[i]]),xlab='x',ylim=c(-2,2),type='l',col='blue',lty=3,lwd=2)
lines(input_plot,output[[i]]-Pred_gasp_stack@mean,xlab='x',ylim=c(-2,2),type='l',col='green',lty=4,lwd=2)
leg <- c("Truth",  "GaSP","S-GaSP","GaSP Stack")
legend("bottomright", legend = leg, col = c("black", "red", "blue","green"),
       lty = c(1,2,3,4),cex=.8,lwd=c(1.5,2,2,2))
#dev.off()


i=2
model_sgasp_LB_95_post_measurement_bias_i=rep(NA,n)
model_sgasp_UB_95_post_measurement_bias_i=rep(NA,n)

for(j in 1:n){
  model_sgasp_LB_95_post_measurement_bias_i[j]=quantile(model_sgasp@post_measurement_bias[[i]][,j],c(0.025))
  model_sgasp_UB_95_post_measurement_bias_i[j]=quantile(model_sgasp@post_measurement_bias[[i]][,j],c(0.975))
  
}

#pdf("measurement_bias_2_k_10.pdf",height=4,width=6)
plot(input_plot,delta_measurement[[i]],xlab='x',ylab=expression(delta[2](x)),ylim=c(-2,2),
     main='Measurement Bias',mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)
polygon( c(input_plot,rev(input_plot)),c(model_sgasp_LB_95_post_measurement_bias_i,
                                         rev(model_sgasp_UB_95_post_measurement_bias_i)),col = "grey80", border = F)
lines(input_plot,delta_measurement[[i]],xlab='x',ylab=expression(delta[1](x)),ylim=c(-2,2),main='Measurement Bias',mgp=c(2.5,1,0),
      cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)
lines(input_plot,colMeans(model_gasp@post_measurement_bias[[i]]),xlab='x',ylim=c(-2,2),type='l',col='red',lty=2,lwd=2)
lines(input_plot,colMeans(model_sgasp@post_measurement_bias[[i]]),xlab='x',ylim=c(-2,2),type='l',col='blue',lty=3,lwd=2)
lines(input_plot,output[[i]]-Pred_gasp_stack@mean,xlab='x',ylim=c(-2,2),type='l',col='green',lty=4,lwd=2)
leg <- c("Truth",  "GaSP","S-GaSP","GaSP Stack")
legend("bottomright", legend = leg, col = c("black", "red", "blue","green"),
       lty = c(1,2,3,4),cex=.8,lwd=c(1.5,2,2,2))
#dev.off()



model_sgasp_LB_95_post_delta=rep(NA,n)
model_sgasp_UB_95_post_delta=rep(NA,n)

for(j in 1:n){
  model_sgasp_LB_95_post_delta[j]=quantile(model_sgasp@post_delta[,j],c(0.025))
  model_sgasp_UB_95_post_delta[j]=quantile(model_sgasp@post_delta[,j],c(0.975))
  
}

#pdf("model_discrepancy_k_10.pdf",height=4,width=6)
plot(input_plot,delta_model,xlab='x',ylab=expression(delta(x)),ylim=c(-2,2),main='Model Discrepancy',
     mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)
polygon( c(input_plot,rev(input_plot)),c(model_sgasp_LB_95_post_delta,
                                         rev(model_sgasp_UB_95_post_delta)),col = "grey80", border = F)
lines(input_plot,delta_model,xlab='x',ylab=expression(delta(x)),ylim=c(-2,2),main='Model Discrepancy',
     mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lwd=1.5)

lines(input_plot,colMeans(model_gasp@post_delta),xlab='x',ylim=c(-2,2),type='l',col='red',lty=2,lwd=2)
lines(input_plot,colMeans(model_sgasp@post_delta),xlab='x',ylim=c(-2,2),type='l',col='blue',lty=3,lwd=2)
lines(input_plot,Pred_gasp_stack@mean-Pred_gasp_stack@math_model_mean,xlab='x',ylim=c(-2,2),type='l',col='green',lty=4,lwd=2)
leg <- c("Truth",  "GaSP","S-GaSP","GaSP Stack")
legend("bottomright", legend = leg, col = c("black", "red", "blue","green"),
       lty = c(1,2,3,4),cex=.8,lwd=c(1.5,2,2,2))
#dev.off()


model_sgasp_LB_95_post_reality=rep(NA,n)
model_sgasp_UB_95_post_reality=rep(NA,n)

pred_reality_sgasp_invididual_all=matrix(NA,M,n)
for(i_M in 1:M){
  pred_reality_sgasp_invididual_all[i_M,]=math_model_sin(input_real,model_sgasp@post_theta[i_M,])+model_sgasp@post_delta[i_M,]
}

for(j in 1:n){
  model_sgasp_LB_95_post_reality[j]=quantile(pred_reality_sgasp_invididual_all[,j],c(0.025))
  model_sgasp_UB_95_post_reality[j]=quantile(pred_reality_sgasp_invididual_all[,j],c(0.975))
  
}

#pdf("reality_k_10.pdf",height=4,width=6)
plot(input_plot,delta_model+math_model_sin(input_real,theta_real),
     xlab='x',ylab=expression({y^R}(x)),mgp=c(2.5,1,0),ylim=c(-2.3,2.3),main='Reality',
     mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lty=1,lwd=1.5)
polygon( c(input_plot,rev(input_plot)),c(model_sgasp_LB_95_post_reality,
                                         rev(model_sgasp_UB_95_post_reality)),col = "grey80", border = F)
lines(input_plot,delta_model+math_model_sin(input_real,theta_real),
     xlab='x',ylab=expression({y^R}(x)),mgp=c(2.5,1,0),ylim=c(-2.3,2.3),main='Reality',
     mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=1.5,type='l',lty=1,lwd=1.5)

lines(input_plot,output[[1]],xlab='x',ylim=c(-2,2),type='p',col='black',pch=17,cex=.6)
lines(input_plot,output[[2]],xlab='x',ylim=c(-2,2),type='p',col='black',pch=19,cex=.6)
lines(input_plot,pred_reality_gasp_invididual,xlab='x',ylim=c(-2,2),type='l',col='red',lty=2,lwd=2)
lines(input_plot,pred_reality_sgasp_invididual,xlab='x',ylim=c(-2,2),type='l',col='blue',lty=3,lwd=2)
lines(input_plot,Pred_gasp_stack@mean,xlab='x',ylim=c(-2,2),type='l',col='green',lty=4,lwd=2)
leg <- c("Truth","Obs Source 1","Obs Source 2",  "GaSP","S-GaSP","GaSP Stack")
legend("bottomright", legend = leg, col = c("black","black","black", "red", "blue","green"),
       pch = c(NA,17,19,NA,NA,NA),lty=c(1,NA,NA,2,3,4),lwd=c(1.5,NA,NA,2,2,2), cex=.8)
#dev.off()


