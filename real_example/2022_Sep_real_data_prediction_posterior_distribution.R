library(RobustCalibration)
library(plot3D)
library(fields)
library(ggplot2)
library(gridExtra)





data_ascending_1=read.csv(file=paste('real_data/2018_sep_real_data_ascending1.csv',sep=''),header=F)


output_ascending_1=as.matrix(data_ascending_1[,1])/100   ##need to go back to meter
input_ascending_1=as.matrix(data_ascending_1[,2:3])


data_ascending_4=read.csv(file=paste('real_data/2018_sep_real_data_ascending4.csv',sep=''),header=F)

#data_ascending_stack=read.csv(file='real_data/2018_sep_real_data_ascending_stack.csv',header=F)

output_ascending_4=as.matrix(data_ascending_4[,1])/100   ##need to go back to meter
input_ascending_4=as.matrix(data_ascending_4[,2:3])


data_descending_1=read.csv(file=paste('real_data/2018_sep_real_data_descending1.csv',sep=''),header=F)

#data_descending_stack=read.csv(file='real_data/2018_sep_real_data_descending_stack.csv',header=F)
output_descending_1=as.matrix(data_descending_1[,1])/100
input_descending_1=as.matrix(data_descending_1[,2:3])


data_descending_2=read.csv(file=paste('real_data/2018_sep_real_data_descending2.csv',sep=''),header=F)

#data_descending_stack=read.csv(file='real_data/2018_sep_real_data_descending_stack.csv',header=F)
output_descending_2=as.matrix(data_descending_2[,1])/100
input_descending_2=as.matrix(data_descending_2[,2:3])


data_descending_4=read.csv(file=paste('real_data/2018_sep_real_data_descending4.csv',sep=''),header=F)

#data_descending_stack=read.csv(file='real_data/2018_sep_real_data_descending_stack.csv',header=F)
output_descending_4=as.matrix(data_descending_4[,1])/100
input_descending_4=as.matrix(data_descending_4[,2:3])


##

theta_range=matrix(0,5,2)
theta_range[1,]=c(-2000, 3000)  
theta_range[2,]=c(-2000, 5000)  
theta_range[3,]=c(500, 6000)  
theta_range[4,]=c(0, .15)   ##not sure I should time 2 here
#theta_range[4,]=c(0, .15)*2   
theta_range[5,]=c(.25, .33)  

mogihammer_ascending<-function(x,theta){
  Mogihammer(x,theta,2) ##this is ascending
}

mogihammer_descending<-function(x,theta){
  Mogihammer(x,theta,3) ##this is descending
}


math_model=list()
math_model[[1]]=mogihammer_ascending
math_model[[2]]=mogihammer_ascending
math_model[[3]]=mogihammer_descending
math_model[[4]]=mogihammer_descending
math_model[[5]]=mogihammer_descending

design=list()
design[[1]]=as.matrix(input_ascending_1)
design[[2]]=as.matrix(input_ascending_4)
design[[3]]=as.matrix(input_descending_1)
design[[4]]=as.matrix(input_descending_2)
design[[5]]=as.matrix(input_descending_4)


observations=list()
observations[[1]]=as.matrix(output_ascending_1)
observations[[2]]=as.matrix(output_ascending_4)
observations[[3]]=as.matrix(output_descending_1)
observations[[4]]=as.matrix(output_descending_2)
observations[[5]]=as.matrix(output_descending_4)

shared_design=(design[[1]]+design[[2]]+design[[3]]+design[[4]]+design[[5]])/5
n=(length(output_ascending_1))
X=list()
for(i in 1:5){
  X[[i]]=matrix(1,n,1)
}


####some proposal 

sd_proposal_cov_par_here=list()
for(i in 1:(5)){
  sd_proposal_cov_par_here[[i]]=rep(0.2,dim(design[[i]])[2]+1)
}
sd_proposal_cov_par_here[[6]]=rep(0.1,dim(design[[i]])[2])


set.seed(0)


time_model_gasp=system.time(
  for(ii in 1:1){
    model_gasp=rcalibration_MS(design=design,observations = observations,p_theta=5,
                               X=X, have_trend=rep(T,length(design)),
                               S=50000,S_0=10000, thinning=10,simul_type=rep(1,5),math_model=math_model,theta_range=theta_range,
                               sd_proposal_theta=c(rep(0.025,4),0.04),sd_proposal_cov_par=sd_proposal_cov_par_here,
                               measurement_bias=T, shared_design=shared_design,have_measurement_bias_recorded=T,
                               discrepancy_type=c(rep('GaSP',length(design)),'GaSP'))
  }
)


model_gasp@accept_S_theta
model_gasp@accept_S_beta


par(mfrow=c(2,3))
quilt.plot(x=(model_gasp@input[[1]][,1]), y=(model_gasp@input[[1]][,2]),
           z=colMeans(model_gasp@post_delta),nrow = 80, ncol = 80)

for(i in 1:5){
  quilt.plot(x=(model_gasp@input[[1]][,1]), y=(model_gasp@input[[1]][,2]),
             z=colMeans(model_gasp@post_measurement_bias[[i]]),nrow = 80, ncol = 80)
}
par(mfrow=c(2,3))
for(i in 1:5){
  plot(model_gasp@post_theta[,i])
}

time_model_sgasp=system.time(
  for(ii in 1:1){
    model_sgasp=rcalibration_MS(design=design,observations = observations,p_theta=5,
                                X=X, have_trend=rep(T,length(design)),
                                S=50000,S_0=10000, thinning=10,simul_type=rep(1,5),math_model=math_model,theta_range=theta_range,
                                sd_proposal_theta=c(rep(0.025,4),0.04),sd_proposal_cov_par=sd_proposal_cov_par_here,
                                measurement_bias=T, shared_design=shared_design,have_measurement_bias_recorded=T,
                                discrepancy_type=c(rep('GaSP',length(design)),'S-GaSP'))
  }
)

model_sgasp@accept_S_theta
model_sgasp@accept_S_beta

par(mfrow=c(2,3))
quilt.plot(x=(model_sgasp@input[[1]][,1]), y=(model_sgasp@input[[1]][,2]),
           z=colMeans(model_sgasp@post_delta),nrow = 80, ncol = 80)

for(i in 1:5){
  quilt.plot(x=(model_sgasp@input[[1]][,1]), y=(model_sgasp@input[[1]][,2]),
             z=colMeans(model_sgasp@post_measurement_bias[[i]]),nrow = 80, ncol = 80)
}
par(mfrow=c(2,3))
for(i in 1:5){
  plot(model_sgasp@post_theta[,i])
}

plot(colMeans(model_sgasp@post_delta) )


M=dim(model_gasp@post_individual_par[[1]])[1]
###plot some 

#pdf("posterior_calibration_histogram_Oct_2018.pdf",height=5,width=9)
par(mfrow=c(2,3))

index=1
param <- data.frame(Posterior = factor(rep(c("GaSP","S-GaSP"), each=M)), 
                    theta = c(model_gasp@post_theta[,index],model_sgasp@post_theta[,index]))

plot1<-ggplot(param, aes(x=theta, fill=Posterior)) +
  geom_histogram(binwidth=(theta_range[index,2]-theta_range[index,1])/300,alpha=.5, position="identity",show.legend = FALSE)+labs(x = expression(theta[1]))

index=2
param <- data.frame(Posterior = factor(rep(c("GaSP","S-GaSP"), each=M)), 
                    theta = c(model_gasp@post_theta[,index],model_sgasp@post_theta[,index]))

plot2<-ggplot(param, aes(x=theta, fill=Posterior)) +
  geom_histogram(binwidth=(theta_range[index,2]-theta_range[index,1])/300,alpha=.5, position="identity",show.legend = FALSE)+labs(x =expression(theta[2]))

index=3
param <- data.frame(Posterior = factor(rep(c("GaSP","S-GaSP"), each=M)), 
                    theta = c(model_gasp@post_theta[,index],model_sgasp@post_theta[,index]))

plot3<-ggplot(param, aes(x=theta, fill=Posterior)) +
  geom_histogram(binwidth=(theta_range[index,2]-theta_range[index,1])/300,alpha=.5, 
                 position="identity",show.legend = FALSE)+labs(x = expression(theta[3]),show.legend = FALSE)

index=4
param <- data.frame(Posterior = factor(rep(c("GaSP","S-GaSP"), each=M)), 
                    theta = c(model_gasp@post_theta[,index],model_sgasp@post_theta[,index]))

plot4<-ggplot(param, aes(x=theta, fill=Posterior)) +
  geom_histogram(binwidth=(theta_range[index,2]-theta_range[index,1])/300,alpha=.5, position="identity",show.legend = FALSE)+labs(x = expression(theta[4]))

index=5
param <- data.frame(Posterior = factor(rep(c("GaSP","S-GaSP"), each=M)), 
                    theta = c(model_gasp@post_theta[,index],model_sgasp@post_theta[,index]))

plot5<-ggplot(param, aes(x=theta, fill=Posterior)) +
  geom_histogram(binwidth=(theta_range[index,2]-theta_range[index,1])/100,alpha=.5, position="identity",show.legend = FALSE)+labs(x = expression(theta[5]))

grid.arrange(plot1, plot2,plot3,plot4,plot5, ncol=3, nrow = 2)

#dev.off()



##prediction
t3_ascending_obs_matrix1_full=as.matrix(read.csv(file='real_data/t3_ascending_obs_matrix1_full.csv',header=F))
t3_ascending_x_coordinate1_full=(read.csv(file='real_data/t3_ascending_x_coordinate1_full.csv',header=F))[[1]]
t3_ascending_y_coordinate1_full=(read.csv(file='real_data/t3_ascending_y_coordinate1_full.csv',header=F))[[1]]


image2D(t(t3_ascending_obs_matrix1_full))

sum(is.na(t3_ascending_obs_matrix1_full))

t3_ascending_obs_matrix4_full=as.matrix(read.csv(file='real_data/t3_ascending_obs_matrix4_full.csv',header=F))
t3_ascending_x_coordinate4_full=(read.csv(file='real_data/t3_ascending_x_coordinate4_full.csv',header=F))[[1]]
t3_ascending_y_coordinate4_full=(read.csv(file='real_data/t3_ascending_y_coordinate4_full.csv',header=F))[[1]]


image2D(t3_ascending_obs_matrix4_full)

t3_descending_obs_matrix1_full=as.matrix(read.csv(file='real_data/t3_descending_obs_matrix1_full.csv',header=F))
t3_descending_x_coordinate1_full=(read.csv(file='real_data/t3_descending_x_coordinate1_full.csv',header=F))[[1]]
t3_descending_y_coordinate1_full=(read.csv(file='real_data/t3_descending_y_coordinate1_full.csv',header=F))[[1]]


image2D(t(t3_descending_obs_matrix1_full))


t3_descending_obs_matrix2_full=as.matrix(read.csv(file='real_data/t3_descending_obs_matrix2_full.csv',header=F))
t3_descending_x_coordinate2_full=(read.csv(file='real_data/t3_descending_x_coordinate2_full.csv',header=F))[[1]]
t3_descending_y_coordinate2_full=(read.csv(file='real_data/t3_descending_y_coordinate2_full.csv',header=F))[[1]]

image2D(t(t3_descending_obs_matrix2_full))

t3_descending_obs_matrix4_full=as.matrix(read.csv(file='real_data/t3_descending_obs_matrix4_full.csv',header=F))
t3_descending_x_coordinate4_full=(read.csv(file='real_data/t3_descending_x_coordinate4_full.csv',header=F))[[1]]
t3_descending_y_coordinate4_full=(read.csv(file='real_data/t3_descending_y_coordinate4_full.csv',header=F))[[1]]

image2D(t(t3_descending_obs_matrix4_full))
testing_output=list()
testing_output[[1]]=t3_ascending_obs_matrix1_full
testing_output[[2]]=t3_ascending_obs_matrix4_full
testing_output[[3]]=t3_descending_obs_matrix1_full
testing_output[[4]]=t3_descending_obs_matrix2_full
testing_output[[5]]=t3_descending_obs_matrix4_full

testing_input_separable=list()
for(i in 1:6){
  testing_input_separable[[i]]=list()
}
testing_input_separable[[1]][[1]]=t3_ascending_x_coordinate1_full
testing_input_separable[[1]][[2]]=t3_ascending_y_coordinate1_full

testing_input_separable[[2]][[1]]=t3_ascending_x_coordinate4_full
testing_input_separable[[2]][[2]]=t3_ascending_y_coordinate4_full


testing_input_separable[[3]][[1]]=t3_descending_x_coordinate1_full
testing_input_separable[[3]][[2]]=t3_descending_y_coordinate1_full

testing_input_separable[[4]][[1]]=t3_descending_x_coordinate2_full
testing_input_separable[[4]][[2]]=t3_descending_y_coordinate2_full

testing_input_separable[[5]][[1]]=t3_descending_x_coordinate4_full
testing_input_separable[[5]][[2]]=t3_descending_y_coordinate4_full

aa=0
bb=0
for(i in 1:5){
  aa=aa+testing_input_separable[[i]][[1]]
  bb=bb+testing_input_separable[[i]][[2]]
}
testing_input_separable[[6]][[1]]=aa/5
testing_input_separable[[6]][[2]]=bb/5





X_testing=as.list(rep(1,length(design)))
for(i in 1:length(design)){
  X_testing[[i]]=matrix(1,length(testing_input_separable[[i]][[1]])*length(testing_input_separable[[i]][[2]]),1)
}



pred_time_sgasp=system.time(
  for(i_pred in 1:1){
    predict_sgasp=predict_separable_2dim_MS(model_sgasp, testing_input_separable,
                                            X_testing=X_testing,math_model=math_model)
  }
)
pred_time_gasp=system.time(
  for(i_pred in 1:1){
      predict_gasp=predict_separable_2dim_MS(model_gasp, testing_input_separable,
                                       X_testing=X_testing,math_model=math_model)
  }
)


#  

max_value=-Inf
min_value=Inf


for(i in 1:5){
  index_not_na=which(is.na(testing_output[[i]])==F)
  max_value=max(max_value,max(predict_sgasp@mean[[i]][index_not_na]/100, testing_output[[i]][index_not_na]/100))
  min_value=min(min_value,min(predict_sgasp@mean[[i]][index_not_na]/100,testing_output[[i]][index_not_na]/100))

}



index_plot_1=(1: (541/3) )*3
index_plot_2=(1: (469/3) )*3

avg_image=0
for(i in 1:5){
  avg_image=avg_image+testing_output[[i]]
}
avg_image=avg_image/5

###some comparison 
#pdf("real_data_5_images_1.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
i=1
image2D((t(matrix(testing_output[[i]]/100,469,541)))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        zlim=c(min_value,max_value),main=paste('interferogram',i),
        mgp=c(2.5,1,0),
        cex.axis=2,cex.lab=2,cex.main=2,colkey = FALSE)
#colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5)
for(i in 2:3){
  image2D((t(matrix(testing_output[[i]]/100,469,541)))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],
          y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          zlim=c(min_value,max_value),main=paste('interferogram',i),
          mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = FALSE)
}


#dev.off()
#pdf("real_data_5_images_2.pdf",width=10,height=3.8)

par(mfrow=c(1,3))
for(i in 4:5){
  image2D((t(matrix(testing_output[[i]]/100,469,541)))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],
          y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          ,zlim=c(min_value,max_value),main=paste('interferogram',i),cex.axis=2,cex.lab=2,cex.main=2,mgp=c(2.5,1,0),colkey = FALSE) #colkey = list(plot = T, side = 4,width = 0.5)
}

image2D((t(matrix(avg_image/100,469,541)))[index_plot_1,index_plot_2],
        x=testing_input_separable[[1]][[1]][index_plot_1],
        y=testing_input_separable[[1]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        zlim=c(min_value,max_value),main='average',mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = FALSE)


#dev.off()


#pdf("sgasp_measurement_bias_model_discrepancy_1_with_bar.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
i=1
index_not_na=which(is.na(testing_output[[i]])==F)
image_plot=t(matrix(predict_sgasp@measurement_bias_mean[[i]],541,469))

image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        y=testing_input_separable[[i]][[2]][index_plot_2],
        zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
 #colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5)

for(i in 2:3){
  index_not_na=which(is.na(testing_output[[i]])==F)
  image_plot=t(matrix(predict_sgasp@measurement_bias_mean[[i]],541,469))
  
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}
#dev.off()

#pdf("sgasp_measurement_bias_model_discrepancy_2.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
for(i in 4:5){
  index_not_na=which(is.na(testing_output[[i]])==F)
  image_plot=t(matrix(predict_sgasp@measurement_bias_mean[[i]],541,469))
  
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}

image_plot=t(matrix(predict_sgasp@delta_mean,541,469))
image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main='model discrepancy',mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
#dev.off()

#pdf("gasp_measurement_bias_model_discrepancy_1_with_bar.pdf",width=12,height=3.8)
par(mfrow=c(1,3))
i=1
index_not_na=which(is.na(testing_output[[i]])==F)
image_plot=t(matrix(predict_gasp@measurement_bias_mean[[i]],541,469))

image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2, colkey = F)
# colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5))

for(i in 2:3){
  index_not_na=which(is.na(testing_output[[i]])==F)
  image_plot=t(matrix(predict_gasp@measurement_bias_mean[[i]],541,469))
  
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],
          y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}
#dev.off()


#pdf("gasp_measurement_bias_model_discrepancy_2.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
for(i in 4:5){
  index_not_na=which(is.na(testing_output[[i]])==F)
  image_plot=t(matrix(predict_gasp@measurement_bias_mean[[i]],541,469))
  
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],
          y=testing_input_separable[[i]][[2]][index_plot_2],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          zlim=c(min_value,max_value),main=paste('measurement bias',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}

image_plot=t(matrix(predict_gasp@delta_mean,541,469))
image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        main='model discrepancy',mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
#dev.off()


#pdf("computer_model.pdf",width=14,height=7)

par(mfrow=c(1,2))

index_not_na=which(is.na(testing_output[[i]])==F)
image_plot=t(matrix(predict_gasp@math_model_mean_no_trend[[5]],541,469))

image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        main='GaSP',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,mgp=c(2.5,1,0),colkey=list(side = 4,width = 0.5,
                                                                              cex.axis = 1.5, cex.clab = 1.5,mgp=c(2.5,0.5,0)))

index_not_na=which(is.na(testing_output[[i]])==F)
image_plot=t(matrix(predict_sgasp@math_model_mean_no_trend[[5]],541,469))

image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],x=testing_input_separable[[i]][[1]][index_plot_1],
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        main='S-GaSP',cex.axis=1.5,cex.lab=1.5,cex.main=1.5,mgp=c(2.5,1,0),colkey=list(side = 4,width = 0.5,
                                                                                       cex.axis = 1.5, cex.clab = 1.5,mgp=c(2.5,0.5,0)))

#dev.off()

#pdf("gasp_computer_model.pdf",width=10,height=6)
par(mfrow=c(2,3))
i=1
image2D(matrix(predict_gasp@math_model_mean[[i]],541,469),zlim=c(-0.072,0.078),main='geophysical model with mean',cex.axis=1.5,cex.lab=1.5,colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5))

for(i in 2:5){
  image2D(matrix(predict_gasp@math_model_mean[[i]],541,469),zlim=c(-0.072,0.078),main='geophysical model with mean',cex.axis=1.5,cex.lab=1.5,colkey = F)
}
image2D(matrix(predict_gasp@math_model_mean_no_trend[[5]],541,469),zlim=c(-0.072,0.078),main='geophysical model without mean',cex.axis=1.5,cex.lab=1.5,colkey = F)
#dev.off()

#pdf("sgasp_computer_model.pdf",width=10,height=6)
par(mfrow=c(2,3))
i=1
image2D(matrix(predict_sgasp@math_model_mean[[i]],541,469),zlim=c(-0.072,0.078),main='geophysical model with mean',cex.axis=1.5,cex.lab=1.5,colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5))

for(i in 2:5){
  image2D(matrix(predict_sgasp@math_model_mean[[i]],541,469),zlim=c(-0.072,0.078),main='geophysical model with mean',cex.axis=1.5,cex.lab=1.5,colkey = F)
}
image2D(matrix(predict_sgasp@math_model_mean_no_trend[[5]],541,469),zlim=c(-0.072,0.078),main='geophysical model without mean',cex.axis=1.5,cex.lab=1.5,colkey = F)
#dev.off()


max(predict_sgasp@mean[[5]])

for(i in 1:5){
  plot(model_gasp@post_theta[,i],type='l',ylim=c(min(model_gasp@post_theta[,i],model_sgasp@post_theta[,i]),
                                                 max(model_gasp@post_theta[,i],model_sgasp@post_theta[,i])),
       col='red')
  lines(model_sgasp@post_theta[,i],col='blue')
}


for(i in 1:3){
  image2D(t(matrix(testing_output[[i]]/100,469,541)))
}
image2D(matrix(predict_sgasp@mean[[1]],541,469))

for(i in 1:3){
  image2D(matrix(predict_sgasp@math_model_mean,541,469))
}


image2D(matrix(predict_sgasp@delta_mean,541,469))
image2D(matrix(predict_gasp@delta_mean,541,469))


quilt.plot(x=(model_gasp@input[[1]][,1]), y=(model_gasp@input[[1]][,2]),
           z=colMeans(model_sgasp@post_measurement_bias[[3]]),nrow = 80, ncol = 80)


quilt.plot(x=(model_gasp@input[[1]][,1]), y=(model_gasp@input[[1]][,2]),
           z=colMeans(model_gasp@post_measurement_bias[[4]]),nrow = 80, ncol = 80)






image2D(t(t3_ascending_obs_matrix1_full))

record_mse=matrix(0,5,2)
record_mse_math_model=matrix(0,5,2)
for(i in 1:5){
  index_not_na=which(is.na(testing_output[[i]])==F)
  pred_matrix=t(matrix(predict_gasp@mean[[i]],541,469))
  
  record_mse[i,1]=mean((pred_matrix[index_not_na]-testing_output[[i]][index_not_na]/100)^2)
  
  pred_matrix=t(matrix(predict_sgasp@mean[[i]],541,469))
  
  record_mse[i,2]=mean((pred_matrix[index_not_na]-testing_output[[i]][index_not_na]/100)^2)
  
  pred_matrix=t(matrix(predict_gasp@math_model_mean[[i]],541,469))
  record_mse_math_model[i,1]=mean((pred_matrix[index_not_na]-testing_output[[i]][index_not_na]/100)^2)
  
  
  pred_matrix=t(matrix(predict_sgasp@math_model_mean[[i]],541,469))
  record_mse_math_model[i,2]=mean((pred_matrix[index_not_na]-testing_output[[i]][index_not_na]/100)^2)
  
}


record_post_theta=matrix(0,5,2)
record_post_theta[,1]=colMeans(model_gasp@post_theta)
record_post_theta[,2]=colMeans(model_sgasp@post_theta)


avg_pred_sgasp=0
for(i in 1:5){
  avg_pred_sgasp=avg_pred_sgasp+(predict_sgasp@mean[[i]])
}
avg_pred_sgasp=avg_pred_sgasp/5
avg_pred_sgasp=t(matrix(avg_pred_sgasp,541,469))

avg_pred_gasp=0
for(i in 1:5){
  avg_pred_gasp=avg_pred_gasp+(predict_gasp@mean[[i]])
}
avg_pred_gasp=avg_pred_gasp/5
avg_pred_gasp=t(matrix(avg_pred_gasp,541,469))

index_not_na_here=which(is.na(avg_image)==F)
mean((avg_pred_gasp[index_not_na_here]-avg_image[index_not_na_here]/100)^2)
mean((avg_pred_sgasp[index_not_na_here]-avg_image[index_not_na_here]/100)^2)


avg_pred_sgasp_cm=0
for(i in 1:5){
  avg_pred_sgasp_cm=avg_pred_sgasp_cm+(predict_sgasp@math_model_mean[[i]])
}
avg_pred_sgasp_cm=avg_pred_sgasp_cm/5
avg_pred_sgasp_cm=t(matrix(avg_pred_sgasp_cm,541,469))

avg_pred_gasp_cm=0
for(i in 1:5){
  avg_pred_gasp_cm=avg_pred_gasp_cm+(predict_gasp@math_model_mean[[i]])
}
avg_pred_gasp_cm=avg_pred_gasp_cm/5
avg_pred_gasp_cm=t(matrix(avg_pred_gasp_cm,541,469))

mean((avg_pred_gasp_cm[index_not_na_here]-avg_image[index_not_na_here]/100)^2)
mean((avg_pred_sgasp_cm[index_not_na_here]-avg_image[index_not_na_here]/100)^2)



#plot predictive mean 

#pdf("gasp_prediction_1.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
i=1
index_not_na=which(is.na(testing_output[[i]])==F)

image_plot=t(matrix(predict_gasp@mean[[i]],541,469))
image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),
        mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F) #colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5)

for(i in 2:3){
  index_not_na=which(is.na(testing_output[[i]])==F)

  image_plot=t(matrix(predict_gasp@mean[[i]],541,469))
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}
#dev.off()

#pdf("gasp_prediction_2.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
for(i in 4:5){
  index_not_na=which(is.na(testing_output[[i]])==F)

  image_plot=t(matrix(predict_gasp@mean[[i]],541,469))
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}

image_plot=avg_pred_gasp
image_plot[-index_not_na]=NA

image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main='prediction average',mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)

#dev.off()




#pdf("sgasp_prediction_1.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
i=1
index_not_na=which(is.na(testing_output[[i]])==F)

image_plot=t(matrix(predict_sgasp@mean[[i]],541,469))
image_plot[-index_not_na]=NA
image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),mgp=c(2.5,1,0),
        cex.axis=2,cex.lab=2,cex.main=2,colkey=F) #colkey = list(plot = T, side = 4,width = 0.5,cex.clab=1.5,cex.axis = 1.5)

for(i in 2:3){
  index_not_na=which(is.na(testing_output[[i]])==F)
  #image_plot=t(matrix(predict_sgasp@measurement_bias_mean[[i]],541,469))
  #image_plot= image_plot+ t(matrix(predict_sgasp@math_model_mean[[i]],541,469))
  
  image_plot=t(matrix(predict_sgasp@mean[[i]],541,469))
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}
#dev.off()

#pdf("sgasp_prediction_2.pdf",width=10,height=3.8)
par(mfrow=c(1,3))
for(i in 4:5){
  index_not_na=which(is.na(testing_output[[i]])==F)

  image_plot=t(matrix(predict_sgasp@mean[[i]],541,469))
  image_plot[-index_not_na]=NA
  image2D((t(image_plot))[index_plot_1,index_plot_2],
          x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
          y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main=paste('prediction',i),mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)
}


image_plot=avg_pred_sgasp
image_plot[-index_not_na]=NA

image2D((t(image_plot))[index_plot_1,index_plot_2],
        x=testing_input_separable[[i]][[1]][index_plot_1],xlab=expression(x[1]),ylab=expression(x[2]),clab='m/yr',
        y=testing_input_separable[[i]][[2]][index_plot_2],zlim=c(min_value,max_value),main='prediction average',mgp=c(2.5,1,0),cex.axis=2,cex.lab=2,cex.main=2,colkey = F)

#dev.off()




###trace plot 

theta_name=c(expression(theta[1]),expression(theta[2]),expression(theta[3]),
             expression(theta[4]),expression(theta[5]) )
mu_name=c(expression(mu[1]),expression(mu[2]),expression(mu[3]),
          expression(mu[4]),expression(mu[5]) )
beta_name=c(expression(beta[list(1,1)]),expression(beta[list(1,2)]),expression(beta[list(2,1)]),expression(beta[list(2,2)]),
            expression(beta[list(3,1)]),expression(beta[list(3,2)]),expression(beta[list(4,1)]),expression(beta[list(4,2)]),
            expression(beta[list(5,1)]),expression(beta[list(5,2)]))

eta_name=c(expression(eta[1]),expression(eta[2]),expression(eta[3]),
           expression(eta[4]),expression(eta[5]) )
sigma_2_name=c(expression(sigma[1]^2),expression(sigma[2]^2),expression(sigma[3]^2),
               expression(sigma[4]^2),expression(sigma[5]^2) )


beta_shared_name=c(expression(beta[1]),expression(beta[2]))


#pdf('trace_plot_theta_gasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot(model_gasp@post_theta[,i],type='l',xlab='s',ylab=theta_name[i] )
}
#dev.off()
##mean par
#pdf('trace_plot_mu_gasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot(model_gasp@post_individual_par[[i]][,5],type='l',xlab='s',ylab=mu_name[i])
}
#dev.off()

#pdf('trace_plot_beta_l_gasp_1.pdf',height=6,width=10)
par(mfrow=c(1,5))
for(i in 1:2){
  for(j in 1:2){
    plot(exp(model_gasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
  }
}
i=3
j=1
plot(exp(model_gasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
#dev.off()

#pdf('trace_plot_beta_l_gasp_2.pdf',height=6,width=10)
par(mfrow=c(1,5))
i=3
j=2
plot(exp(model_gasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
for(i in 4:5){
  for(j in 1:2){
    plot(exp(model_gasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
  }
}
#dev.off()


#pdf('trace_plot_eta_l_gasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot(exp(model_gasp@post_individual_par[[i]][,3]),type='l',xlab='s',ylab=eta_name[i])
}
#dev.off()

#pdf('trace_plot_sigma_2_l_gasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot((model_gasp@post_individual_par[[i]][,4]),type='l',xlab='s',ylab=sigma_2_name[i],mgp=c(2.5,1,0))
}
#dev.off()

#pdf('trace_rest_gasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
plot(exp(model_gasp@post_individual_par[[6]][,1]),type='l',xlab='s',ylab=beta_shared_name[1])
plot(exp(model_gasp@post_individual_par[[6]][,2]),type='l',xlab='s',ylab=beta_shared_name[2])
plot(exp(model_gasp@post_individual_par[[6]][,3]),type='l',xlab='s',ylab=expression(tau^2))

#dev.off()



###sgasp

#pdf('trace_plot_theta_sgasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot(model_sgasp@post_theta[,i],type='l',xlab='s',ylab=theta_name[i] )
}
#dev.off()
##mean par
#pdf('trace_plot_mu_sgasp.pdf',height=3,width=10)
par(mfrow=c(1,5))

for(i in 1:5){
  plot(model_sgasp@post_individual_par[[i]][,5],type='l',xlab='s',ylab=mu_name[i])
}
#dev.off()

#pdf('trace_plot_beta_l_sgasp_1.pdf',height=6,width=10)
par(mfrow=c(1,5))
for(i in 1:2){
  for(j in 1:2){
    plot(exp(model_sgasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
  }
}
i=3
j=1
plot(exp(model_sgasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
#dev.off()

#pdf('trace_plot_beta_l_sgasp_2.pdf',height=6,width=10)
par(mfrow=c(1,5))
i=3
j=2
plot(exp(model_sgasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
for(i in 4:5){
  for(j in 1:2){
    plot(exp(model_sgasp@post_individual_par[[i]][,j]),type='l',xlab='s',ylab=beta_name[(i-1)*2+j])
  }
}
#dev.off()


#pdf('trace_plot_eta_l_sgasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot(exp(model_sgasp@post_individual_par[[i]][,3]),type='l',xlab='s',ylab=eta_name[i])
}
#dev.off()

#pdf('trace_plot_sigma_2_l_sgasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
for(i in 1:5){
  plot((model_sgasp@post_individual_par[[i]][,4]),type='l',xlab='s',ylab=sigma_2_name[i],mgp=c(2.5,1,0))
}
#dev.off()

#pdf('trace_rest_sgasp.pdf',height=3,width=10)
par(mfrow=c(1,5))
plot(exp(model_sgasp@post_individual_par[[6]][,1]),type='l',xlab='s',ylab=beta_shared_name[1])
plot(exp(model_sgasp@post_individual_par[[6]][,2]),type='l',xlab='s',ylab=beta_shared_name[2])
plot(exp(model_sgasp@post_individual_par[[6]][,3]),type='l',xlab='s',ylab=expression(tau^2))

#dev.off()


