#####################
# Bikeshare Dataset
#####################
#####################
# Load the data
#####################
bike.all=read.csv('hour.csv')  

# choose non-working day
ind = (bike.all$workingday==0) 
bike.tmp = bike.all[ind,]

# remove the dates which does not have full hours
tmp=split(bike.tmp[,c("hr")],bike.tmp$dteday)
rm.list = names(which(unlist(lapply(tmp,length))!=24))
bike = bike.tmp[!(bike.tmp$dteday %in% rm.list),]


#temp: Normalized temperature in Celsius. The values are divided to 41 (max)
bike$dteday = as.Date(bike$dteday,format='%Y-%m-%d')
dates = (unique(bike$dteday))
dates = sort(as.Date(dates, format="%m/%d/%Y"))


hr = unique(bike$hr)


n = length(dates)
nt = length(hr)
tt=seq(0,1,len=nt)

dat=array(0,c(nt,n,4))
dimnames(dat) = list(hr/24,dates, c('temp','count', 'humid','windspeed'))

for(i in 1:n){
  index = which(bike$dteday==dates[i])
  dat[,i,1] = bike$temp[index]*41
  dat[,i,2] = bike$cnt[index]
  dat[,i,3] = bike$hum[index]*100
  dat[,i,4] = bike$windspeed[index]*67
  bike$temp_high[index] = max(bike$temp[index]*41)
  bike$temp_low[index] = min(bike$temp[index]*41)
  
  bike$temp_avg[index] = mean(bike$temp[index]*41)
  bike$logcnt[index] = log(bike$cnt[index]+1)
}

#####################
# Plot
#####################
library(ggplot2)
ttt = 1:24

dat.df=bike
dat.df$temp=bike$temp*41

index.date = (dat.df$dteday %in% dates[1:100])   # first 50 days

p <- ggplot(dat.df[index.date,], aes(x=hr, y=logcnt,group=dteday, color=temp_avg)) + 
  stat_smooth(method = lm, formula = y ~ splines::bs(x, 11), se=F)+
  labs(title = paste0('Daily bike rental pattern (w/ avg-temp)'), x='Hour', y='log(Count+1)')+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=20))
# 9 X 7 
p
p <- ggplot(dat.df[index.date,], aes(x=hr, y=cnt,group=dteday, color=temp_low)) + 
  stat_smooth(method = lm, formula = y ~ splines::bs(x, 11), se=F)+
  labs(title = paste0('Bike rental pattern (w/ low-temp)'), x='Hour', y='Count (bike rental)')+
  theme_bw() 
  
p


p <- ggplot(dat.df[index.date,], aes(x=hr, y=cnt,group=dteday, color=-temp_high)) + 
  stat_smooth(method = lm, formula = y ~ splines::bs(x, 11), se=F)+
  labs(title = paste0('Bike rental pattern (w/ low-temp)'), x='Hour', y='Count (bike rental)')+
  theme_bw() 

p

p <- ggplot(dat.df[index.date,], aes(x=hr, y=cnt,group=dteday, color=dteday)) + 
  stat_smooth(method = lm, formula = y ~ splines::bs(x, 11), se=F)   +
  labs(title = paste0('Bike rental pattern (w/ dates)'), x='Hour', y='Count (bike rental)')+
  theme_bw() 
p



p <- ggplot(dat.df[index.date,], aes(x=hr, y=temp,group=dteday, colour=temp_avg)) + 
  stat_smooth(method = lm, formula = y ~ splines::bs(x, 11), se=F)   +
  labs(title = paste0('Daily temperature (w/ avg-temp)'), x='Hour', y='Temperature in Celsius')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20))

p



pp <- p + geom_line() + theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Model I with Y'), y='X(t)')

#####################
# Load the results
#####################
load('bike_nonlinear_gauss_01242019.RData')
load('bikeshare_11_12.RData')

colnames(mse.1112)=1:10
table.df=as.data.frame(as.table(mse.1112))
names(table.df) = c("method", "dimension", "mse")
p = ggplot(table.df, aes(x=dimension,y=mse, group=method, colour=method)) + geom_point()
table.df



rownames(mse.1112)= c("linear", "nonlinear", "FPCA/Linear", "FPCA/Nonlinear", "WIRE/linar","WAVE/linear","WDR/linear", "WIRE","WAVE","WDR")
rownames(mse.1112.train)= c("linear", "nonlinear", "FPCA/Linear", "FPCA/Nonlinear", "WIRE/linar","WAVE/linear","WDR/linear", "WIRE","WAVE","WDR")
mse.1112.arr =array(0,c(2,5,10))  # linear/nonlinear X reg/fpca/wire/wave/wdr X 1:10 dimension
mse.1112.train.arr =array(0,c(2,5,10))  # linear/nonlinear X reg/fpca/wire/wave/wdr X 1:10 dimension
mse.1112.arr[1,,] = mse.1112[c(1,3,5,6,7),]
mse.1112.arr[2,,] = mse.1112[-c(1,3,5,6,7),]
mse.1112.train.arr[1,,] = mse.1112.train[c(1,3,5,6,7),]
mse.1112.train.arr[2,,] = mse.1112.train[-c(1,3,5,6,7),]


source('~/Dropbox (UNC Charlotte)/Codes/for_Paper/maketable.R')
make.table(mse.all,way=2, showsd=TRUE, rdnum=2, showname=TRUE)
make.table(mse.1112[,1:7],way=2, showsd=TRUE, rdnum=2, showname=TRUE)
make.table(mse.1112.train[,1:7],way=2, showsd=TRUE, rdnum=2, showname=TRUE)

