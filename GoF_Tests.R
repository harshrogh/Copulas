# Required Packages
require(rugarch)
require(zoo)
require(xts)
require(copula)
require(ggplot2)
require(gridExtra)
require("CopulaTests")

data <-read.csv("snp500-adj_close_2004-2018.csv")
# extract GSPC and VIX data as well to understand general market movements
GSPC <- zoo(data$GSPC,as.Date(as.character(data$Date), format = c("%Y-%m-%d")))
GSPC_ret <- log(GSPC/lag(GSPC,-1)) #get log returns
VIX <- zoo(data$VIX,as.Date(as.character(data$Date), format = c("%Y-%m-%d")))

# plot GSPC and VIX movement through the years
par(mfrow = c(1,3),mar = c(3,4,2,.2)+.5, oma = c(2,0,0,1))
plot(GSPC,main = "GSPC Price 2007-2017", xlab = "", ylab = "Price")
plot(GSPC_ret,main = "GSPC Log Price Returns 2007-2017", xlab = "", ylab = "Log Returns")
plot(VIX,main="VIX Volatility",xlab="",ylab="Volatility")
mtext(text = "Date", cex = 1.2, outer = TRUE, side = 1, line = .5, font = 2)

data <- data[,2:ncol(data)-2]
# Sectors
sec = c(2:4,14:16,23:25,40:42)  #[financials, tech, industrials, health]

# create indexed time series object
prices <- zoo(data[,sec],as.Date(as.character(data$Date), format = c("%Y-%m-%d")))
lreturns <- log(prices/lag(prices,-1)) #get log returns

### Find returns of assets in high volatility, low volatality, period for backtesting
ret_high <- as.xts(lreturns)['20070101/20101231']
ret <- as.xts(lreturns)['20120601/20170531']

# ARMA(1,1)-GARCH(1,1) with skewed t-distribution 
uspec.st <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
                      mean.model=list(armaOrder=c(1,0)), 
                      distribution.model = "sstd")

####### FIT THE GARCHSPEC on just one asset for testing
GARCH.WFC.st <- ugarchfit(spec = uspec.st, data = ret$WFC)  # pvalue for goodness-of-test > 0.05, use skewed because we assume financial data to be heavy tailed? 

#### Find residuals of each asset after fitting GARCH(1,1), ARMA(1,0) with skewed distribution on (2007-2017)
dim <- ncol(ret)
n <- nrow(ret)
GARCH.fit <- sapply(1:dim, function(i) {ugarchfit(spec = uspec.st, data = ret[,i])})  #fit GARCH
eps <- sapply(GARCH.fit, residuals, standardize=T) #find residual returns
colnames(eps) = colnames(ret)
eps_mat <- as.matrix(eps)
U_all <- pobs(eps_mat) # empirical cdf of overall asset prices

##### Save index of assets in each sector ######
# 3 EACH
Fin_ind<-c(1:3)
Tech_ind<-c(4:6)
Ind_ind<-c(7:9)
Health_ind<-c(10:12)

#### Fit GARCH(1,1) on high volatility period
dim.h <- ncol(ret_high)
n.h <- nrow(ret_high)
GARCH.fit.h <- sapply(1:dim.h, function(i) {ugarchfit(spec = uspec.st, data = ret_high[,i])})  #fit GARCH
eps.h <- sapply(GARCH.fit.h, residuals, standardize=T) #find residual returns
colnames(eps.h) = colnames(ret_high)
eps_mat.h <- as.matrix(eps.h)
U.h <- pobs(eps_mat.h) # empirical cdf of overall asset prices

# Start fitting models

#### NORMAL COPULA FIT ###
# Entire period
normFit <- fitCopula(normalCopula(dim = dim),data = U_all, method = "irho")
normCop <- normFit@copula # The fitted normal copula

# One copula per sector
# Finance
normFit_Fin <- fitCopula(normalCopula(dim = length(Fin_ind)),data = U_all[,Fin_ind], method = "irho")
normCop_Fin <- normFit_Fin@copula # The fitted normal copula
# Tech sector
normFit_Tech <- fitCopula(normalCopula(dim = length(Tech_ind)),data = U_all[,Tech_ind], method = "irho")
normCop_Tech <- normFit_Tech@copula # The fitted normal copula
# Industrials sector
normFit_Ind <- fitCopula(normalCopula(dim = length(Ind_ind)),data = U_all[,Ind_ind], method = "irho")
normCop_Ind <- normFit_Ind@copula # The fitted normal copula
# Health sector
normFit_Health <- fitCopula(normalCopula(dim = length(Health_ind)),data = U_all[,Health_ind], method = "irho")
normCop_Health <- normFit_Health@copula # The fitted normal copula

# High volatility
normFit.h <- fitCopula(normalCopula(dim = ncol(ret_high)),data = U.h, method = "irho")
normCop.h <- normFit.h@copula # The fitted normal copula

##### T COPULA FIT ###
# We use itau.mpl as it is the suggested estimator for t copulas in the copula package
tFit <- fitCopula(tCopula(dim = dim, dispstr = "un"), data = U_all, method = "itau.mpl")
tCop <- tFit@copula # The fitted t-copula

# One copula per sector
# Finance
tFit_Fin <- fitCopula(tCopula(dim = length(Fin_ind), dispstr = "un"), U_all[,Fin_ind], method = "itau.mpl")
tCop_Fin <- tFit_Fin@copula # The fitted t-copula
# Tech sector
tFit_Tech <- fitCopula(tCopula(dim = length(Tech_ind), dispstr = "un"), U_all[,Tech_ind], method = "itau.mpl")
tCop_Tech <- tFit_Tech@copula # The fitted t-copula
# Industrials sector
tFit_Ind <- fitCopula(tCopula(dim = length(Ind_ind), dispstr = "un"), U_all[,Ind_ind], method = "itau.mpl")
tCop_Ind <- tFit_Ind@copula # The fitted t-copula
# Health sector
tFit_Health <- fitCopula(tCopula(dim = length(Health_ind), dispstr = "un"), U_all[,Health_ind], method = "itau.mpl")
tCop_Health <- tFit_Health@copula # The fitted t-copula

# High volatility
tFit.h <- fitCopula(tCopula(dim = ncol(ret_high), dispstr = "un"), U.h, method = "itau.mpl")
tCop.h <- tFit.h@copula # The fitted t-copula

###### GUMBEL COPULA FIT ####
gumFit <- fitCopula(gumbelCopula(dim = dim),data = U_all, method = "irho")
gumCop <- gumFit@copula # The fitted gumbel copula

# One copula per sector
# Finance
gumFit_Fin <- fitCopula(gumbelCopula(dim = length(Fin_ind)),data = U_all[,Fin_ind], method = "irho")
gumCop_Fin <- gumFit_Fin@copula # The fitted gumbel copula
# Tech sector
gumFit_Tech <- fitCopula(gumbelCopula(dim = length(Tech_ind)),data = U_all[,Tech_ind], method = "irho")
gumCop_Tech <- gumFit_Tech@copula # The fitted gumbel copula
# Industrials sector
gumFit_Ind <- fitCopula(gumbelCopula(dim = length(Ind_ind)),data = U_all[,Ind_ind], method = "irho")
gumCop_Ind <- gumFit_Ind@copula # The fitted gumbel copula
# Health sector
gumFit_Health <- fitCopula(gumbelCopula(dim = length(Health_ind)),data = U_all[,Health_ind], method = "irho")
gumCop_Health <- gumFit_Health@copula # The fitted gumbel copula

# High volatility
gumFit.h <- fitCopula(gumbelCopula(dim = ncol(ret_high)),data = U.h, method = "irho")
gumCop.h <- gumFit.h@copula # The fitted gumbel copula

###### Clayton COPULA FIT ####
clayFit <- fitCopula(claytonCopula(dim = dim),data = U_all, method = "itau")
clayCop <- clayFit@copula # The fitted clayton copula

# One copula per sector
# Finance
clayFit_Fin <- fitCopula(claytonCopula(dim = length(Fin_ind)),data = U_all[,Fin_ind], method = "itau")
clayCop_Fin <- clayFit_Fin@copula # The fitted clayton copula
# Tech sector
clayFit_Tech <- fitCopula(claytonCopula(dim = length(Tech_ind)),data = U_all[,Tech_ind], method = "itau")
clayCop_Tech <- clayFit_Tech@copula # The fitted clayton copula
# Industrials sector
clayFit_Ind <- fitCopula(claytonCopula(dim = length(Ind_ind)),data = U_all[,Ind_ind], method = "itau")
clayCop_Ind <- clayFit_Ind@copula # The fitted clayton copula
# Health sector
clayFit_Health <- fitCopula(claytonCopula(dim = length(Health_ind)),data = U_all[,Health_ind], method = "itau")
clayCop_Health <- clayFit_Health@copula # The fitted clayton copula

# High volatility
clayFit.h <- fitCopula(claytonCopula(dim = ncol(ret_high)),data = U.h, method = "itau")
clayCop.h <- clayFit.h@copula # The fitted clayton copula

### Goodness of Fit
norm = GoFcopula(dim=dim,n=n,U=U_all,copula="normal",copulaClass=gumCop,ntest=100)
norm_Finance = GoFcopula(dim=length(Fin_ind),n=n,U=U_all[,Fin_ind],copula="normal",copulaClass=normCop_Fin,ntest=100)
norm_Tech = GoFcopula(dim=length(Tech_ind),n=n,U=U_all[,Tech_ind],copula="normal",copulaClass=normCop_Tech,ntest=100)
norm_Industry = GoFcopula(dim=length(Ind_ind),n=n,U=U_all[,Ind_ind],copula="normal",copulaClass=normCop_Ind,ntest=100)
norm_Health = GoFcopula(dim=length(Health_ind),n=n,U=U_all[,Health_ind],copula="normal",copulaClass=normCop_Health,ntest=100)
norm_HV = GoFcopula(dim=dim,n=nrow(U.h),U=U.h,copula="normal",copulaClass=normCop.h,ntest=100)

tdist = GoFcopula(dim=dim,n=n,U=U_all,copula="std",copulaClass=gumCop,ntest=100)
tdist_Finance = GoFcopula(dim=length(Fin_ind),n=n,U=U_all[,Fin_ind],copula="std",copulaClass=tCop_Fin,ntest=100)
tdist_Tech = GoFcopula(dim=length(Tech_ind),n=n,U=U_all[,Tech_ind],copula="std",copulaClass=tCop_Tech,ntest=100)
tdist_Industry = GoFcopula(dim=length(Ind_ind),n=n,U=U_all[,Ind_ind],copula="std",copulaClass=tCop_Ind,ntest=100)
tdist_Health = GoFcopula(dim=length(Health_ind),n=n,U=U_all[,Health_ind],copula="std",copulaClass=tCop_Health,ntest=100)
tdist_HV = GoFcopula(dim=dim,n=nrow(U.h),U=U.h,copula="std",copulaClass=tCop.h,ntest=100)

gumbel = GoFcopula(dim=dim,n=n,U=U_all,copula="gumbel",copulaClass=gumCop,ntest=100)
gumbel_Finance = GoFcopula(dim=length(Fin_ind),n=n,U=U_all[,Fin_ind],copula="gumbel",copulaClass=gumCop_Fin,ntest=100)
gumbel_Tech = GoFcopula(dim=length(Tech_ind),n=n,U=U_all[,Tech_ind],copula="gumbel",copulaClass=gumCop_Tech,ntest=100)
gumbel_Industry = GoFcopula(dim=length(Ind_ind),n=n,U=U_all[,Ind_ind],copula="gumbel",copulaClass=gumCop_Ind,ntest=100)
gumbel_Health = GoFcopula(dim=length(Health_ind),n=n,U=U_all[,Health_ind],copula="gumbel",copulaClass=gumCop_Health,ntest=100)
gumbel_HV = GoFcopula(dim=dim,n=nrow(U.h),U=U.h,copula="gumbel",copulaClass=gumCop.h,ntest=100)

clay = GoFcopula(dim=dim,n=n,U=U_all,copula="clayton",copulaClass=clayCop,ntest=100)
clay_Finance = GoFcopula(dim=length(Fin_ind),n=n,U=U_all[,Fin_ind],copula="clayton",copulaClass=clayCop_Fin,ntest=100)
clay_Tech = GoFcopula(dim=length(Tech_ind),n=n,U=U_all[,Tech_ind],copula="clayton",copulaClass=clayCop_Tech,ntest=100)
clay_Industry = GoFcopula(dim=length(Ind_ind),n=n,U=U_all[,Ind_ind],copula="clayton",copulaClass=clayCop_Ind,ntest=100)
clay_Health = GoFcopula(dim=length(Health_ind),n=n,U=U_all[,Health_ind],copula="clayton",copulaClass=clayCop_Health,ntest=100)
clay_HV = GoFcopula(dim=dim,n=nrow(U.h),U=U.h,copula="clayton",copulaClass=clayCop.h,ntest=100)


# Check VaR violations, if violations > 5%, proposed model is not a good fit
# We look at the returns of a portfolio of all stocks with equal weights

# Testing data
backtest <- as.xts(lreturns)['20170601/20180330'] 
act_returns <- rowMeans(backtest)
dates <- index(backtest)  #keep a list of dates, will be used later to plot
bt_highVol <- as.xts(lreturns)['20110101/20111231'] # testing data for high volatility models
highVol_returns <- rowMeans(bt_highVol)
highVol_dates <- index(bt_highVol)


# values for testing
alpha <- 0.05 #if violations deviated very far from alpha, proposed model is not a good fit
days <- nrow(backtest) # days we are backtesting against
days_highVol <- nrow(bt_highVol)
dates <- index(backtest)  #keep a list of dates, will be used later to plot
sims <- 500   # number of simulations that will be done to get the VaR

### normal copula + GARCH(1,1)
# Test overall period (2012-2017)
normVaR <- simulateVaR(GARCH.fit,normCop,alpha,sims,days)
normViol <- violatedVaRs(normVaR,act_returns)
# Sector based copula model
normVaR_sec <- (simulateVaR(GARCH.fit[Fin_ind],normCop_Fin,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Tech_ind],normCop_Tech,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Ind_ind],normCop_Ind,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Health_ind],normCop_Health,alpha,sims,days))/4 # assuming equal weights in portfolio
normViol_sec <- violatedVaRs(normVaR_sec,act_returns)
# high volatility period
normVaR_high <- simulateVaR(GARCH.fit.h,normCop.h,alpha,sims,days_highVol)
normViol_high <- violatedVaRs(normVaR_high,highVol_returns)

### tdist copula + GARCH(1,1)
# Test overall period (2012-2017)
tdistVaR <- simulateVaR(GARCH.fit,tCop,alpha,sims,days)
tdistViol <- violatedVaRs(tdistVaR,act_returns)
# Sector based copula model
tdistVaR_sec <- (simulateVaR(GARCH.fit[Fin_ind],tCop_Fin,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Tech_ind],tCop_Tech,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Ind_ind],tCop_Ind,alpha,sims,days)+
                 simulateVaR(GARCH.fit[Health_ind],tCop_Health,alpha,sims,days))/4
tdistViol_sec <- violatedVaRs(tdistVaR_sec,act_returns)
# high volatility period
tVaR_high <- simulateVaR(GARCH.fit.h,tCop.h,alpha,sims,days_highVol)
tViol_high <- violatedVaRs(tVaR_high,highVol_returns)

### gumbel copula + GARCH(1,1)
# Test overall period (2012-2017)
gumVaR <- simulateVaR(GARCH.fit,gumCop,alpha,sims,days)
gumViol <- violatedVaRs(gumVaR,act_returns)
# Sector based copula model
gumVaR_sec <- (simulateVaR(GARCH.fit[Fin_ind],gumCop_Fin,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Tech_ind],gumCop_Tech,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Ind_ind],gumCop_Ind,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Health_ind],gumCop_Health,alpha,sims,days))/4
gumViol_sec <- violatedVaRs(gumVaR_sec,act_returns)
# high volatility period
gumVaR_high <- simulateVaR(GARCH.fit.h,gumCop.h,alpha,sims,days_highVol)
gumViol_high <- violatedVaRs(gumVaR_high,highVol_returns)

### clayton copula + GARCH(1,1)
# Test overall period (2012-2017)
clayVaR <- simulateVaR(GARCH.fit,clayCop,alpha,sims,days)
clayViol <- violatedVaRs(clayVaR,act_returns)
# Sector based copula model
clayVaR_sec <- (simulateVaR(GARCH.fit[Fin_ind],clayCop_Fin,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Tech_ind],clayCop_Tech,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Ind_ind],clayCop_Ind,alpha,sims,days)+
                  simulateVaR(GARCH.fit[Health_ind],clayCop_Health,alpha,sims,days))/4
clayViol_sec <- violatedVaRs(clayVaR_sec,act_returns)
# high volatility period
clayVaR_high <- simulateVaR(GARCH.fit.h,clayCop.h,alpha,sims,days_highVol)
clayViol_high <- violatedVaRs(clayVaR_high,highVol_returns)

# plot actual returns against simulated VaR

# plot
theme_update(plot.title = element_text(hjust = 0.5),
             panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.x = element_line(color="black",size=0.2),
             axis.line.y = element_line(color="black", size = 0.2))

p1 <- ggplot(as.data.frame(normVaR), aes(dates)) +
  geom_line(aes(y=normVaR), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") +
  labs(x="",y="Normal",title="Normal Volatility (2012-2017)")

p2 <- ggplot(as.data.frame(normVaR_sec), aes(dates)) +
  geom_line(aes(y=normVaR_sec), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") + 
  labs(x="",y="",title="Sector Level")

p3 <- ggplot(as.data.frame(normVaR_high), aes(highVol_dates)) +
  geom_line(aes(y=normVaR_high), colour="red") + 
  geom_point(aes(y=highVol_returns), colour="gray40") + 
  labs(x="",y="",title="High Volatility Period")

p4 <- ggplot(as.data.frame(tdistVaR), aes(dates)) +
  geom_line(aes(y=tdistVaR), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") + 
  labs(x="",y="t")

p5 <- ggplot(as.data.frame(tdistVaR_sec), aes(dates)) +
  geom_line(aes(y=tdistVaR_sec), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") +
  labs(x="",y="")

p6 <- ggplot(as.data.frame(tVaR_high), aes(highVol_dates)) +
  geom_line(aes(y=tVaR_high), colour="red") + 
  geom_point(aes(y=highVol_returns), colour="gray40") +
  labs(x="",y="")

p7 <- ggplot(as.data.frame(gumVaR), aes(dates)) +
  geom_line(aes(y=gumVaR), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") + 
  labs(x="",y="Gumbel")

p8 <- ggplot(as.data.frame(gumVaR_sec), aes(dates)) +
  geom_line(aes(y=gumVaR_sec), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") + 
  labs(x="",y="")

p9 <- ggplot(as.data.frame(gumVaR_high), aes(highVol_dates)) +
  geom_line(aes(y=gumVaR_high), colour="red") + 
  geom_point(aes(y=highVol_returns), colour="gray40") +
  labs(x="",y="")

p10 <- ggplot(as.data.frame(clayVaR), aes(dates)) +
  geom_line(aes(y=clayVaR), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") +
  labs(x="",y="Clayton")

p11 <- ggplot(as.data.frame(clayVaR_sec), aes(dates)) +
  geom_line(aes(y=clayVaR_sec), colour="red") + 
  geom_point(aes(y=act_returns), colour="gray40") +
  labs(x="",y="")

p12 <- ggplot(as.data.frame(clayVaR_high), aes(highVol_dates)) +
  geom_line(aes(y=clayVaR_high), colour="red") + 
  geom_point(aes(y=highVol_returns), colour="gray40") + 
  labs(x="",y="")

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,nrow=4,ncol=3)
