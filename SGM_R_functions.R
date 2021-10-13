#Methods for the stochastic growth model
#Robin Ristl
#Center for Medical Statistics, Informatics and Intelligent Systems
#Medical University of Vienna
#robin.ristl@meduniwien.ac.at
#10 October 2021

#Function to fit the stochastic growth model to data
#y ... vector of maximum diameter values
#time ... vector of measurement time points 
#id ... vector of subject id codes
#jackknife ... logical, use jackknife for parameter variance estimation, default is FALSE
#laplace ... logical, if TRUE use laplace appproximation, otherwise directly use likelihood, default is TRUE
#M ... in numeric integration over normal distribution, integrate over mean +/- M times the standarddeviation default is 8
#showplot ... logical, show plots regarding model fit, default is TRUE
sgm_fun<-function(y,time,id,jackknife=FALSE,laplace=TRUE,M=8,showplot=TRUE) {
	dat0<-prepare_dat_fun(y,time,id)
	dat<-na.omit(dat0)[,c("id","delta_log_y","delta_time")]
	r<-dat$delta_log_y
	t<-dat$delta_time
	id<-dat$id

	ids<-unique(id)
	#Starting values:
	indi_lambda<-indi_sd<-rep(NA,length(ids))
	for(i in 1:length(indi_lambda)) {
		mod_i<-lm(r~t-1,weights=1/t,subset=id==ids[i])
		indi_lambda[i]<-coef(mod_i)
		indi_sd[i]<-summary(mod_i)$sigma
	}
	
	set_indi_lambda<-indi_lambda>0
	indi_sd<-indi_sd[set_indi_lambda]
	indi_lambda<-indi_lambda[set_indi_lambda]

	if(showplot) {
		par(mfrow=c(3,2))
		boxplot(indi_lambda)
		boxplot(log(indi_lambda))
		qqnorm(indi_lambda)
		qqnorm(log(indi_lambda))
	}
	theta<-mean(log(indi_lambda))
	nu<-sd(log(indi_lambda))

	modo<-lm(indi_sd~indi_lambda-1)
	if(showplot) {
		plot(indi_lambda,indi_sd)
		abline(0,coef(modo))
	}
	k<-as.numeric(coef(modo))
	k*(exp(theta))
	
	modo2<-lm(log(indi_sd)~log(indi_lambda)-1)
	if(showplot) {
		plot(log(indi_lambda),log(indi_sd))
		abline(0,coef(modo2))
	}
	k<-as.numeric(coef(modo2))
	exp(k*theta)
		
	likelihood_i<-function(theta,log_nu,log_k,pat,tol=10^-6,maxiter=40) {
		nu<-exp(log_nu)
		nu2<-nu^2

		k<-exp(log_k)
	
		set<-id==pat
		y<-r[set]
		z<-t[set]
		N<-length(y)		

		if(!laplace) {
		intfun<-function(beta) {
			lambda<-exp(beta)
			sigma<-k*lambda # *sqrt(zeit)
			sigma2<-sigma^2
			mu<-(lambda-sigma2/2)*z
			#c(prod(dnorm(y,mu,sigma))*dnorm(beta,theta,nu)
			exp(sum(dnorm(y,mu,sigma*sqrt(z),log=TRUE))+dnorm(beta,theta,nu,log=TRUE))
			#exp(sum(dnorm(y,mu,sigma,log=TRUE))+dnorm(beta,theta,nu,log=TRUE))/M
		}
		intfun_vek<-function(vek) sapply(vek,intfun)
		
		LOW<- theta-M*nu
		UP<- theta+M*nu

		Int<-integrate(intfun_vek,lower=LOW,upper=UP)
		return(log(Int$value))
		
		} else {


		#Laplace Approximation
				
		f0<-function(x) sum( -(y-(exp(x)-exp(2*x)*k^2/2)*z)^2 /(2*k^2*exp(2*x)*z) )  -(x-theta)^2/(2*nu2)-x*N
     		f1<-function(x) sum( -(exp(-2*x)*(k^2*z*exp(2*x)-2*y)*(k^2*z*exp(2*x)-2*z*exp(x)+2*y))/(4*k^2*z) )  -(x-theta)/nu2-N
		f2<-function(x) sum( -(exp(-2*x)*(k^4*z^2*exp(4*x)-k^2*z^2*exp(3*x)-2*y*z*exp(x)+4*y^2))/(2*k^2*z) ) -1/nu2

		LOW<- theta-M*nu
		UP<- theta+M*nu

		#x0 as argmax of f0 over beta
		x0<-optimize(f0,interval=c(LOW,UP),maximum=TRUE)$maximum
		
		logIntLaplace<- -N/2 *log(2*pi)  -0.5*sum(log(z)) -0.5*log(2*pi)+ 0.5*log(2*pi) -N*log(k) -log(nu) -0.5*log(-f2(x0))+f0(x0)
	
		return(logIntLaplace)
		}
	}

	likelihood_all<-function(param,...) {
		loglik<-rep(NA,length(ids))
		for(j in 1:length(ids)) loglik[j]<-likelihood_i(param[1],param[2],param[3],ids[j],...)
		-sum(loglik)
	}

	likelihood_all_but_1<-function(param,ii,...) {
		ids_ii<-ids[-ii]
		loglik<-rep(NA,length(ids_ii))
		for(j in 1:length(ids_ii)) loglik[j]<-likelihood_i(param[1],param[2],param[3],ids_ii[j],...)
		-sum(loglik)
	}
	
	para0<-c(theta,log(nu),log(k))
	optL<-optim(para0,likelihood_all,hessian=TRUE)
	loglik<- -optL$value
	
	est<-c(optL$par[1],exp(optL$par[2:3]))
	
	V<-diag(c(1,est[2],est[3]))%*%solve(optL$hessian)%*%diag(c(1,est[2],est[3]))
	SE<-sqrt(diag(V))

	if(jackknife) {
		nn<-length(ids)
		est_jn<-matrix(NA,nrow=nn,ncol=3)
		for(ii1 in 1:length(ids)) {
			print(paste("Jackknife",ii1,"of",nn))
			flush.console()
			#optL$par
			#para0
			optLii<-optim(optL$par,likelihood_all_but_1,ii=ii1)#,method="BFGS")#,lower=c(-Inf,0,0),upper=c(Inf,Inf,Inf),method="L-BFGS-B",control=list(trace=3,pgtol=0))
			est_jn[ii1,]<-optLii$par
		}
		colMeans(est_jn)-optL$par
		est_jn[,2:3]<-exp(est_jn[,2:3])
		zent<-est_jn-rep(colMeans(est_jn),each=dim(est_jn)[1])
		Vjn<-(nn-1)/nn*t(zent)%*%zent
		SEjn<-sqrt(diag(Vjn))
	} else {
		Vjn<-NA
		SEjn<-NA
	}

	out<-list(
		estimate=data.frame(name=c("theta","nu","k"),est=est,SE=SE,SEjn=SEjn),
		loglik=loglik,
		V=V,
		Vjn=Vjn,
		data=dat0,
		indi_sd=indi_sd,
		indi_lambda=indi_lambda,
		n=length(ids),
		id=id
	)
	class(out)<-"sgm"
	return(out)
}


#Helper function for sgm_fun
prepare_dat_fun<-function(y,time,id) {
	dat<-data.frame(y,time,id)
	dat<-na.omit(dat)
	dat<-dat[order(dat$id,dat$time),]
	ids<-unique(dat$id)

	dat$log_y<-log(dat$y)

	dat$delta_log_y<-NA
	dat$delta_time<-NA
	dat$ratio_from_baseline<-NA

	for(i in ids) {
		set<-dat$id==i
		dat$delta_log_y[set][-1]<-diff(dat$log_y[set])
		dat$delta_time[set][-1]<-diff(dat$time[set])
		dat$ratio_from_baseline[set]<-dat$y[set]/(dat$y[set][1])

	}
	dat
}


#Generic functions for models of class sgm as produced by sgm_fun()
#Print a summary of a fitted model object of class sgm
#x ... a fitted model object of class sgm
print.sgm<-function(x,...) {
	print(x$estimate)
}

#Plot relative growth from raw data and from model fit
#x ... a fitted model object of class sgm
#p ... a vector of probabilities defining which quantiles to include in the plot
#steps ... number of points over the time span to plot a quantile line
#LTYqu ... numeric vector of line types for the quantile lines
#COL ... colour of plotted data
#COLqu ... colour of plotted quantile lines
#LWDqu ... line width for quantile lines
#XLIM ... xlim argument for the plot window, if NULL the default plot window is used 
#YLIM ... ylim argument for the plot window, if NULL the default plot window is used
#percent ... logical, defines if y-axis is in relative growth should be multiplied by 100 to be plotted in percent, default is FALSE
#XLAB ... label for x-axis
#YLAB ... label for y-asix
#PLOTTYPE ... type argument for the generic plot function for plotting the raw data trajectories
#PLOTPCH ... plotting character for raw data points
#add ... logical, defines whether the plot should be added to an existing plot window, defaults to FALSE
plot.sgm<-function(x,p=c(0.975,0.9,0.75,0.5,0.25,0.1,0.025),steps=20,LTYqu=c(3,4,2,1,2,4,3),
COL=1,COLqu=2,LWDqu=2,XLIM=NULL,YLIM=NULL,percent=FALSE,XLAB=NULL,YLAB=NULL,PLOTTYPE="o",PLOTPCH=1,add=FALSE,...) {	
	if(percent) faktor<-100 else faktor<-1
	if(!add) {
		plot.new()
		if(is.null(XLIM)) XLIM<-c(0,max(x$data$time))
		if(is.null(YLIM)) YLIM<-range(x$data$ratio_from_baseline)*faktor

		plot.window(xlim=XLIM,ylim=YLIM)
		axis(1)
		axis(2)
		if(is.null(YLAB)) YLAB="Relative size"
		if(is.null(XLAB)) XLAB="Years from initial observation"

		title(xlab=XLAB,ylab=YLAB)
		box()
	}
	for(ii in unique(x$data$id)) {
		set<-x$data$id==ii
		lines(x$data$time[set],faktor*x$data$ratio_from_baseline[set],type=PLOTTYPE,pch=PLOTPCH,col=COL)
	}
	if(length(p)>0) {
		Qu<-quantile(x,p=p,steps=steps)
		matlines(Qu[,1],faktor*exp(Qu[,-1]),lty=LTYqu,col=COLqu,lwd=LWDqu)
	}
}


#Function to calculate the vaule of a quantile of the estimated growth distribution
#zeit ... time interval for which to predict the growth distribution
#p ... which quantile (values between 0 and 1), default is 0.5 corresponding to the median
#theta_est ... model parameter theta
#nu_est ... model parameter nu_est
#k_est ... model paramter k
#modell ... an fitted model object of class sgm, if this is provided, theta_est, nu_est and k_est are extracted from the model object
qfun<-function(zeit,p=0.5,theta_est=NULL,nu_est=NULL,k_est=NULL,modell=NULL) {
	if(!is.null(modell)) {
		theta_est<-modell$estimate$est[modell$estimate$name=="theta"]
		nu_est<-modell$estimate$est[modell$estimate$name=="nu"]
		k_est<-modell$estimate$est[modell$estimate$name=="k"]

	}
	if(zeit==0) {
		ret<-0
	} else {
		rootfun<-function(q) {
			faltfun(q,zeit,theta_est,nu_est,k_est)-p
		}
		ret<-uniroot(rootfun,interval=c(log(0.05),log(20)))$root
	}
	ret
}

#Function to calculate the vaule of a quantile of the estimated growth distribution conditional on two ore more previous measurements
#zeit ... time interval for which to predict the growth distribution
#p ... which quantile (values between 0 and 1), default is 0.5 corresponding to the median
#y ... vector of previous maximum diameter values
#t ... vector of measurment time points for the previous measurementes
#theta_est ... model parameter theta
#nu_est ... model parameter nu_est
#k_est ... model paramter k
#modell ... an fitted model object of class sgm, if this is provided, theta_est, nu_est and k_est are extracted from the model object
qfun_update<-function(zeit,p=0.5,y,t,theta_est=NULL,nu_est=NULL,k_est=NULL,modell=NULL,...) {
	if(!is.null(modell)) {
		theta_est<-modell$estimate$est[modell$estimate$name=="theta"]
		nu_est<-modell$estimate$est[modell$estimate$name=="nu"]
		k_est<-modell$estimate$est[modell$estimate$name=="k"]

	}
	if(zeit==0) {
		ret<-0
	} else {
		rootfun<-function(q) {
			faltfun_update(q,zeit,y,t,theta_est,nu_est,k_est,...)-p
		}
		ret<-uniroot(rootfun,interval=c(log(0.05),log(20)))$root
	}
	ret
}


#Function to caculate the probability that the maximum diameter stays below a given threshold value within a given time interval
#mod ... fitted model object of class sgm
#y0 ... current maximum diameter value
#y_crit ... the threshold value
#years ... time interval for growth
guide_fun<-function(mod,y0,y_crit,years) {
	log_relative_change<-log(y_crit/y0)
	prob<-faltfun(log_relative_change,years,mod$estimate$est[1],mod$estimate$est[2],mod$estimate$est[3])
	prob
}

#Function to caculate the probability that the maximum diameter stays below a given threshold value within a given time interval conditional on two ore more previous measurements
#mod ... fitted model object of class sgm
#y0 ... current maximum diameter value (y0 will typically be the same value as the last value in y, but other values may be used)
#y_crit ... the threshold value
#years ... time interval for growth
#y ... vector of previous maximum diameter values
#t ... vector of measurment time points for the previous measurementes
guide_fun_update<-function(mod,y0,y_crit,years,y,t) {
	log_relative_change<-log(y_crit/y0)
	prob<-faltfun_update(q=log_relative_change,zeit=years,y=y,t=t,mod$estimate$est[1],mod$estimate$est[2],mod$estimate$est[3])
	prob
}

#Helper functions
faltfun<-function(q,zeit,theta_est,nu_est,k_est) {
	intfun<-function(beta,q) {
		lambda<-exp(beta)
		sigma<-k_est*lambda #hier war in 1-3 noch *sqrt(zeit), gehört aber nicht hierher, steht jetzt in pnorm.
		sigma2<-sigma^2
		pnorm(q,(lambda-sigma2/2)*zeit,sigma*sqrt(zeit))*dnorm(beta,theta_est,nu_est)
	}
	LOW<- theta_est-6*nu_est
	UP<- theta_est+6*nu_est
	Int<-integrate(intfun,lower=LOW,upper=UP,q=q)
	Int$value
}




faltfun_update<-function(q,zeit,y,t,theta_est,nu_est,k_est,nn=300,M=6) {
	z<-diff(t)
	r<-diff(log(y))
	beta<-theta_est

	k<-k_est
	theta<-theta_est
	nu<-nu_est
	nu2<-nu^2
	N<-length(r)
	#BLUP
	f0<-function(x) sum( -(r-(exp(x)-exp(2*x)*k^2/2)*z)^2 /(2*k^2*exp(2*x)*z) )  -(x-theta)^2/(2*nu2)-x*N

	LOW<- theta-M*nu
	UP<- theta+M*nu

	x0<-optimize(f0,interval=c(LOW,UP),maximum=TRUE)$maximum

	bed_dichte_zaehler_fun<-function(beta) {
		lambda<-exp(beta)
		sigma<-k_est*lambda # *sqrt(zeit)
		sigma2<-sigma^2
		mu<-(lambda-sigma2/2)*z
		#c(prod(dnorm(y,mu,sigma))*dnorm(beta,theta,nu)
		exp(sum(dnorm(r,mu,sigma*sqrt(z),log=TRUE))+dnorm(beta,theta_est,nu_est,log=TRUE))
	}

	LOW<- x0-M*nu_est
	UP<- x0+M*nu_est
	bed_dichte_zaehler_fun(x0)
	bed_dichte_zaehler_fun(x0-nu_est)
	bed_dichte_zaehler_fun(x0+nu_est)

	s<-seq(LOW,UP,length.out=nn)
	w<-mean(diff(s))
	f_s<-sapply(s,bed_dichte_zaehler_fun)
		
	intfun<-function(beta,q) {
		lambda<-exp(beta)
		sigma<-k_est*lambda #hier war in 1-3 noch *sqrt(zeit), gehört aber nicht hierher, steht jetzt in pnorm.
		sigma2<-sigma^2
		#pnorm(q,(lambda-sigma2/2)*zeit,sigma*sqrt(zeit))*dnorm(beta,theta_est,nu_est)
		#bedingte_dichte<-dnorm(beta,theta_est,nu_est)
		pnorm(q,(lambda-sigma2/2)*zeit,sigma*sqrt(zeit))*bed_dichte_zaehler_fun(beta)
		#dnorm(beta,theta_est,nu_est)

	}
	g_s<-sapply(s,intfun,q=q)
	sum(g_s*w)/sum(f_s*w)
}

