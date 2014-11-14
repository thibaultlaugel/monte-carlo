#############################################################################################################
#                                VARIANCE REDUCTION USING ROMBERG'S METHOD                                  #
#                                                                                                           #
#                            Thibault LAUGEL - Leonor CHABANNES - Chloe LACOMBE                             #
#############################################################################################################





############################### Third part: Applications



MouvBrownien1 = function(n){
  Normk =rnorm(n,sd=sqrt(1/n)) # k simulations de loi normale (0,t(i)-t(i-1))
  VecBrownien = cumsum(Normk) # cumulated sums vector =  n Brownian motions
  VecBrownien
}
MouvBrownien1(100)


#Simulation of N Brownian motion vectors of size n
#output: matrix N x n
MatriceBrownienMC = function(N,n){
  mat = sapply(rep(n,length=N),MouvBrownien1)
  return(t(mat))
}
MatriceBrownienMC(100,10)




############################ 1. Theoritical example



################################################
# Classic Monte Carlo method
################################################



# Euler's method for discretization:
SchemaEuler1=function(X0,Y0,B){
	l=length(B)
	X=c(1:l)
	Y=c(1:l)
	X[1]=X0-X0*1/(2*l)-Y0*(B[1]-0)
	Y[1]=Y0-Y0*1/(2*l)+X0*(B[1]-0)
	for(i in 1:(l-1)){
		X[i+1]=X[i]-X[i]*1/(2*l)-Y[i]*(B[i+1]-B[i])
		Y[i+1]=Y[i]-Y[i]*1/(2*l)+X[i]*(B[i+1]-B[i])
}		
	X=c(X0,X)
	Y=c(Y0,Y)
	return (matrix(c(X,Y),ncol=2,byrow=F,dimnames=list(c(0:l),c("X","Y"))))	
}
SchemaEuler1(0.5,0.5,MouvBrownien1(100))
##Construction of the process Z discretized, of Brownian motion B and initial value (X0,Y0)



# Definition of the function f_alpha
f=function(x,alpha){
	(	abs((x[1]^2) + (x[2]^2)-1)	)^(2*alpha)+x[1]
}


# Approximation :
MC1= function(n,Z0){
		alpha=0.5
		NombreSimu=n^(2*alpha)
		B=MatriceBrownienMC(NombreSimu,n)
		Z=apply(B,1,function(x,y,B) SchemaEuler1(x,y,B),x=Z0[1],y=Z0[2])
		XT=Z[dim(Z)[1]/2,]
		YT=Z[dim(Z)[1],]		
		ZT=cbind(XT,YT)
		F=apply(ZT,1,f,alpha=0.5)
		return (mean(F))
}
MC1(20,c(0.5,0.5))
## Construction, using MC method, of the expectancy of the image of our process Z under f_alpha (with alpha=0.5)



param1=function(M){
	theta=runif(M,min=0, max=2*pi)
	return(cbind(cos(theta),sin(theta)))
}
param1(20)
##Simulation of M values for Z0

#To evaluate the complexity of the alorithm, we measure, for a certain values of (n,M) the RMS and calculation time of the approximation by Monte Carlo
Calcul_RMS1=function(n,M){
		Z0=param1(M)
		t1=Sys.time()
		Estim=apply(Z0,1,function(x,Z0) MC1(x,Z0),x=n)
		t2=Sys.time()
		dt=as.numeric(t2-t1, units="secs")
		rms=sqrt(mean((Estim-exp(-0.5)*Z0[,1])^2))
		return(c(dt,rms))
	
}

Calcul_RMS1(10,3)


#plot of the curve (RMS(n),temps de calcul(n)), with 1/n the discretization stem
Graphe1=function(M,vecn){
	Z0=param1(M)
	err=NULL
	time=NULL
	for(i in vecn){
		NombreSimu=i
		tmp=Calcul_RMS1(i,M)
		err=c(err, tmp[2])
		time=c(time, tmp[1])
		}
	speed=M/time
	res=matrix(c(err,speed),ncol=2,byrow=F,dimnames=list(vecn,c("erreur","vitesse")))
	return(res)
}
Graphe1(3,c(10,20))




################################################
#	Romberg's Statistic Method
################################################


#Brownian motions:
MatriceBrownienRS=function (n,m,NombreSimu){
	t1=matrix(c(seq(1/n,1,by=1/n),rep(1,n)),n,2) 
	t2=matrix(c(seq(1/m,1,by=1/m),rep(2,m)),m,2)
	time=rbind(t1,t2)
	A=time[order(time[,1]),]
	dt=diff(c(0,A[,1]))
	Normales=matrix(rnorm((n+m)*NombreSimu,0,sd=sqrt(rep(dt,NombreSimu))),nrow=NombreSimu,ncol=n+m)
	Brown=t(apply(Normales,1,cumsum))
	B=rbind(Brown,t(A)) 
	Bsmall=B[,B[dim(B)[1],]<2] 
	Bsmall=Bsmall[-c(dim(Bsmall)[1]-1,dim(Bsmall)[1]),]
	Bbig=B[,B[dim(B)[1],]>1]
	Bbig=Bbig[-c(dim(Bbig)[1]-1,dim(Bbig)[1]),]
	return (list(Bsmall,Bbig))
}
##Simulation of NombreSimu Brownian motions, of steps 1/n and 1/m.
MatriceBrownienRS(10,10,10)


#Discretization of the process:
SchemaEuler2=function(X0,Y0,B){
	B1=B[[1]]
	B2=B[[2]]
	Xn=apply(B1,1,function(x,y,B1) SchemaEuler1(x,y,B1),x=X0,y=Y0)
	Xm=apply(B2,1,function(x,y,B2) SchemaEuler1(x,y,B2),x=X0,y=Y0)
	return (list(Xn,Xm))
}
SchemaEuler2(0.5,0.5,MatriceBrownienRS(10,10,10)
)



MCRS1=function(n,m,NombreSimu,X0,Y0){
	B=MatriceBrownienRS(n,m,NombreSimu)
	Z=SchemaEuler2(X0,Y0,B)
	XnT=Z[[1]][n+1,]
	YnT=Z[[1]][2*(n+1),]
	XmT=Z[[2]][m+1,]
	YmT=Z[[2]][2*(m+1),]
	ZnT=t(rbind(XnT,YnT))
	ZmT=t(rbind(XmT,YmT))
	t1=apply(ZnT, 1, f ,alpha=0.5)
	t2=apply(ZmT, 1 ,f, alpha=0.5)
	dt=t1-t2
	return(mean(dt))
}
MCRS1(10,10,10,0.5,0.5)


MC2= function(n,NombreSimu,Z0){
	B=MatriceBrownienMC(NombreSimu,n)
	Z=apply(B,1,function(x,y,B) SchemaEuler1(x,y,B),x=Z0[1],y=Z0[2])
	XT=Z[dim(Z)[1]/2,]
	YT=Z[dim(Z)[1],]		
	tmp=cbind(XT,YT)
	time=apply(tmp,1,f,alpha=0.5)
	return (mean(time))
}
MC2(10,10,param1(3))


RS=function(n,Z0){ 
	alpha=0.5
	X0=Z0[1]
	Y0=Z0[2]
	m=floor(n^0.5)
	Nm=floor(n^(2*alpha))
	Nn=floor(n^(2*alpha-0.5))
	m1=MC2(m,Nm,Z0)
	m2=MCRS1(n,m,Nn,X0,Y0)
	return (m1+m2)
}
RS(10,c(0.5,0.5))


#Complexity of the algorithm: RMS and calculation time:

Calcul_RMS2=function(n,M){
	Z0=param1(M)
	t1=Sys.time()
	Estim=apply(Z0,1,function(x,Z0) RS(x,Z0),x=n)
	t2=Sys.time()
	dt=as.numeric(t2-t1, units="secs")
	rms=sqrt(mean((Estim-exp(-0.5)*Z0[,1])^2))
	return(c(dt,rms))	
}

Calcul_RMS2(10,10)


#plot of the curve (RMS(n),temps de calcul(n))
Graphe2=function(M,vecn){
	err=NULL
	time=NULL
	for(i in vecn){
		NombreSimu=i
		tmp=Calcul_RMS2(i,M)
		err=c(err, tmp[2])
		time=c(time, tmp[1])
		}
	vitesse=M/time
	A=matrix(c(err,vitesse),ncol=2,byrow=F,dimnames=list(vecn,c("err","vitesse")))
	return(A)
}
Graphe2(10,c(10,20))


################################################
#			Plots
################################################

M=100

vecn=c(60,90,120,180,220,250)
mc=Graphe1(M,vecn)
rs=Graphe2(M,vecn)

plot(c(min(mc[,1],rs[,1]),max(mc[,1],rs[,1])),c(min(mc[,2],rs[,2]),max(mc[,2],rs[,2])),type="b",col="white",log="xy",xlab="RMS erreur",ylab="Vitesse")
lines(rs[,1],rs[,2],type="b")
lines(mc[,1],mc[,2],type="b",lty=2)









############################### 2. Application to pricing of Asian Options




#################################################
#	Classic Monte Carlo Method
################################################




CallAsiatique=function(B, S0, mu, sigma,  T){  
  n=length(B)
  t=seq(1/n, 1, by=1/n)*T #vector of the discretization steps
  
  S=S0*exp((mu-0.5*sigma^2)*t+sigma*B)
  S=c(S0,S[1:n-1]) #vector of the underlying hedged item prices
  
  Calculation of I
  vecmoinsun=B[1:n-1] # Brownian motions (size p-1)
  vec=B[2:n] # vector of the Brownien Motions indexed+1 (soze p-1)
  difference=c(B[1],vec-vecmoinsun) #initialization at zero, then difference to obtain B(t)-B(t-1)
  Trapeze= S * ( 1 + (mu*T)/(2*n) + sigma*difference/2 ) # trapezoidal scheme
  I=mean(Trapeze)
  return (exp(-mu*T)*max(0,I)) #estimation for N simulations, with n steps for the trapezoidal method
}
CallAsiatique(MouvBrownien1(10),10,0.5,0.2,10)


#estimation by Monte carlo, with m the step for the trapezoidal method, and Nm simulations

monteCarlo= function(m, Nm, S0, mu, sigma,  T){
	B=MatriceBrownienMC(Nm,m)*sqrt(T)
	res=apply(B,1,CallAsiatique, S0=S0, mu=mu, sigma=sigma,  T=T)
	M=mean(res)

	return (M) 
}
monteCarlo(10,10,10,0.5,0.2,10)




################################################
#	Romberg's Statistic Method
################################################


MatriceBrownienRomberg=function (n,m,N){

  tempsn=matrix(c(seq(1/n,1,by=1/n),rep(1,n)),n,2) 
  tempsm=matrix(c(seq(1/m,1,by=1/m),rep(2,m)),m,2)
  temps=rbind(tempsn,tempsm)

  
  triage=temps[order(temps[,1]),] #Sort using the order of the steps
  pas=diff(c(0,triage[,1])) #Differenciation
  
  Norm=t(matrix(rnorm((n+m)*N,0,sd=sqrt(rep(pas,N))),nrow=n+m,ncol=N))
  B=t(apply(Norm,1,cumsum)) # N Brownian motions of size n+m  
  Btriage=rbind(B,t(triage)) #To identify steps of the BM (1/m or 1/n)
  
  MatBrown_n=Btriage[,Btriage[dim(Btriage)[1],]<2] 
  MatBrown_n=MatBrown_n[-c(dim(MatBrown_n)[1]-1,dim(MatBrown_n)[1]),] #BM step 1/n
  MatBrown_m=Btriage[,Btriage[dim(Btriage)[1],]>1]
  MatBrown_m=MatBrown_m[-c(dim(MatBrown_m)[1]-1,dim(MatBrown_m)[1]),] #BM step 1/m
  return (list(MatBrown_n,MatBrown_m))
 
}
MatriceBrownienRomberg(10,10,10)



#First part of En
monteCarlo_pour_Romberg= function(n, m, Nn, S0, mu, sigma, T){
	B=MatriceBrownienRomberg(n,m, Nn)
	B[[1]]=B[[1]]*sqrt(T) #BM (n)
	B[[2]]=B[[2]]*sqrt(T) #BM (m)

	Q1=apply(B[[1]],1,CallAsiatique, S0=S0, mu=mu, sigma=sigma, T=T) #Value of the call ,  n
	Q2=apply(B[[2]],1,CallAsiatique, S0=S0, mu=mu, sigma=sigma, T=T)#value of the call, m

	E=Q1-Q2
	res=mean(E)
	return (res) 
}
monteCarlo_pour_Romberg(101,50,100,35,0.1,0.1,10)


Estimation_Romberg=function(err, S0, mu, sigma, T){
  n=floor(1/err)+1 # Error 1/n, floor to avoid null values
  m=floor(n^(1/3))+1
  Nn=floor(n^(4/3))+1
  Nm=n^2
	E1=monteCarlo(m, Nm, S0, mu, sigma, T)
	E2=monteCarlo_pour_Romberg(n, m, Nn, S0, mu, sigma, T)
	res=E1+E2
	return (res)
}
Estimation_Romberg(0.01,10,0.1,0.1,2)


#real value of the Asian call
valeur=function(S0,mu,T){
  res=(S0/(mu*T))*(1-exp(-mu*T))
  return(res)
}

#Generation of the parameters
param2=function(M){
  
  S0=runif(M, min=20, max=40)
  mu=runif(M, min=0.01, max=0.15)
  sigma=runif(M, min=0.05, max=0.7)
  T=runif(M, min=0.05, max=5)
  valeursreelles=mapply(valeur, S0=S0,mu=mu,T=T)

  return (cbind(S0, mu, sigma, T, valeursreelles))  
}
param2(50)

Calcul_RMS_Romberg=function(n, M){	
	#Récuperation of the parameters
	parame=param2(M)	
	S0=parame[,1]
	mu=parame[,2]
	sigma=parame[,3]
	Temps=parame[,4]
	valeursreelles=parame[,5]

	err=array(1/n, M)

#Calculation time
	T1=Sys.time()
	valeursSimulees=mapply(Estimation_Romberg, err=err, S0=S0, mu=mu, sigma=sigma, T=Temps)
	T2=Sys.time()
	time=as.numeric(T2-T1, units="secs")

#Calculation of the RMS
	RMS=sqrt(mean((valeursSimulees-valeursreelles)^2	))
	return(c(time,RMS))
}
Calcul_RMS_Romberg(10,20)

Calcul_RMS_Montecarlo=function(n, M){

#Parameters
	param=param2(M)
	S0=param[,1]
	mu=param[,2]
	sigma=param[,3]
	T=param[,4]
	valeursreelles=param[,5]

	m=array(n,M)
	Nm=m^2
#Calculation time
	T1=Sys.time()
	valeursSimulees=mapply(monteCarlo, m=m, Nm=Nm, S0=S0, mu=mu, sigma=sigma, T=T)
	T2=Sys.time()
	temps=as.numeric(T2-T1, units="secs")

#Calculation of the RMS		
	RMS=sqrt(mean((valeursSimulees-valeursreelles)^2	))
	return(c(temps,RMS))
}

Calcul_RMS_Montecarlo(10,50)

#Calculation of the calculation time and RMS for several values of n, for MC
FinalMC=function(M,vecn){
	#initialisation	
	err=NULL
	temps=NULL
	for(i in vecn){
		temp=Calcul_RMS_Montecarlo(i,M)
		err=c(err,temp[2])
		temps=c(temps,temp[1])
	}
	vitesse=M/temps
	res=matrix(c(err,vitesse, temps),ncol=3,byrow=F,dimnames=list(vecn,c("erreur","vitesse","temps")))
	return(res)
}
FinalMC(50,c(10, 20))	

#Idem for Romgberg

FinalRomberg=function(M,vecn){
	#initialisation		
	err=NULL
	temps=NULL
	for(i in vecn){
		temp=Calcul_RMS_Romberg(i,M)
		err=c(err,temp[2])
		temps=c(temps,temp[1])
	}
	vitesse=M/temps
	res=matrix(c(err,vitesse, temps),ncol=3,byrow=F,dimnames=list(vecn,c("erreur","vitesse","temps")))
	return(res)

}
FinalRomberg(50,c(10, 20))  





################################################
#			Plots
################################################


#Final values
a=FinalRomberg(50,c(10, 20,60,120,180))	
b=FinalMC(50,c(10, 20,60,120,180))	


#Plots
plot(c(min(b[,1],a[,1]),max(b[,1],a[,1])),c(min(b[,2],a[,2]),max(b[,2],a[,2])),type="b",col="white",log="xy",xlab="RMS erreur",ylab="Vitesse")
lines(a[,1],a[,2],type="b")
lines(b[,1],b[,2],type="b",lty=2)










###########################Using Sobol sequence

library(fOptions)
library(randtoolbox)

MouvBrownien1 = function(n){
  Normk =sobol(n , dim = 1, scrambling = 3, normal = TRUE) * sqrt(1/n)
 VecBrownien = cumsum(Normk) # cumulated sums vector
  VecBrownien
}
MouvBrownien1(100)

sobol(n = 10, dim = 1, scrambling = 3, normal = TRUE)

#Simulation of N BM
#output is a matrix N x n
MatriceBrownienMC = function(N,n){
  mat = sapply(rep(n,length=N),MouvBrownien1)
  return(t(mat))
}
MatriceBrownienMC(100,20)




############################ 1. Theoritical Example

################################################
#	Classic Monte Carlo Method
################################################



#Euler Method
SchemaEuler1=function(X0,Y0,B){
	l=length(B)
	X=c(1:l)
	Y=c(1:l)
	X[1]=X0-X0*1/(2*l)-Y0*(B[1]-0)
	Y[1]=Y0-Y0*1/(2*l)+X0*(B[1]-0)
	for(i in 1:(l-1)){
		X[i+1]=X[i]-X[i]*1/(2*l)-Y[i]*(B[i+1]-B[i])
		Y[i+1]=Y[i]-Y[i]*1/(2*l)+X[i]*(B[i+1]-B[i])
}		
	X=c(X0,X)
	Y=c(Y0,Y)
	return (matrix(c(X,Y),ncol=2,byrow=F,dimnames=list(c(0:l),c("X","Y"))))	
}
SchemaEuler1(0.5,0.5,MatriceBrownienMC(100,20)
)


f=function(x,alpha){
	(	abs((x[1]^2) + (x[2]^2)-1)	)^(2*alpha)+x[1]
}


MC1= function(n,Z0){
		alpha=0.5
		NombreSimu=n^(2*alpha)
		B=MatriceBrownienMC(NombreSimu,n)
		Z=apply(B,1,function(x,y,B) SchemaEuler1(x,y,B),x=Z0[1],y=Z0[2])
		XT=Z[dim(Z)[1]/2,]
		YT=Z[dim(Z)[1],]		
		ZT=cbind(XT,YT)
		F=apply(ZT,1,f,alpha=0.5)
		return (mean(F))
}
MC1(20,c(0.5,0.5))


param1=function(M){
	theta=runif(M,min=0, max=2*pi)
	return(cbind(cos(theta),sin(theta)))
}
param1(50)

#Complexity
Calcul_RMS1=function(n,M){
		Z0=param1(M)
		t1=Sys.time()
		Estim=apply(Z0,1,function(x,Z0) MC1(x,Z0),x=n)
		t2=Sys.time()
		dt=as.numeric(t2-t1, units="secs")
		rms=sqrt(mean((Estim-exp(-0.5)*Z0[,1])^2))
		return(c(dt,rms))
	
}

Calcul_RMS1(10,3)


Graphe1=function(M,vecn){
	Z0=param1(M)
	err=NULL
	time=NULL
	for(i in vecn){
		NombreSimu=i
		tmp=Calcul_RMS1(i,M)
		err=c(err, tmp[2])
		time=c(time, tmp[1])
		}
	speed=M/time
	res=matrix(c(err,speed),ncol=2,byrow=F,dimnames=list(vecn,c("erreur","vitesse")))
	return(res)
}
Graphe1(3,c(10,20))




################################################
#	Romberg's Statistic Method
################################################


MatriceBrownienRS=function (n,m,NombreSimu){
	t1=matrix(c(seq(1/n,1,by=1/n),rep(1,n)),n,2) 
	t2=matrix(c(seq(1/m,1,by=1/m),rep(2,m)),m,2)
	time=rbind(t1,t2)
	A=time[order(time[,1]),]
	dt=diff(c(0,A[,1]))
	Normales=matrix(rnorm((n+m)*NombreSimu,0,sd=sqrt(rep(dt,NombreSimu))),nrow=NombreSimu,ncol=n+m)
	Brown=t(apply(Normales,1,cumsum))
	B=rbind(Brown,t(A)) 
	Bsmall=B[,B[dim(B)[1],]<2] 
	Bsmall=Bsmall[-c(dim(Bsmall)[1]-1,dim(Bsmall)[1]),]
	Bbig=B[,B[dim(B)[1],]>1]
	Bbig=Bbig[-c(dim(Bbig)[1]-1,dim(Bbig)[1]),]
	return (list(Bsmall,Bbig))
}
MatriceBrownienRS(10,10,10)


SchemaEuler2=function(X0,Y0,B){
	B1=B[[1]]
	B2=B[[2]]
	Xn=apply(B1,1,function(x,y,B1) SchemaEuler1(x,y,B1),x=X0,y=Y0)
	Xm=apply(B2,1,function(x,y,B2) SchemaEuler1(x,y,B2),x=X0,y=Y0)
	return (list(Xn,Xm))
}
SchemaEuler2(0.5,0.5,MatriceBrownienRS(10,10,10)
)



MCRS1=function(n,m,NombreSimu,X0,Y0){
	B=MatriceBrownienRS(n,m,NombreSimu)
	Z=SchemaEuler2(X0,Y0,B)
	XnT=Z[[1]][n+1,]
	YnT=Z[[1]][2*(n+1),]
	XmT=Z[[2]][m+1,]
	YmT=Z[[2]][2*(m+1),]
	ZnT=t(rbind(XnT,YnT))
	ZmT=t(rbind(XmT,YmT))
	t1=apply(ZnT, 1, f ,alpha=0.5)
	t2=apply(ZmT, 1 ,f, alpha=0.5)
	dt=t1-t2
	return(mean(dt))
}
MCRS1(10,10,10,0.5,0.5)


MC2= function(n,NombreSimu,Z0){
	B=MatriceBrownienMC(NombreSimu,n)
	Z=apply(B,1,function(x,y,B) SchemaEuler1(x,y,B),x=Z0[1],y=Z0[2])
	XT=Z[dim(Z)[1]/2,]
	YT=Z[dim(Z)[1],]		
	tmp=cbind(XT,YT)
	time=apply(tmp,1,f,alpha=0.5)
	return (mean(time))
}
MC2(10,10,param1(3))


RS=function(n,Z0){ 
	alpha=0.5
	X0=Z0[1]
	Y0=Z0[2]
	m=floor(n^0.5)
	Nm=floor(n^(2*alpha))
	Nn=floor(n^(2*alpha-0.5))
	m1=MC2(m,Nm,Z0)
	m2=MCRS1(n,m,Nn,X0,Y0)
	return (m1+m2)
}
RS(10,c(0.5,0.5))


Calcul_RMS2=function(n,M){
	Z0=param1(M)
	t1=Sys.time()
	Estim=apply(Z0,1,function(x,Z0) RS(x,Z0),x=n)
	t2=Sys.time()
	dt=as.numeric(t2-t1, units="secs")
	rms=sqrt(mean((Estim-exp(-0.5)*Z0[,1])^2))
	return(c(dt,rms))	
}
Calcul_RMS2(10,10)


Graphe2=function(M,vecn){
	
	err=NULL
	time=NULL
	for(i in vecn){
		NombreSimu=i
		tmp=Calcul_RMS2(i,M)
		err=c(err, tmp[2])
		time=c(time, tmp[1])
		}
	
	vitesse=M/time
	A=matrix(c(err,vitesse),ncol=2,byrow=F,dimnames=list(vecn,c("err","vitesse")))
	return(A)
}
Graphe2(10,c(10,20))
##On a bien les vecteurs attendus. 


################################################
#			Plots
################################################

M=100


vecn=c(60,90,120,180,220,250)
mc=Graphe1(M,vecn)
rs=Graphe2(M,vecn)

plot(c(min(mc[,1],rs[,1]),max(mc[,1],rs[,1])),c(min(mc[,2],rs[,2]),max(mc[,2],rs[,2])),type="b",col="white",log="xy",xlab="RMS erreur",ylab="Vitesse")
lines(rs[,1],rs[,2],type="b")
lines(mc[,1],mc[,2],type="b",lty=2)







