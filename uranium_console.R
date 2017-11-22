library("Rcpp")
sourceCpp("uranium.cpp")
simulation <- function(n=200, dt=0.01, sleep=0.07, bomb=0.1*n) {
	#n: number of particles
	#dt: time in particles' system between two successive frames
	#sleep: real time between two successive frames
	#in case of problem with performance increase sleep time
	#bomb: lowest limit of mass for particle to explode
	
	l<-10;k<-2;border<-l+0.1*l;limit.q<-l/2;limit.v<-l/2;E<-0.5
	#l: half-length of box
	#k: factor by which size of particle is reduces
	#border: half-length of external border of box
	#limit.q: limit of initial coordinate
	#limit.v: limit of initial velocity
	#E: increase in velocity of particles after explosion
	
	#creating space of particles: ($x - $vy: coordinates and corresponding velocities)
	#							  (m: mass; e (exist?): 1 if particles exists, 0 otherwise)
	space <- list(x=runif(n,-limit.q,limit.q),y=runif(n,-limit.q,limit.q),vx=runif(n,-limit.v,limit.v),vy=runif(n,-limit.v,limit.v),m=rep(1,times=n),e=rep(1,times=n))
	
	#draw layout
	plot(x=c(-border,border),y=c(-border,border),type="n",xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
	polygon(x=c(-border,-border,border,border),y=c(-border,-l,-l,-border),col=rgb(201,254,255,maxColorValue=255),border=NA)
	polygon(x=c(-border,-border,border,border),y=c(border,l,l,border),col=rgb(201,254,255,maxColorValue=255),border=NA)
	polygon(x=c(-border,-border,-l,-l),y=c(-l,l,l,-l),col=rgb(201,254,255,maxColorValue=255),border=NA)
	polygon(x=c(l,l,border,border),y=c(-l,l,l,-l),col=rgb(201,254,255,maxColorValue=255),border=NA)
	abline(v=c(-l,l),h=c(-l,l))
	
	#start
	while(T){
		t1<-Sys.time()
		motion(space$x,space$y,space$vx,space$vy,l,dt,space$m,space$e,bomb,E,k)
		abline(v=c(-l,l),h=c(-l,l))
		points(x=space$x,y=space$y,pch=16,cex=space$m/k,col=rgb(0,space$m,0,bomb/1.5,maxColorValue=bomb))
		while(Sys.time()-t1<sleep){}
		points(x=space$x,y=space$y,pch=16,col="white",cex=1.3*space$m/k)
	}
	return()
}
