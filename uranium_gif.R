library("Rcpp")
sourceCpp("/Users/Mikhail/Documents/ProgrammingI/R/Uranium/uranium.cpp")

simulation <- function(times, folder, n=200, dt=0.01, bomb=0.1*n) {
	#times: number of pictures
	#folder: path to folder where images and gif will be created
		#format: /Users/Documents/Uranium/
		#Warning: last forward slash is necessary!
	#n: number of particles
	#dt: time in particles' system between two successive pictures
	#bomb: lowest limit of mass for particle to explode
	
	l<-10;k<-2;border<-l+0.1*l;limit.q<-l/2;limit.v<-l;E<-0.5
	#l: half-length of box
	#k: factor by which size of particle is reduces
	#border: half-length of external border of box
	#limit.q: limit of initial coordinate
	#limit.v: limit of initial velocity
	#E: increase in velocity of particles after explosion
	
	#creating space of particles: ($x - $vy: coordinates and corresponding velocities)
	#							  (m: mass; e (exist?): T if particles exists, F otherwise)
	space <- list(x=runif(n,-limit.q,limit.q),y=runif(n,-limit.q,limit.q),vx=runif(n,-limit.v,limit.v),vy=runif(n,-limit.v,limit.v),m=rep(1,times=n),e=rep(T,times=n))

	#fill number with zeroes to width 6
	extend <- function(a) {
		n<-length(strsplit(as.character(a),split="")[[1]])
		zeroes<-6-n
		zeroes<-paste(rep(0,times=zeroes),collapse="")
		res<-paste(zeroes,a,sep="")
		return(res)
	}
	
	#draw layout
	draw_plot<- function() {
		plot(x=c(-border,border),y=c(-border,border),type="n",xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
		polygon(x=c(-border,-border,border,border),y=c(-border,-l,-l,-border),col=rgb(201,254,255,maxColorValue=255),border=NA)
		polygon(x=c(-border,-border,border,border),y=c(border,l,l,border),col=rgb(201,254,255,maxColorValue=255),border=NA)
		polygon(x=c(-border,-border,-l,-l),y=c(-l,l,l,-l),col=rgb(201,254,255,maxColorValue=255),border=NA)
		polygon(x=c(l,l,border,border),y=c(-l,l,l,-l),col=rgb(201,254,255,maxColorValue=255),border=NA)
		abline(v=c(-l,l),h=c(-l,l))
	}
	
	#progress br for uploading photos
	progress_bar <- txtProgressBar(min=0,max=times,style=1,char="=")
	
	#creating "myplotXXXXXX.png"s in /Users/Name/Documents/Uranium/test/

	attach(space)
	t<-0
	while(t<times){
		#creating "myplotXXXXXX.png"
		png(filename=paste(folder,"myplot",extend(t),".png",sep=""),width=1000,height=1000)		
		#draw to .png
		draw_plot()
		motion(x,y,vx,vy,l,dt,m,e,bomb,E,k)
		points(x=x[e],y=y[e],pch=16,cex=m[e]/k,col=rgb(0,m[e],0,bomb/1.5,maxColorValue=bomb))
		dev.off()
		t<-t+1
		setTxtProgressBar(progress_bar,value=t)
	}
	
	
	close(progress_bar)
	cat("Assembling gif...\n")
	#to be implementated
	#automatically create gif and delete pictures
	#system("convert -delay 5 *.png uranium.gif")
	#file.remove(list.files(pattern=".png"))
	cat("Finished\n")
	detach("space")
}