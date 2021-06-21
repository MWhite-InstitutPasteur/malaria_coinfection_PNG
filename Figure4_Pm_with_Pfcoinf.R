library(VennDiagram)
library(binom)
library(fields)



ALB_par <- read.csv("Albinama_data_PAR.csv", sep=";")

ALB_par <- ALB_par[which(is.na(ALB_par$qmal_pos)==FALSE),]

N_data <- nrow(ALB_par)


ALB_epi <- read.csv("Albimana_data_EPI.csv")

village_names <- unique(ALB_epi$villagegroup)
N_village <- length(village_names)


ALB_par_village <- rep(NA, N_data)

for(i in 1:N_data)
{
	ALB_par_village[i] <- as.character(ALB_epi$villagegroup_corrected[which(ALB_epi$studyid == ALB_par$studyid[i])][1])
}

ALB_par$village <- ALB_par_village



ALB_par_treat <- rep(NA, N_data)

for(i in 1:N_data)
{
	ALB_par_treat[i] <- as.character(ALB_epi$treat[which(ALB_epi$studyid == ALB_par$studyid[i])][1])
}

#for(i in 1:N_data)
#{
#	if( ALB_par$study_day[i] < 0 )
#	{	
#		ALB_par_treat[i] <- "NONE"
#	}
#}

ALB_par$treat <- ALB_par_treat



tt <- c(-30, 0, 14, 28, 42, 56, 70, 84, 112, 140, 168, 196, 224)  

tt_edges <- c(-50, -5, 7, 21, 35, 49, 63, 77, 98, 126, 154, 182, 210, 240)

N_tt <- length(tt)   


###############
###############
##           ##
##  PANEL 1  ##
##           ##
###############
###############

#############################################
## Asexuals; Placebo; Prevalence

Pm_PLPQ <- matrix(NA, nrow=N_tt, ncol=3)
colnames(Pm_PLPQ) <- c("med", "lwr", "upr")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	x_i <- ALB_par$pm_pos[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pm_PLPQ[i,] <- as.numeric(binom.confint( sum(x_i), length(x_i), method="wilson")[1,4:6])
}


#############################################
## Asexuals; Placebo; Prevalence co-infection

Pm_PLPQ_coinf <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pm_PLPQ_coinf) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	x_i <- ALB_par$pm_pos[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	N_i <- length(x_i)

	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 1 | 
                                ALB_par$pfg_pos[index_i] == 1  ) ]

	x_i <- ALB_par$pm_pos[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pm_PLPQ_coinf[i,1:3] <- as.numeric(binom.confint( sum(x_i), N_i, method="wilson")[1,4:6])
	Pm_PLPQ_coinf[i,4] <- N_i
}


###############
###############
##           ##
##  PANEL 2  ##
##           ##
###############
###############

#############################################
## Gametocytes; Placebo; Prevalence

Pmg_PLPQ <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmg_PLPQ) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	x_i <- ALB_par$pmg_pos[index_i]

	Pmg_PLPQ[i,1:3] <- as.numeric(binom.confint( sum(x_i, na.rm=TRUE), length(x_i), method="wilson")[1,4:6])

	Pmg_PLPQ[i,4] <- length(x_i) 
}



#############################################
## Gametocytes; Placebo; Prevalence co-infection

Pmg_PLPQ_coinf <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmg_PLPQ_coinf) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	N_i <- length(index_i)

	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 1 | 
                                ALB_par$pfg_pos[index_i] == 1  ) ]

	x_i <- ALB_par$pmg_pos[index_i]

	Pmg_PLPQ_coinf[i,1:3] <- as.numeric(binom.confint( sum(x_i, na.rm=TRUE), N_i, method="wilson")[1,4:6])
	Pmg_PLPQ_coinf[i,4] <- N_i 
}





###############
###############
##           ##
##  PANEL 5  ##
##           ##
###############
###############


#############################################
## Asexuals; Placebo; Density

Pmcopy_PLPQ <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmcopy_PLPQ) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] & 
                        ALB_par$pm_pos==1 ) 	

	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 0 & 
                                ALB_par$pfg_pos[index_i] == 0  ) ]

	x_i <- ALB_par$pm_copyno[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pmcopy_PLPQ[i,1:3] <- quantile(x_i, prob=c(0.5, 0.025, 0.975))
	Pmcopy_PLPQ[i,4] <- length(x_i) 
}



#############################################
## Asexuals; Placebo; Density co-infection

Pmcopy_PLPQ_coinf <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmcopy_PLPQ_coinf) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] 
                         & ALB_par$pm_pos == 1 ) 	


	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 1 | 
                                ALB_par$pfg_pos[index_i] == 1  ) ]

	x_i <- ALB_par$pm_copyno[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pmcopy_PLPQ_coinf[i,1:3] <- quantile(x_i, prob=c(0.5, 0.025, 0.975))
	Pmcopy_PLPQ_coinf[i,4] <- length(x_i)
}



###############
###############
##           ##
##  PANEL 7  ##
##           ##
###############
###############



#############################################
## Gametocytes; Placebo; Density

Pmgcopy_PLPQ <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmgcopy_PLPQ) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] & 
                        ALB_par$pmg_pos==1 ) 	

	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 0 & 
                                ALB_par$pfg_pos[index_i] == 0  ) ]

	x_i <- ALB_par$pmg_copyno[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pmgcopy_PLPQ[i,1:3] <- quantile(x_i, prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
	Pmgcopy_PLPQ[i,4] <- length(x_i) 
}



#############################################
## Gametocytes; Placebo; Density co-infection

Pmgcopy_PLPQ_coinf <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pmgcopy_PLPQ_coinf) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] 
                         & ALB_par$pmg_pos == 1 ) 	


	index_i <- index_i[which( ALB_par$pf_pos[index_i] == 1 | 
                                ALB_par$pfg_pos[index_i] == 1  ) ]

	x_i <- ALB_par$pmg_copyno[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pmgcopy_PLPQ_coinf[i,1:3] <- quantile(x_i, prob=c(0.5, 0.025, 0.975))
	Pmgcopy_PLPQ_coinf[i,4] <- length(x_i) 
}




##################################
##################################
##                              ##
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ##
##  ##     ##    ##  ##   ##    ##
##  ##     #####  ####    ##    ##
##                              ##
##################################
##################################


Pf_col <- "red"
Pv_col <- "orange"
Pm_col <- "black"
Po_col <- "dodgerblue"

Coinf_col <- "red"


main.size = 1.3
axis.size = 0.7
lab.size  = 1.5

point.size = 1
line.size  = 1
rr.size1    = 0.005
rr.size2    = 0.0075
arrow.edge = 0.015
arrow.size = 0.5
text.size = 0.5


tiff( file="Figure4_Pm_with_Pfcoinf.tif", width=16, height=12, units="cm", res=500)


lay.mat = rbind( c( 1, 2 ),
                 c( 3, 3 ),
                 c( 4, 5 ), 
                 c( 6, 6 ) )

layout(lay.mat, heights=c(10,1,10,1))
layout.show(6)

par(mar = c(3,3,2,0.5))

par(mgp = c(1.85, 0.6, 0))



########################################
##                                    ##
##  PANEL 1: Asexual prevalence       ##
##          	                    ## 
########################################


line_seq_x <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)

line_seq_y <- c(-50, 0, 50, 100, 150, 200, 250)



plot(x=100, y=100, 
xlim=c(-50,250), ylim=c(0,0.2), 
xlab="time (days)", ylab=" prevalence",
main="(A) P. malariae asexual infection",
xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.lab=lab.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(-1e10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


polygon( x=c(0, 14, 14, 0), y=c(-1, -1, 1, 1),
         col=rgb(229/256, 229/256, 229/256, 0.5), border=NA )

points(x=tt-1, y=Pm_PLPQ[,1], 
type='l', lwd=line.size, col=Pm_col)

points(x=tt-1, y=Pm_PLPQ[,1], 
pch=19, cex=point.size, col=Pm_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]-1, y0=Pm_PLPQ[i,2], 
     		 x1=tt[i]-1, y1=Pm_PLPQ[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}



points(x=tt+1, y=Pm_PLPQ_coinf[,1], 
type='l', lwd=line.size, col=Coinf_col)

points(x=tt+1, y=Pm_PLPQ_coinf[,1], 
pch=17, cex=point.size, col=Coinf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pm_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pm_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Coinf_col, lwd=arrow.size)	
}





axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=axis.size )


axis(2, las=2, at = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2), 
     labels = c("0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%", "20%"),
cex.axis=axis.size )


########################################
##                                    ##
##  PANEL 2: Gametocyte prevalence    ##
##          	                    ## 
########################################


line_seq_x <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12)

line_seq_y <- c(-50, 0, 50, 100, 150, 200, 250)



plot(x=100, y=100, 
xlim=c(-50,250), ylim=c(0,0.12), 
xlab="time (days)", ylab=" prevalence",
main="(B) P. malariae gametocyte infection",
xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.lab=lab.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(-1e10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


polygon( x=c(0, 14, 14, 0), y=c(-1, -1, 1, 1),
         col=rgb(229/256, 229/256, 229/256, 0.5), border=NA )

points(x=tt-1, y=Pmg_PLPQ[,1], 
type='l', lwd=line.size, col=Pm_col)

points(x=tt-1, y=Pmg_PLPQ[,1], 
pch=19, cex=point.size, col=Pm_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]-1, y0=Pmg_PLPQ[i,2], 
     		 x1=tt[i]-1, y1=Pmg_PLPQ[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}



points(x=tt+1, y=Pmg_PLPQ_coinf[,1], 
type='l', lwd=line.size, col=Coinf_col)

points(x=tt+1, y=Pmg_PLPQ_coinf[,1], 
pch=17, cex=point.size, col=Coinf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pmg_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pmg_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Coinf_col, lwd=arrow.size)	
}





axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=axis.size )



axis(2, las=2, at = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12), 
     labels = c("0%", "2%", "4%", "6%", "8%", "10%", "12%"),
cex.axis=axis.size )






########################################
##                                    ##
##  PANEL 3: Legend                   ##
##          	                    ## 
########################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("P. malariae (all infections)", "P. malariae & P. falciparum co-infection"),
       #fill = c(Pm_col, Coinf_col), 
       #border = c(Pm_col, Coinf_col), 
       col = c(Pm_col, Coinf_col), 
       pch=c(19,17), ncol=2, cex=1.25, bty="n" )


Pm_col <- "forestgreen"

########################################
##                                    ##
##  PANEL 4: Asexual density          ##
##          	                    ## 
########################################


par(mar = c(3,3,2,0.5))

par(mgp = c(1.85, 0.6, 0))


line_seq_x <- c(1, 3, 10, 30, 100, 300, 1000, 3000)

line_seq_y <- c(-50, 0, 50, 100, 150, 200, 250)



plot(x=10000, y=10000, 
xlim=c(-50,250), ylim=c(1,3000), log="y",
xlab="time (days)", ylab="copy number",
main="(C) P. malariae asexual density",
xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.lab=lab.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


polygon( x=c(0, 14, 14, 0), y=c(0.001, 0.001, 10000, 10000),
         col=rgb(229/256, 229/256, 229/256, 0.5), border=NA )

points(x=tt-1, y=Pmcopy_PLPQ[,1], 
type='l', lwd=line.size, col=Pm_col)

points(x=tt-1, y=Pmcopy_PLPQ[,1], 
pch=19, cex=point.size, col=Pm_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]-1, y0=Pmcopy_PLPQ[i,2], 
     		 x1=tt[i]-1, y1=Pmcopy_PLPQ[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}



points(x=tt+1, y=Pmcopy_PLPQ_coinf[,1], 
type='l', lwd=line.size, col=Coinf_col)

points(x=tt+1, y=Pmcopy_PLPQ_coinf[,1], 
pch=17, cex=point.size, col=Coinf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pmcopy_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pmcopy_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Coinf_col, lwd=arrow.size)	
}


for(i in 1:N_tt)
{
	if( Pmcopy_PLPQ[i,4] > 0 )
	{	
		text(x=tt[i]-4, y=2000, 
	     	labels = Pmcopy_PLPQ[i,4], cex=text.size, col=Pm_col)	
	}
}


for(i in 1:N_tt)
{
	if( Pmcopy_PLPQ_coinf[i,4] > 0 )
	{
		text(x=tt[i]+4, y=2000, 
     		labels = Pmcopy_PLPQ_coinf[i,4], cex=text.size, col=Coinf_col)	
	}
}



axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=axis.size )


axis(2, las=2, at = c(1, 3, 10, 30, 100, 300, 1000, 3000), 
     labels = c(1, 3, 10, 30, 100, 300, 1000, 3000),
cex.axis=axis.size )






########################################
##                                    ##
##  PANEL 5: Gametocyte density       ##
##          	                    ## 
########################################


line_seq_x <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300)

line_seq_y <- c(-50, 0, 50, 100, 150, 200, 250)



plot(x=10000, y=10000, 
xlim=c(-50,250), ylim=c(0.01,300), log="y",
xlab="time (days)", ylab="copy number",
main="(D) P. malariae gametocyte density",
xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.lab=lab.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


polygon( x=c(0, 14, 14, 0), y=c(0.001, 0.001, 10000, 10000),
         col=rgb(229/256, 229/256, 229/256, 0.5), border=NA )

#points(x=tt-1, y=Pmgcopy_PLPQ[,1], 
#type='l', lwd=line.size, col=Pm_col)

points(x=tt-1, y=Pmgcopy_PLPQ[,1], 
pch=19, cex=point.size, col=Pm_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]-1, y0=Pmgcopy_PLPQ[i,2], 
     		 x1=tt[i]-1, y1=Pmgcopy_PLPQ[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}



#points(x=tt+1, y=Pmgcopy_PLPQ_coinf[,1], 
#type='l', lwd=line.size, col=Coinf_col)

points(x=tt+1, y=Pmgcopy_PLPQ_coinf[,1], 
pch=17, cex=point.size, col=Coinf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pmgcopy_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pmgcopy_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Coinf_col, lwd=arrow.size)	
}




for(i in 1:N_tt)
{
	if( Pmgcopy_PLPQ[i,4] > 0 )
	{	
		text(x=tt[i]-4, y=200, 
	     	labels = Pmgcopy_PLPQ[i,4], cex=text.size, col=Pm_col)	
	}
}


for(i in 1:N_tt)
{
	if( Pmgcopy_PLPQ_coinf[i,4] > 0 )
	{
		text(x=tt[i]+4, y=200, 
     		labels = Pmgcopy_PLPQ_coinf[i,4], cex=text.size, col=Coinf_col)	
	}
}



axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=axis.size )


axis(2, las=2, at = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300), 
     labels = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300),
cex.axis=axis.size )




########################################
##                                    ##
##  PANEL 6: Legend                   ##
##          	                    ## 
########################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("P. malariae (without Pf co-infection)", "P. malariae (Pf co-infection only)"),
       #fill = c(Pm_col, Coinf_col), 
       #border = c(Pm_col, Coinf_col), 
       col = c(Pm_col, Coinf_col), 
       pch=c(19,17), ncol=2, cex=1.25, bty="n" )


dev.off()








