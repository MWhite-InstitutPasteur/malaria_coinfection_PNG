
library(binom)




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



Pf_PLPQ <- matrix(NA, nrow=N_tt, ncol=3)
colnames(Pf_PLPQ) <- c("med", "lwr", "upr")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	x_i <- ALB_par$pf_pos[index_i]
	x_i <- x_i[which(is.na(x_i)==FALSE)]

	Pf_PLPQ[i,] <- as.numeric(binom.confint( sum(x_i), length(x_i), method="wilson")[1,4:6])
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


#############################################
## Asexuals; Proportion co-infected

Pm_Pcoinf <- matrix(NA, nrow=N_tt, ncol=3)
colnames(Pm_Pcoinf) <- c("med", "lwr", "upr")


for(i in 1:N_tt)
{
	Pm_Pcoinf[i,1:3] <- as.numeric(binom.confint( Pm_PLPQ_coinf[i,4]*Pm_PLPQ_coinf[i,1], 
                                                    Pm_PLPQ_coinf[i,4]*Pm_PLPQ[i,1], 
                                                    method="wilson")[1,4:6])
}





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



Pfg_PLPQ <- matrix(NA, nrow=N_tt, ncol=4)
colnames(Pfg_PLPQ) <- c("med", "lwr", "upr", "n")

for(i in 1:N_tt)
{
	index_i <- which( ALB_par$study_day > tt_edges[i] & ALB_par$study_day <= tt_edges[i+1] ) 	

	x_i <- ALB_par$pfg_pos[index_i]

	Pfg_PLPQ[i,1:3] <- as.numeric(binom.confint( sum(x_i, na.rm=TRUE), length(x_i), method="wilson")[1,4:6])

	Pfg_PLPQ[i,4] <- length(x_i) 
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




#############################################
## Asexuals; Proportion co-infected

Pmg_Pcoinf <- matrix(NA, nrow=N_tt, ncol=3)
colnames(Pmg_Pcoinf) <- c("med", "lwr", "upr")


for(i in 1:N_tt)
{
	Pmg_Pcoinf[i,1:3] <- as.numeric(binom.confint( Pmg_PLPQ_coinf[i,4]*Pmg_PLPQ_coinf[i,1], 
                                                     Pmg_PLPQ_coinf[i,4]*Pmg_PLPQ[i,1], 
                                                     method="wilson")[1,4:6])
}





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


##################################
## Bespoke circle plotting

semi_circle_plotter <- function(xx, yy, rr, left_colour, right_colour)
{
	N_plot <- 1000

	theta <- seq(-pi, pi, length=N_plot)
	
	x_plot <- xx + rr*cos(theta) 
	y_plot <- yy + rr*sin(theta) 

	polygon( x_plot[1:(N_plot/2)], y_plot[1:(N_plot/2)], 
               col=left_colour, border=left_colour)

	polygon( x_plot[(N_plot/2+1):N_plot], y_plot[(N_plot/2+1):N_plot], 
               col=right_colour, border=right_colour)

}




Pf_col <- "red"
Pv_col <- "orange"
Pm_col <- "forestgreen"
Po_col <- "dodgerblue"




main.size = 1.3
axis.size = 0.7
lab.size  = 1.25

point.size = 1
line.size  = 1
rr.size1    = 0.005
rr.size2    = 0.0075
arrow.edge = 0.015
arrow.size = 0.5
text.size = 0.5


tiff( file="Figure5_Pm_with_Pfcoinf.tif", width=10, height=16, units="cm", res=500)


lay.mat = rbind( c( 1 ),
                 c( 2 ),
                 c( 3 ), 
                 c( 4 ), 
                 c( 5 ) )

layout(lay.mat, heights=c(10,10,1,10,1))
layout.show(5)

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
type='l', lwd=line.size, col=Pf_col)



points(x=tt+1, y=Pm_PLPQ_coinf[,1], 
pch=19, cex=point.size, col=Pf_col)

for(i in 1:N_tt)
{
	semi_circle_plotter( xx=tt[i]+1, yy=Pm_PLPQ_coinf[i,1], rr=rr.size1, 
      	               left_colour=Pf_col, right_colour=Pm_col)
}

for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pm_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pm_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	
}


axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=1.25*axis.size )


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
type='l', lwd=line.size, col=Pf_col)

points(x=tt+1, y=Pmg_PLPQ_coinf[,1], 
pch=19, cex=point.size, col=Pf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pmg_PLPQ_coinf[i,2], 
     		 x1=tt[i]+1, y1=Pmg_PLPQ_coinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	
}





axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=1.25*axis.size )



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
       legend = c("P. malariae (all infections)", "Pm & Pf co-infection"),
       #fill = c(Pm_col, Pf_col), 
       #border = c(Pm_col, Pf_col), 
       col = c(Pm_col, Pf_col), 
       pch=c(19,19), ncol=2, cex=1, bty="n" )



########################################
##                                    ##
##  PANEL 3: Gametocyte prevalence    ##
##          	                    ## 
########################################

par(mar = c(3,3,2,0.5))

par(mgp = c(1.95, 0.6, 0))

line_seq_x <- c(0, 0.25, 0.5, 0.75, 1)

line_seq_y <- c(-50, 0, 50, 100, 150, 200, 250)



plot(x=100, y=100, 
xlim=c(-50,250), ylim=c(0,1), 
xlab="time (days)", ylab="proportion co-infected",
main="(C) P. malariae asexual co-infection",
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

points(x=tt-1, y=Pm_Pcoinf[,1], 
type='l', lwd=line.size, col=Pf_col)

points(x=tt-1, y=Pm_Pcoinf[,1], 
pch=19, cex=point.size, col=Pf_col)


for(i in 1:N_tt)
{
	arrows(x0=tt[i]-1, y0=Pm_Pcoinf[i,2], 
     		 x1=tt[i]-1, y1=Pm_Pcoinf[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	
}



points(x=tt+1, y=Pf_PLPQ[,1], 
type='l', lwd=line.size, col="black")

points(x=tt+1, y=Pf_PLPQ[,1], 
pch=19, cex=point.size, col="black")


for(i in 1:N_tt)
{
	arrows(x0=tt[i]+1, y0=Pf_PLPQ[i,2], 
     		 x1=tt[i]+1, y1=Pf_PLPQ[i,3], 
    	 	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}








axis(1,  at = c(-50, 0, 50, 100, 150, 200, 250), 
     labels = c(-50, 0, 50, 100, 150, 200, 250),
cex.axis=1.25*axis.size )



axis(2, las=2, at = c(0, 0.25, 0.50, 0.75, 1), 
     labels = c("0%", "25%", "50%", "75%", "100%"),
cex.axis=axis.size )






########################################
##                                    ##
##  PANEL 3: Legend                   ##
##          	                    ## 
########################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Pm with Pf co-infection: observed", "Pm with Pf co-infection: expected"),
       #fill = c(Pm_col, Coinf_col), 
       #border = c(Pm_col, Coinf_col), 
       col = c(Pf_col, "black"), 
       pch=c(19,19), ncol=2, cex=1, bty="n" )


dev.off()








