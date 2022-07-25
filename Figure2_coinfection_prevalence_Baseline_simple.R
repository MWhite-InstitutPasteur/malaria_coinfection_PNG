library(binom)
library(LaplacesDemon)

###################################
###################################
##                               ## 
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
################################### 


ALB_par <- read.csv("Albinama_data_PAR.csv", sep=";")

ALB_par <- ALB_par[which(is.na(ALB_par$qmal_pos)==FALSE),]

study_id <- unique(ALB_par$studyid)
N_data <- length(study_id)

index_baseline <- rep(NA, length=N_data)

for(i in 1:N_data)
{
	index <- which( ALB_par$studyid == study_id[i] & ALB_par$study_day < 0 )

	if( length(index) == 1	)
	{
		index_baseline[i] <- index	
	}

	if( length(index) > 1	)
	{
		index_baseline[i] <- index[which.min(ALB_par[index,]$study_day)]
	}
}


ALB_par <- ALB_par[index_baseline,]





ALB_epi <- read.csv("Albimana_data_EPI.csv")

village_names <- unique(ALB_epi$villagegroup)
N_village <- length(village_names)


ALB_par_village <- rep(NA, N_data)

for(i in 1:N_data)
{
	ALB_par_village[i] <- as.character(ALB_epi$villagegroup[which(ALB_epi$studyid == ALB_par$studyid[i])][1])
}

ALB_par$village <- ALB_par_village



ALB_par_treat <- rep(NA, N_data)

for(i in 1:N_data)
{
	ALB_par_treat[i] <- as.character(ALB_epi$treat[which(ALB_epi$studyid == ALB_par$studyid[i])][1])
}

ALB_par$treat <- ALB_par_treat



#########################################
#########################################
##                                     ##
##  #     #  ####  ####   ##### ##     ##
##  ##   ## ##  ## ## ##  ##    ##     ##
##  ####### ##  ## ##  ## ####  ##     ##
##  ## # ## ##  ## ## ##  ##    ##     ##
##  ##   ##  ####  ####   ##### #####  ##
##                                     ## 
#########################################
#########################################

N_seq <- 10000

prev_seq <- exp(seq(from=log(1e-10), to=log(1), length=N_seq))


#########################################
#########################################
#########################################
##                                     ##
##  1. Asexual parasites               ##
##                                     ##
#########################################
#########################################
#########################################

#########################################
#########################################
##                                     ##
##  1.1. Baseline data                 ##
##                                     ##
#########################################
#########################################

N_data_base <- length(which(ALB_par$study_day<0))

N_f_base <- length( which(ALB_par$pf_pos==1 & ALB_par$study_day<0) )
N_v_base <- length( which(ALB_par$pv_pos==1 & ALB_par$study_day<0) )
N_m_base <- length( which(ALB_par$pm_pos==1 & ALB_par$study_day<0) )
N_o_base <- length( which(ALB_par$po_pos==1 & ALB_par$study_day<0) )

N_fv_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$study_day<0) )
N_fm_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$study_day<0) )
N_fo_base <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )
N_vm_base <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$study_day<0) )
N_vo_base <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )
N_mo_base <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )

N_fvm_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$study_day<0) )
N_fvo_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )
N_fmo_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )
N_vmo_base <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )

N_fvmo_base <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$study_day<0) )




#########################################
##                                     ##
##  1.11. Baseline data - double       ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double_base     <- matrix(NA, nrow=4, ncol=4)
ASX_double_obs_base <- matrix(NA, nrow=4, ncol=4)
ASX_double_exp_base <- matrix(NA, nrow=4, ncol=4)

ASX_double_obs_low_base  <- matrix(NA, nrow=4, ncol=4)
ASX_double_obs_high_base <- matrix(NA, nrow=4, ncol=4)


ASX_double_base <- rbind( c(N_f_base,  N_fv_base, N_fm_base, N_fo_base),
     	 	              c(N_fv_base, N_v_base,  N_vm_base, N_vo_base),
           		        c(N_fm_base, N_vm_base, N_m_base,  N_mo_base),
           	              c(N_fo_base, N_vo_base, N_mo_base, N_o_base ) )

colnames(ASX_double_base) <- c("Pf", "Pv", "Po", "Pm")
rownames(ASX_double_base) <- c("Pf", "Pv", "Po", "Pm")

ASX_double_obs_base <- ASX_double_base/N_data_base

ASX_double_exp_base <- (diag(ASX_double_base)/N_data_base) %o% (diag(ASX_double_base)/N_data_base)
diag(ASX_double_exp_base) <- diag(ASX_double_obs_base)


ASX_double_obs_low_base  <- ASX_double_base
ASX_double_obs_high_base <- ASX_double_base

for(i in 1:4)
{
	for(j in 1:4)
	{
		ASX_double_obs_low_base[i,j]  <- binom.confint( ASX_double_base[i,j], N_data_base, method="wilson")[1,5]
		ASX_double_obs_high_base[i,j] <- binom.confint( ASX_double_base[i,j], N_data_base, method="wilson")[1,6]
	}
}




#########################################
## Regression model

x_exp <- c( ASX_double_exp_base[1,2:4], ASX_double_exp_base[2,3:4], ASX_double_exp_base[3,4] )
y_obs <- c( ASX_double_obs_base[1,2:4], ASX_double_obs_base[2,3:4], ASX_double_obs_base[3,4] )

ASX_double_base_mod <- lm( y_obs ~ 0 + x_exp )




ASX_double_base_mod_fit <- predict.lm( ASX_double_base_mod, data.frame(x_exp=prev_seq), interval="prediction", se.fit=TRUE)

ASX_double_base_mod_fit$fit[which(ASX_double_base_mod_fit$fit[,2] < 0),2] <- 1e-10



c( summary(ASX_double_base_mod)$coef[1],
   summary(ASX_double_base_mod)$coef[1] - 1.96*summary(ASX_double_base_mod)$coef[2], 
   summary(ASX_double_base_mod)$coef[1] + 1.96*summary(ASX_double_base_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_base_mod)$coef[1], 
         sigma = summary(ASX_double_base_mod)$coef[2], 
         nu    = summary(ASX_double_base_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_base_mod)$coef[1], 
         sigma = summary(ASX_double_base_mod)$coef[2], 
         nu    = summary(ASX_double_base_mod)$df[2])




#########################################
#########################################
##                                     ##
##  1.2. Baseline data                 ##
##       Village stratified            ##
##                                     ##
#########################################
#########################################

N_data_vil_base <- rep(NA, N_village)

N_f_vil_base <- rep(NA, N_village)
N_v_vil_base <- rep(NA, N_village)
N_m_vil_base <- rep(NA, N_village)
N_o_vil_base <- rep(NA, N_village)

N_fv_vil_base <- rep(NA, N_village)
N_fm_vil_base <- rep(NA, N_village)
N_fo_vil_base <- rep(NA, N_village)
N_vm_vil_base <- rep(NA, N_village)
N_vo_vil_base <- rep(NA, N_village)
N_mo_vil_base <- rep(NA, N_village)

N_fvm_vil_base <- rep(NA, N_village)
N_fvo_vil_base <- rep(NA, N_village)
N_fmo_vil_base <- rep(NA, N_village)
N_vmo_vil_base <- rep(NA, N_village)

N_fvmo_vil_base <- rep(NA, N_village)

for(v in 1:N_village)
{
	N_data_vil_base[v] <- length(which(ALB_par$village == village_names[v] & ALB_par$study_day<0))

	N_f_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_v_vil_base[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_m_vil_base[v] <- length( which(ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_o_vil_base[v] <- length( which(ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fv_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fm_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fo_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vm_vil_base[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vo_vil_base[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_mo_vil_base[v] <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fvm_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fvo_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fmo_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vmo_vil_base[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fvmo_vil_base[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
}



#########################################
##                                     ##
##  1.2.1. Baseline data - double      ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double_vil_base <- array(NA, dim=c(N_village, 4, 4))
ASX_double_obs_vil_base <- array(NA, dim=c(N_village, 4, 4))
ASX_double_exp_vil_base <- array(NA, dim=c(N_village, 4, 4))

ASX_double_obs_low_vil_base <- array(NA, dim=c(N_village, 4, 4))
ASX_double_obs_high_vil_base <- array(NA, dim=c(N_village, 4, 4))


for(v in 1:N_village)
{
	ASX_double_vil_base[v,,] <- rbind( c(N_f_vil_base[v],  N_fv_vil_base[v], N_fm_vil_base[v], N_fo_vil_base[v]),
      	 	                       c(N_fv_vil_base[v], N_v_vil_base[v],  N_vm_vil_base[v], N_vo_vil_base[v]),
            		                 c(N_fm_vil_base[v], N_vm_vil_base[v], N_m_vil_base[v],  N_mo_vil_base[v]),
              	                       c(N_fo_vil_base[v], N_vo_vil_base[v], N_mo_vil_base[v], N_o_vil_base[v] ) )

	colnames(ASX_double_vil_base[v,,]) <- c("Pf", "Pv", "Po", "Pm")
	rownames(ASX_double_vil_base[v,,]) <- c("Pf", "Pv", "Po", "Pm")

	ASX_double_obs_vil_base[v,,] <- ASX_double_vil_base[v,,]/N_data_vil_base[v]

	ASX_double_exp_vil_base[v,,] <- (diag(ASX_double_vil_base[v,,])/N_data_vil_base[v]) %o% (diag(ASX_double_vil_base[v,,])/N_data_vil_base[v])
	diag(ASX_double_exp_vil_base[v,,]) <- diag(ASX_double_obs_vil_base[v,,])


	ASX_double_obs_low_vil_base[v,,]  <- ASX_double_vil_base[v,,]
	ASX_double_obs_high_vil_base[v,,] <- ASX_double_vil_base[v,,]

	for(i in 1:4)
	{
		for(j in 1:4)
		{
			ASX_double_obs_low_vil_base[v,i,j]  <- binom.confint( ASX_double_vil_base[v,i,j], N_data_vil_base[v], method="wilson")[1,5]
			ASX_double_obs_high_vil_base[v,i,j] <- binom.confint( ASX_double_vil_base[v,i,j], N_data_vil_base[v], method="wilson")[1,6]
		}
	}

}


#########################################
## Regression model

x_exp <- c( ASX_double_exp_vil_base[,1,2:4], ASX_double_exp_vil_base[,2,3:4], ASX_double_exp_vil_base[,3,4] )
y_obs <- c( ASX_double_obs_vil_base[,1,2:4], ASX_double_obs_vil_base[,2,3:4], ASX_double_obs_vil_base[,3,4] )

ASX_double_vil_base_mod <- lm( y_obs ~ 0 + x_exp )




ASX_double_vil_base_mod_fit <- predict.lm( ASX_double_vil_base_mod, data.frame(x_exp=prev_seq), interval="prediction", se.fit=TRUE)

ASX_double_vil_base_mod_fit$fit[which(ASX_double_vil_base_mod_fit$fit[,2] < 0),2] <- 1e-10



c( summary(ASX_double_vil_base_mod)$coef[1],
   summary(ASX_double_vil_base_mod)$coef[1] - 1.96*summary(ASX_double_vil_base_mod)$coef[2], 
   summary(ASX_double_vil_base_mod)$coef[1] + 1.96*summary(ASX_double_vil_base_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])





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



tri_circle_plotter <- function(xx, yy, rr, colour_1, colour_2, colour_3)
{
	N_plot <- 1000

	theta <- seq(-pi, pi, length=N_plot)
	
	x_plot <- xx + rr*cos(theta) 
	y_plot <- yy + rr*sin(theta) 

	polygon( c(xx, x_plot[1:(N_plot/3)]), 
               c(yy, y_plot[1:(N_plot/3)]), 
               col=colour_1, border=colour_1)

	polygon( c(xx, x_plot[(N_plot/3+1):(2*N_plot/3)]), 
               c(yy, y_plot[(N_plot/3+1):(2*N_plot/3)]), 
               col=colour_2, border=colour_2)

	polygon( c(xx, x_plot[(2*N_plot/3+1):N_plot]), 
               c(yy, y_plot[(2*N_plot/3+1):N_plot]), 
               col=colour_3, border=colour_3)
}





##################################
##################################
##                              ## 
##  Co-infection plots          ##
##                              ##
##################################
##################################


Pf_col <- "red"
Pv_col <- "orange"
Pm_col <- "forestgreen"
Po_col <- "dodgerblue"



main.size = 1.2
axis.size = 0.65
lab.size  = 1.2

point.size = 1.5
line.size  = 1
rr.size1    = 0.0075
rr.size2    = 0.01
rr.size3    = 0.0125
arrow.edge = 0.015
arrow.size = 0.75



tiff( file="Figure2_coinfection_prevalence_Baseline_simple.tif", width=14, height=18, units="cm", res=500)


lay.mat = rbind( c( 1, 2 ),
                 c( 3, 3 ) )

layout(lay.mat, heights=c(1,2))
layout.show(3)



########################################
##                                    ##
##  PANEL 1: Asexual co-infection     ##
##          	                    ## 
########################################

par(mar = c(3,3,2,0.5))

par(mgp = c(1.7, 0.6, 0))


X_Pf <- N_f_base/N_data_base
X_Pv <- N_v_base/N_data_base

X_Pfv <- N_fv_base/N_data_base


plot(x=100, y=100, 
xlim=c(0,1), ylim=c(0,1), 
xlab="P. vivax prevalence", ylab="P. falciparum prevalence",
main="(A) Expected co-infection",
xaxt="n", yaxt="n", xaxs='i', yaxs='i', bty='n',
cex.main=main.size, cex.lab=0.8*lab.size)

polygon( x = c(0, X_Pv, X_Pv, 0),
         y = c(0, 0, 1, 1), col=Pv_col, border=NA ) 

polygon( x = c(0, 1, 1, 0),
         y = c(0, 0, X_Pf, X_Pf), 
	   col=Pf_col, border=NA ) 


polygon( x = c(0, X_Pv, X_Pv, 0),
         y = c(0, 0, X_Pf, X_Pf), 
         density = 10, lwd = 5, angle = 45,
         col=Pv_col, border=NA ) 

axis(1,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%", "40%", "60%", "80%", "100%" ),
cex.axis=axis.size )

axis(2,  las=2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%", "40%", "60%", "80%", "100%" ),
cex.axis=axis.size )



########################################
##                                    ##
##  PANEL 2: Asexual co-infection     ##
##          	                    ## 
########################################
#
#Z_v <- sqrt( X_Pfv*X_Pv/X_Pf )
#Z_f <- sqrt( X_Pfv*X_Pf/X_Pv )
#
#Y_v = 2*(X_Pv - X_Pfv)/( 1 - Z_f ) - Z_v
#Y_f = 2*(X_Pf - X_Pfv)/( 1 - Z_v ) - Z_f
#
#
#
#plot(x=100, y=100, 
#xlim=c(0,1), ylim=c(0,1), 
#xlab="expected prevalence", ylab="observed prevalence",
#main="(A) Malaria co-infection",
#xaxt="n", yaxt="n", xaxs='i', yaxs='i', bty='n',
#cex.main=main.size, cex.lab=lab.size)
#
#
#polygon( x = c(0, Z_v, Z_v, 0),
#         y = c(0, 0, Z_f, Z_f), 
#         col=Pf_col, border=NA ) 
#
#
#polygon( x = c(0, Z_v, Z_v, 0),
#         y = c(0, 0, Z_f, Z_f),
#         density = 10, lwd = 5, angle = 45,
#         col=Pv_col, border=NA ) 
#
#
#polygon( x = c( Z_v, 1, 1, Z_v),
#         y = c(0, 0, Y_f, Z_f), 
#         col=Pf_col, border=NA ) 
#
#polygon( x = c( 0, Z_v, Y_v, 0),
#         y = c(Z_f, Z_f, 1, 1), 
#         col=Pv_col, border=NA ) 


########################################
##                                    ##
##  PANEL 2: Asexual co-infection     ##
##          	                    ## 
########################################

Z_v <- sqrt( X_Pfv*X_Pv/X_Pf )
Z_f <- sqrt( X_Pfv*X_Pf/X_Pv )

Y_v = (X_Pv - X_Pfv)/( 1 - Z_f )
Y_f = (X_Pf - X_Pfv)/( 1 - Z_v ) 


plot(x=100, y=100, 
xlim=c(0,1), ylim=c(0,1), 
xlab="P. vivax prevalence", ylab="P. falciparum prevalence",
main="(B) Observed co-infection",
xaxt="n", yaxt="n", xaxs='i', yaxs='i', bty='n',
cex.main=main.size, cex.lab=0.8*lab.size)


polygon( x = c(0, Z_v, Z_v, 0),
         y = c(0, 0, Z_f, Z_f), 
         col=Pf_col, border=NA ) 


polygon( x = c(0, Z_v, Z_v, 0),
         y = c(0, 0, Z_f, Z_f),
         density = 10, lwd = 5, angle = 45,
         col=Pv_col, border=NA ) 


polygon( x = c( Z_v, 1, 1, Z_v),
         y = c(0, 0, Y_f, Y_f), 
         col=Pf_col, border=NA ) 

polygon( x = c( 0, Y_v, Y_v, 0),
         y = c(Z_f, Z_f, 1, 1), 
         col=Pv_col, border=NA )

axis(1,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%", "40%", "60%", "80%", "100%" ),
cex.axis=axis.size )

axis(2,  las=2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%", "40%", "60%", "80%", "100%" ),
cex.axis=axis.size )

########################################
##                                    ##
##  PANEL 3: Asexual co-infection     ##
##          	                    ## 
########################################



par(mar = c(3,3,2,0.5))

par(mgp = c(1.6, 0.5, 0))


ff <- 0.25  ## scaling factor for plotting



line_seq_x <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)^ff

line_seq_y <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)^ff



plot(x=100, y=100, 
xlim=c(0,1^ff), ylim=c(0,1^ff), 
xlab="expected prevalence", ylab="observed prevalence",
main="(C) Multi-species co-infection",
xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.lab=lab.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(-1e10,1), type='l', lwd=0.5, col="grey", lty="dashed")
}


points(x=c(0,1), y=c(0,1), 
type='l', lwd=1.5, lty="dotted")




points(x=prev_seq^ff, y=ASX_double_base_mod_fit$fit[,1]^ff, type='l', lwd=1.5, lty="solid")

polygon(x=c(prev_seq^ff, rev(prev_seq^ff)), 
	  y=c( ASX_double_base_mod_fit$fit[,2]^ff, rev(ASX_double_base_mod_fit$fit[,3]^ff) ),
	  col=rgb(100/256,100/256,100/256,0.1), border=NA)


################################
## Pf single infection

semi_circle_plotter( xx=ASX_double_exp_base[1,1]^ff, yy=ASX_double_obs_base[1,1]^ff, rr=rr.size1, 
                     left_colour=Pf_col, right_colour=Pf_col)

arrows(x0=ASX_double_exp_base[1,1]^ff, y0=ASX_double_obs_low_base[1,1]^ff, 
     	 x1=ASX_double_exp_base[1,1]^ff, y1=ASX_double_obs_high_base[1,1]^ff, 
    	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	


################################
## Pv single infection

semi_circle_plotter( xx=ASX_double_exp_base[2,2]^ff, yy=ASX_double_obs_base[2,2]^ff, rr=rr.size1, 
     	               left_colour=Pv_col, right_colour=Pv_col)

arrows(x0=ASX_double_exp_base[2,2]^ff, y0=ASX_double_obs_low_base[2,2]^ff, 
     	 x1=ASX_double_exp_base[2,2]^ff, y1=ASX_double_obs_high_base[2,2]^ff, 
    	 length=arrow.edge, angle=90, code=3, col=Pv_col, lwd=arrow.size)	


################################
## Pm single infection

semi_circle_plotter( xx=ASX_double_exp_base[3,3]^ff, yy=ASX_double_obs_base[3,3]^ff, rr=rr.size1, 
                     left_colour=Pm_col, right_colour=Pm_col)

arrows(x0=ASX_double_exp_base[3,3]^ff, y0=ASX_double_obs_low_base[3,3]^ff, 
     	 x1=ASX_double_exp_base[3,3]^ff, y1=ASX_double_obs_high_base[3,3]^ff, 
    	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	


################################
## Po single infection

semi_circle_plotter( xx=ASX_double_exp_base[4,4]^ff, yy=ASX_double_obs_base[4,4]^ff, rr=rr.size1, 
                     left_colour=Po_col, right_colour=Po_col)

arrows(x0=ASX_double_exp_base[4,4]^ff, y0=ASX_double_obs_low_base[4,4]^ff, 
     	 x1=ASX_double_exp_base[4,4]^ff, y1=ASX_double_obs_high_base[4,4]^ff, 
    	 length=arrow.edge, angle=90, code=3, col=Po_col, lwd=arrow.size)	



################################
## Pf + Pv co-infection


arrows(x0=ASX_double_exp_base[1,2]^ff, y0=ASX_double_obs_low_base[1,2]^ff, 
     	 x1=ASX_double_exp_base[1,2]^ff, y1=ASX_double_obs_high_base[1,2]^ff, 
    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[1,2]^ff, yy=ASX_double_obs_base[1,2]^ff, rr=rr.size2, 
                     left_colour=Pf_col, right_colour=Pv_col)


################################
## Pf + Pm co-infection

arrows(x0=ASX_double_exp_base[1,3]^ff, y0=ASX_double_obs_low_base[1,3]^ff, 
     	 x1=ASX_double_exp_base[1,3]^ff, y1=ASX_double_obs_high_base[1,3]^ff, 
    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[1,3]^ff, yy=ASX_double_obs_base[1,3]^ff, rr=rr.size2, 
                     left_colour=Pf_col, right_colour=Pm_col)



################################
## Pf + Po co-infection

arrows(x0=ASX_double_exp_base[1,4]^ff, y0=ASX_double_obs_low_base[1,4]^ff, 
     	 x1=ASX_double_exp_base[1,4]^ff, y1=ASX_double_obs_high_base[1,4]^ff, 
    	 length=0.03, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[1,4]^ff, yy=ASX_double_obs_base[1,4]^ff, rr=rr.size2, 
                     left_colour=Pf_col, right_colour=Po_col)


################################
## Pv + Pm co-infection

arrows(x0=ASX_double_exp_base[2,3]^ff, y0=ASX_double_obs_low_base[2,3]^ff, 
     	 x1=ASX_double_exp_base[2,3]^ff, y1=ASX_double_obs_high_base[2,3]^ff, 
    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[2,3]^ff, yy=ASX_double_obs_base[2,3]^ff, rr=rr.size2, 
                     left_colour=Pv_col, right_colour=Pm_col)


################################
## Pv + Po co-infection


arrows(x0=ASX_double_exp_base[2,4]^ff, y0=ASX_double_obs_low_base[2,4]^ff, 
     	 x1=ASX_double_exp_base[2,4]^ff, y1=ASX_double_obs_high_base[2,4]^ff, 
    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[2,4]^ff, yy=ASX_double_obs_base[2,4]^ff, rr=rr.size2, 
                     left_colour=Pv_col, right_colour=Po_col)


################################
## Pm + Po co-infection

arrows(x0=ASX_double_exp_base[3,4]^ff, y0=ASX_double_obs_low_base[3,4]^ff, 
	 x1=ASX_double_exp_base[3,4]^ff, y1=ASX_double_obs_high_base[3,4]^ff, 
    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

semi_circle_plotter( xx=ASX_double_exp_base[3,4]^ff, yy=ASX_double_obs_base[3,4]^ff, rr=rr.size2, 
                     left_colour=Pm_col, right_colour=Po_col)



axis(1,  at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%" ),
cex.axis=axis.size )

axis(1,  at = c(0.001, 0.01)^ff, 
     labels = c("0.1%", "1%" ),
cex.axis=axis.size )

axis(2, las=2, at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%", "100%" ),
cex.axis=axis.size )







########################################
##                                    ##
## Legend 


legend(x='bottomright', 
       legend = c("P. falciparum", "P. vivax", "P. malariae", "P. ovale"),
       fill = c(Pf_col, Pv_col, Pm_col, Po_col), 
       border = c(Pf_col, Pv_col, Pm_col, Po_col), 
       ncol=1, cex=1.5, bty="n" )



dev.off()


