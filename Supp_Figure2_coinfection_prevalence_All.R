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

N_data <- nrow(ALB_par)


ALB_epi <- read.csv("Albimana_data_EPI.csv")

village_names <- unique(ALB_epi$villagegroup_corrected)
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
##  1.2. All data                      ##
##       Village stratified            ##
##                                     ##
#########################################
#########################################

N_data_vil <- rep(NA, N_village)

N_f_vil <- rep(NA, N_village)
N_v_vil <- rep(NA, N_village)
N_m_vil <- rep(NA, N_village)
N_o_vil <- rep(NA, N_village)

N_fv_vil <- rep(NA, N_village)
N_fm_vil <- rep(NA, N_village)
N_fo_vil <- rep(NA, N_village)
N_vm_vil <- rep(NA, N_village)
N_vo_vil <- rep(NA, N_village)
N_mo_vil <- rep(NA, N_village)

N_fvm_vil <- rep(NA, N_village)
N_fvo_vil <- rep(NA, N_village)
N_fmo_vil <- rep(NA, N_village)
N_vmo_vil <- rep(NA, N_village)

N_fvmo_vil <- rep(NA, N_village)

for(v in 1:N_village)
{
	N_data_vil[v] <- length(which(ALB_par$village == village_names[v]))

	N_f_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$village==village_names[v]) )
	N_v_vil[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$village==village_names[v]) )
	N_m_vil[v] <- length( which(ALB_par$pm_pos==1 & ALB_par$village==village_names[v]) )
	N_o_vil[v] <- length( which(ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )

	N_fv_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$village==village_names[v]) )
	N_fm_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v]) )
	N_fo_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )
	N_vm_vil[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v]) )
	N_vo_vil[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )
	N_mo_vil[v] <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )

	N_fvm_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v]) )
	N_fvo_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )
	N_fmo_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )
	N_vmo_vil[v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )

	N_fvmo_vil[v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v]) )
}



#########################################
##                                     ##
##  1.2.1. All data - double           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double_vil <- array(NA, dim=c(N_village, 4, 4))
ASX_double_obs_vil <- array(NA, dim=c(N_village, 4, 4))
ASX_double_exp_vil <- array(NA, dim=c(N_village, 4, 4))

ASX_double_obs_low_vil <- array(NA, dim=c(N_village, 4, 4))
ASX_double_obs_high_vil <- array(NA, dim=c(N_village, 4, 4))


for(v in 1:N_village)
{
	ASX_double_vil[v,,] <- rbind( c(N_f_vil[v],  N_fv_vil[v], N_fm_vil[v], N_fo_vil[v]),
      	 	                  c(N_fv_vil[v], N_v_vil[v],  N_vm_vil[v], N_vo_vil[v]),
            		            c(N_fm_vil[v], N_vm_vil[v], N_m_vil[v],  N_mo_vil[v]),
              	                  c(N_fo_vil[v], N_vo_vil[v], N_mo_vil[v], N_o_vil[v] ) )

	colnames(ASX_double_vil[v,,]) <- c("Pf", "Pv", "Po", "Pm")
	rownames(ASX_double_vil[v,,]) <- c("Pf", "Pv", "Po", "Pm")

	ASX_double_obs_vil[v,,] <- ASX_double_vil[v,,]/N_data_vil[v]

	ASX_double_exp_vil[v,,] <- (diag(ASX_double_vil[v,,])/N_data_vil[v]) %o% (diag(ASX_double_vil[v,,])/N_data_vil[v])
	diag(ASX_double_exp_vil[v,,]) <- diag(ASX_double_obs_vil[v,,])


	ASX_double_obs_low_vil[v,,]  <- ASX_double_vil[v,,]
	ASX_double_obs_high_vil[v,,] <- ASX_double_vil[v,,]

	for(i in 1:4)
	{
		for(j in 1:4)
		{
			ASX_double_obs_low_vil[v,i,j]  <- binom.confint( ASX_double_vil[v,i,j], N_data_vil[v], method="wilson")[1,5]
			ASX_double_obs_high_vil[v,i,j] <- binom.confint( ASX_double_vil[v,i,j], N_data_vil[v], method="wilson")[1,6]
		}
	}

}


#########################################
## Regression model

x_exp <- c( ASX_double_exp_vil[,1,2:4], ASX_double_exp_vil[,2,3:4], ASX_double_exp_vil[,3,4] )
y_obs <- c( ASX_double_obs_vil[,1,2:4], ASX_double_obs_vil[,2,3:4], ASX_double_obs_vil[,3,4] )

ASX_double_vil_mod <- lm( y_obs ~ 0 + x_exp )




ASX_double_vil_mod_fit <- predict.lm( ASX_double_vil_mod, data.frame(x_exp=prev_seq), interval="prediction", se.fit=TRUE)

ASX_double_vil_mod_fit$fit[which(ASX_double_vil_mod_fit$fit[,2] < 0),2] <- 1e-10



c( summary(ASX_double_vil_mod)$coef[1],
   summary(ASX_double_vil_mod)$coef[1] - 1.96*summary(ASX_double_vil_mod)$coef[2], 
   summary(ASX_double_vil_mod)$coef[1] + 1.96*summary(ASX_double_vil_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_mod)$coef[1], 
         sigma = summary(ASX_double_vil_mod)$coef[2], 
         nu    = summary(ASX_double_vil_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_mod)$coef[1], 
         sigma = summary(ASX_double_vil_mod)$coef[2], 
         nu    = summary(ASX_double_vil_mod)$df[2])



#########################################
##                                     ##
##  1.2.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

ASX_triple_vil <- c(N_vmo_vil, N_fmo_vil, N_fvo_vil, N_fvm_vil)
names(ASX_triple_vil) <- c( rep("P_vmo", N_village), rep("P_fmo", N_village), rep("P_fvo", N_village), rep("P_fvm", N_village) )

ASX_triple_vil_obs      <- ASX_triple_vil
ASX_triple_vil_obs_low  <- ASX_triple_vil
ASX_triple_vil_obs_high <- ASX_triple_vil

for(v in 1:N_village)
{
	ASX_triple_vil_obs[0*N_village + v] <- ASX_triple_vil[0*N_village + v]/N_data_vil[v]
	ASX_triple_vil_obs[1*N_village + v] <- ASX_triple_vil[1*N_village + v]/N_data_vil[v]
	ASX_triple_vil_obs[2*N_village + v] <- ASX_triple_vil[2*N_village + v]/N_data_vil[v]
	ASX_triple_vil_obs[3*N_village + v] <- ASX_triple_vil[3*N_village + v]/N_data_vil[v]
}


for(v in 1:N_village)
{
	ASX_triple_vil_obs_low[0*N_village + v] <- binom.confint( ASX_triple_vil[0*N_village + v], N_data_vil[v], method="wilson")[1,5]
	ASX_triple_vil_obs_low[1*N_village + v] <- binom.confint( ASX_triple_vil[1*N_village + v], N_data_vil[v], method="wilson")[1,5]
	ASX_triple_vil_obs_low[2*N_village + v] <- binom.confint( ASX_triple_vil[2*N_village + v], N_data_vil[v], method="wilson")[1,5]
	ASX_triple_vil_obs_low[3*N_village + v] <- binom.confint( ASX_triple_vil[3*N_village + v], N_data_vil[v], method="wilson")[1,5]

	ASX_triple_vil_obs_high[0*N_village + v] <- binom.confint( ASX_triple_vil[0*N_village + v], N_data_vil[v], method="wilson")[1,6]
	ASX_triple_vil_obs_high[1*N_village + v] <- binom.confint( ASX_triple_vil[1*N_village + v], N_data_vil[v], method="wilson")[1,6]
	ASX_triple_vil_obs_high[2*N_village + v] <- binom.confint( ASX_triple_vil[2*N_village + v], N_data_vil[v], method="wilson")[1,6]
	ASX_triple_vil_obs_high[3*N_village + v] <- binom.confint( ASX_triple_vil[3*N_village + v], N_data_vil[v], method="wilson")[1,6]
}





ASX_triple_vil_exp <- ASX_triple_vil_obs

for(v in 1:N_village)
{
	ASX_triple_vil_exp[0*N_village + v] <- N_v_vil[v]*N_m_vil[v]*N_o_vil[v]/( N_data_vil[v]^3 )
	ASX_triple_vil_exp[1*N_village + v] <- N_f_vil[v]*N_m_vil[v]*N_o_vil[v]/( N_data_vil[v]^3 )
	ASX_triple_vil_exp[2*N_village + v] <- N_f_vil[v]*N_v_vil[v]*N_o_vil[v]/( N_data_vil[v]^3 )
	ASX_triple_vil_exp[3*N_village + v] <- N_f_vil[v]*N_v_vil[v]*N_m_vil[v]/( N_data_vil[v]^3 )
}


#########################################
## Regression model

ASX_triple_vil_mod <- lm( ASX_triple_vil_obs ~ 0 + ASX_triple_vil_exp )



ASX_triple_vil_mod_fit <- predict.lm( ASX_triple_vil_mod, data.frame(ASX_triple_vil_exp=prev_seq), interval="prediction", se.fit=TRUE)

ASX_triple_vil_mod_fit$fit[which(ASX_triple_vil_mod_fit$fit[,2] < 0),2] <- 1e-10


c( summary(ASX_double_vil_mod)$coef[1],
   summary(ASX_double_vil_mod)$coef[1] - 1.96*summary(ASX_double_vil_mod)$coef[2], 
   summary(ASX_double_vil_mod)$coef[1] + 1.96*summary(ASX_double_vil_mod)$coef[2] )^2

c( summary(ASX_triple_vil_mod)$coef[1],
   summary(ASX_triple_vil_mod)$coef[1] - 1.96*summary(ASX_triple_vil_mod)$coef[2], 
   summary(ASX_triple_vil_mod)$coef[1] + 1.96*summary(ASX_triple_vil_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_mod)$df[2])



#########################################
#########################################
#########################################
##                                     ##
##  2. Gametocytes                     ##
##                                     ##
#########################################
#########################################
#########################################

#########################################
#########################################
##                                     ##
##  2.2. All data                      ##
##       Village stratified            ##
##                                     ##
#########################################
#########################################

N_data_vil <- rep(NA, N_village)

N_fg_vil <- rep(NA, N_village)
N_vg_vil <- rep(NA, N_village)
N_mg_vil <- rep(NA, N_village)
N_og_vil <- rep(NA, N_village)

N_fvg_vil <- rep(NA, N_village)
N_fmg_vil <- rep(NA, N_village)
N_fog_vil <- rep(NA, N_village)
N_vmg_vil <- rep(NA, N_village)
N_vog_vil <- rep(NA, N_village)
N_mog_vil <- rep(NA, N_village)

N_fvmg_vil <- rep(NA, N_village)
N_fvog_vil <- rep(NA, N_village)
N_fmog_vil <- rep(NA, N_village)
N_vmog_vil <- rep(NA, N_village)

N_fvmog_vil <- rep(NA, N_village)

for(v in 1:N_village)
{
	N_data_vil[v] <- length(which(ALB_par$village == village_names[v]))

	N_fg_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$village==village_names[v]) )
	N_vg_vil[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$village==village_names[v]) )
	N_mg_vil[v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$village==village_names[v]) )
	N_og_vil[v] <- length( which(ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )

	N_fvg_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$village==village_names[v]) )
	N_fmg_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v]) )
	N_fog_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )
	N_vmg_vil[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v]) )
	N_vog_vil[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )
	N_mog_vil[v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )

	N_fvmg_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v]) )
	N_fvog_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )
	N_fmog_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )
	N_vmog_vil[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )

	N_fvmog_vil[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v]) )
}



#########################################
##                                     ##
##  2.2.1. All data - double           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

GAM_double_vil <- array(NA, dim=c(N_village, 4, 4))
GAM_double_obs_vil <- array(NA, dim=c(N_village, 4, 4))
GAM_double_exp_vil <- array(NA, dim=c(N_village, 4, 4))

GAM_double_obs_low_vil <- array(NA, dim=c(N_village, 4, 4))
GAM_double_obs_high_vil <- array(NA, dim=c(N_village, 4, 4))


for(v in 1:N_village)
{
	GAM_double_vil[v,,] <- rbind( c(N_fg_vil[v],  N_fvg_vil[v], N_fmg_vil[v], N_fog_vil[v]),
      	 	                  c(N_fvg_vil[v], N_vg_vil[v],  N_vmg_vil[v], N_vog_vil[v]),
            		            c(N_fmg_vil[v], N_vmg_vil[v], N_mg_vil[v],  N_mog_vil[v]),
              	                  c(N_fog_vil[v], N_vog_vil[v], N_mog_vil[v], N_og_vil[v] ) )

	colnames(GAM_double_vil[v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")
	rownames(GAM_double_vil[v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")

	GAM_double_obs_vil[v,,] <- GAM_double_vil[v,,]/N_data_vil[v]

	GAM_double_exp_vil[v,,] <- (diag(GAM_double_vil[v,,])/N_data_vil[v]) %o% (diag(GAM_double_vil[v,,])/N_data_vil[v])
	diag(GAM_double_exp_vil[v,,]) <- diag(GAM_double_obs_vil[v,,])


	GAM_double_obs_low_vil[v,,]  <- GAM_double_vil[v,,]
	GAM_double_obs_high_vil[v,,] <- GAM_double_vil[v,,]

	for(i in 1:4)
	{
		for(j in 1:4)
		{
			GAM_double_obs_low_vil[v,i,j]  <- binom.confint( GAM_double_vil[v,i,j], N_data_vil[v], method="wilson")[1,5]
			GAM_double_obs_high_vil[v,i,j] <- binom.confint( GAM_double_vil[v,i,j], N_data_vil[v], method="wilson")[1,6]
		}
	}

}


#########################################
## Regression model

x_exp <- c( GAM_double_exp_vil[,1,2:4], GAM_double_exp_vil[,2,3:4], GAM_double_exp_vil[,3,4] )
y_obs <- c( GAM_double_obs_vil[,1,2:4], GAM_double_obs_vil[,2,3:4], GAM_double_obs_vil[,3,4] )

GAM_double_vil_mod <- lm( y_obs ~ 0 + x_exp )




GAM_double_vil_mod_fit <- predict.lm( GAM_double_vil_mod, data.frame(x_exp=prev_seq), interval="prediction", se.fit=TRUE)

GAM_double_vil_mod_fit$fit[which(GAM_double_vil_mod_fit$fit[,2] < 0),2] <- 1e-10



c( summary(GAM_double_vil_mod)$coef[1],
   summary(GAM_double_vil_mod)$coef[1] - 1.96*summary(GAM_double_vil_mod)$coef[2], 
   summary(GAM_double_vil_mod)$coef[1] + 1.96*summary(GAM_double_vil_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_mod)$coef[1], 
         sigma = summary(GAM_double_vil_mod)$coef[2], 
         nu    = summary(GAM_double_vil_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_mod)$coef[1], 
         sigma = summary(GAM_double_vil_mod)$coef[2], 
         nu    = summary(GAM_double_vil_mod)$df[2])



#########################################
##                                     ##
##  2.2.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

GAM_triple_vil <- c(N_vmog_vil, N_fmog_vil, N_fvog_vil, N_fvmg_vil)
names(GAM_triple_vil) <- c( rep("P_vmog", N_village), rep("P_fmog", N_village), rep("P_fvog", N_village), rep("P_fvmg", N_village) )

GAM_triple_vil_obs      <- GAM_triple_vil
GAM_triple_vil_obs_low  <- GAM_triple_vil_obs
GAM_triple_vil_obs_high <- GAM_triple_vil_obs

for(v in 1:N_village)
{
	GAM_triple_vil_obs[0*N_village + v] <- GAM_triple_vil[0*N_village + v]/N_data_vil[v]
	GAM_triple_vil_obs[1*N_village + v] <- GAM_triple_vil[1*N_village + v]/N_data_vil[v]
	GAM_triple_vil_obs[2*N_village + v] <- GAM_triple_vil[2*N_village + v]/N_data_vil[v]
	GAM_triple_vil_obs[3*N_village + v] <- GAM_triple_vil[3*N_village + v]/N_data_vil[v]
}



for(v in 1:N_village)
{
	GAM_triple_vil_obs_low[0*N_village + v] <- binom.confint( GAM_triple_vil[0*N_village + v], N_data_vil[v], method="wilson")[1,5]
	GAM_triple_vil_obs_low[1*N_village + v] <- binom.confint( GAM_triple_vil[1*N_village + v], N_data_vil[v], method="wilson")[1,5]
	GAM_triple_vil_obs_low[2*N_village + v] <- binom.confint( GAM_triple_vil[2*N_village + v], N_data_vil[v], method="wilson")[1,5]
	GAM_triple_vil_obs_low[3*N_village + v] <- binom.confint( GAM_triple_vil[3*N_village + v], N_data_vil[v], method="wilson")[1,5]

	GAM_triple_vil_obs_high[0*N_village + v] <- binom.confint( GAM_triple_vil[0*N_village + v], N_data_vil[v], method="wilson")[1,6]
	GAM_triple_vil_obs_high[1*N_village + v] <- binom.confint( GAM_triple_vil[1*N_village + v], N_data_vil[v], method="wilson")[1,6]
	GAM_triple_vil_obs_high[2*N_village + v] <- binom.confint( GAM_triple_vil[2*N_village + v], N_data_vil[v], method="wilson")[1,6]
	GAM_triple_vil_obs_high[3*N_village + v] <- binom.confint( GAM_triple_vil[3*N_village + v], N_data_vil[v], method="wilson")[1,6]
}





GAM_triple_vil_exp <- GAM_triple_vil_obs

for(v in 1:N_village)
{
	GAM_triple_vil_exp[0*N_village + v] <- N_vg_vil[v]*N_mg_vil[v]*N_og_vil[v]/( N_data_vil[v]^3 )
	GAM_triple_vil_exp[1*N_village + v] <- N_fg_vil[v]*N_mg_vil[v]*N_og_vil[v]/( N_data_vil[v]^3 )
	GAM_triple_vil_exp[2*N_village + v] <- N_fg_vil[v]*N_vg_vil[v]*N_og_vil[v]/( N_data_vil[v]^3 )
	GAM_triple_vil_exp[3*N_village + v] <- N_fg_vil[v]*N_vg_vil[v]*N_mg_vil[v]/( N_data_vil[v]^3 )
}


#########################################
## Regression model

GAM_triple_vil_mod <- lm( GAM_triple_vil_obs ~ 0 + GAM_triple_vil_exp )



GAM_triple_vil_mod_fit <- predict.lm( GAM_triple_vil_mod, data.frame(GAM_triple_vil_exp=prev_seq), interval="prediction", se.fit=TRUE)

GAM_triple_vil_mod_fit$fit[which(GAM_triple_vil_mod_fit$fit[,2] < 0),2] <- 1e-10


c( summary(GAM_double_vil_mod)$coef[1],
   summary(GAM_double_vil_mod)$coef[1] - 1.96*summary(GAM_double_vil_mod)$coef[2], 
   summary(GAM_double_vil_mod)$coef[1] + 1.96*summary(GAM_double_vil_mod)$coef[2] )^2

c( summary(GAM_triple_vil_mod)$coef[1],
   summary(GAM_triple_vil_mod)$coef[1] - 1.96*summary(GAM_triple_vil_mod)$coef[2], 
   summary(GAM_triple_vil_mod)$coef[1] + 1.96*summary(GAM_triple_vil_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_mod)$df[2])





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



main.size = 1.4
axis.size = 0.7
lab.size  = 1.5

point.size = 1
line.size  = 1
rr.size1    = 0.005
rr.size2    = 0.0075
rr.size3    = 0.01
arrow.edge = 0.015
arrow.size = 0.5



tiff( file="Supp_Figure2_coinfection_prevalence_All.tif", width=18, height=18, units="cm", res=500)


lay.mat = rbind( c( 1,  3 ),
                 c( 2,  4 ),
                 c( 5,  5 ) )

layout(lay.mat, heights=c(12,12,1))
layout.show(5)

par(mar = c(3,3,2,0.5))

par(mgp = c(1.85, 0.6, 0))


ff <- 0.25  ## scaling factor for plotting


########################################
##                                    ##
##  PANEL 1: Asexual co-infection     ##
##          	                    ## 
########################################


line_seq_x <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff

line_seq_y <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff



plot(x=100, y=100, 
xlim=c(0,0.5^ff), ylim=c(0,0.5^ff), 
xlab="expected prevalence", ylab="observed prevalence",
main="(A) Asexual co-infection",
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


points(x=c(0,1), y=c(0,1), 
type='l', lwd=1.5, lty="dotted")




points(x=prev_seq^ff, y=ASX_double_vil_mod_fit$fit[,1]^ff, type='l', lwd=1.5, lty="solid")

polygon(x=c(prev_seq^ff, rev(prev_seq^ff)), 
	  y=c( ASX_double_vil_mod_fit$fit[,2]^ff, rev(ASX_double_vil_mod_fit$fit[,3]^ff) ),
	  col=rgb(100/256,100/256,100/256,0.1), border=NA)


################################
## Pf single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=ASX_double_exp_vil[v,1,1]^ff, yy=ASX_double_obs_vil[v,1,1]^ff, rr=rr.size1, 
	                     left_colour=Pf_col, right_colour=Pf_col)

	arrows(x0=ASX_double_exp_vil[v,1,1]^ff, y0=ASX_double_obs_low_vil[v,1,1]^ff, 
	     	 x1=ASX_double_exp_vil[v,1,1]^ff, y1=ASX_double_obs_high_vil[v,1,1]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	
}

################################
## Pv single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=ASX_double_exp_vil[v,2,2]^ff, yy=ASX_double_obs_vil[v,2,2]^ff, rr=rr.size1, 
      	               left_colour=Pv_col, right_colour=Pv_col)

	arrows(x0=ASX_double_exp_vil[v,2,2]^ff, y0=ASX_double_obs_low_vil[v,2,2]^ff, 
	     	 x1=ASX_double_exp_vil[v,2,2]^ff, y1=ASX_double_obs_high_vil[v,2,2]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pv_col, lwd=arrow.size)	
}

################################
## Pm single infection

for(v in 1:N_village)
{	
	semi_circle_plotter( xx=ASX_double_exp_vil[v,3,3]^ff, yy=ASX_double_obs_vil[v,3,3]^ff, rr=rr.size1, 
	                     left_colour=Pm_col, right_colour=Pm_col)

	arrows(x0=ASX_double_exp_vil[v,3,3]^ff, y0=ASX_double_obs_low_vil[v,3,3]^ff, 
	     	 x1=ASX_double_exp_vil[v,3,3]^ff, y1=ASX_double_obs_high_vil[v,3,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}

################################
## Po single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=ASX_double_exp_vil[v,4,4]^ff, yy=ASX_double_obs_vil[v,4,4]^ff, rr=rr.size1, 
	                     left_colour=Po_col, right_colour=Po_col)

	arrows(x0=ASX_double_exp_vil[v,4,4]^ff, y0=ASX_double_obs_low_vil[v,4,4]^ff, 
	     	 x1=ASX_double_exp_vil[v,4,4]^ff, y1=ASX_double_obs_high_vil[v,4,4]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Po_col, lwd=arrow.size)	
}


################################
## Pf + Pv co-infection

for(v in 1:N_village)
{
	arrows(x0=ASX_double_exp_vil[v,1,2]^ff, y0=ASX_double_obs_low_vil[v,1,2]^ff, 
	     	 x1=ASX_double_exp_vil[v,1,2]^ff, y1=ASX_double_obs_high_vil[v,1,2]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,1,2]^ff, yy=ASX_double_obs_vil[v,1,2]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Pv_col)
}

################################
## Pf + Pm co-infection

for(v in 1:N_village)
{	
	arrows(x0=ASX_double_exp_vil[v,1,3]^ff, y0=ASX_double_obs_low_vil[v,1,3]^ff, 
	     	 x1=ASX_double_exp_vil[v,1,3]^ff, y1=ASX_double_obs_high_vil[v,1,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,1,3]^ff, yy=ASX_double_obs_vil[v,1,3]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Pm_col)
}


################################
## Pf + Po co-infection

for(v in 1:N_village)
{	
	arrows(x0=ASX_double_exp_vil[v,1,4]^ff, y0=ASX_double_obs_low_vil[v,1,4]^ff, 
	     	 x1=ASX_double_exp_vil[v,1,4]^ff, y1=ASX_double_obs_high_vil[v,1,4]^ff, 
	    	 length=0.03, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,1,4]^ff, yy=ASX_double_obs_vil[v,1,4]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Po_col)
}

################################
## Pv + Pm co-infection

for(v in 1:N_village)
{
	arrows(x0=ASX_double_exp_vil[v,2,3]^ff, y0=ASX_double_obs_low_vil[v,2,3]^ff, 
	     	 x1=ASX_double_exp_vil[v,2,3]^ff, y1=ASX_double_obs_high_vil[v,2,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,2,3]^ff, yy=ASX_double_obs_vil[v,2,3]^ff, rr=rr.size2, 
	                     left_colour=Pv_col, right_colour=Pm_col)
}

################################
## Pv + Po co-infection

for(v in 1:N_village)
{
	arrows(x0=ASX_double_exp_vil[v,2,4]^ff, y0=ASX_double_obs_low_vil[v,2,4]^ff, 
	     	 x1=ASX_double_exp_vil[v,2,4]^ff, y1=ASX_double_obs_high_vil[v,2,4]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,2,4]^ff, yy=ASX_double_obs_vil[v,2,4]^ff, rr=rr.size2, 
	                     left_colour=Pv_col, right_colour=Po_col)
}

################################
## Pm + Po co-infection

for(v in 1:N_village)
{
	arrows(x0=ASX_double_exp_vil[v,3,4]^ff, y0=ASX_double_obs_low_vil[v,3,4]^ff, 
     		 x1=ASX_double_exp_vil[v,3,4]^ff, y1=ASX_double_obs_high_vil[v,3,4]^ff, 
    		 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=ASX_double_exp_vil[v,3,4]^ff, yy=ASX_double_obs_vil[v,3,4]^ff, rr=rr.size2, 
	                     left_colour=Pm_col, right_colour=Po_col)
}


axis(1,  at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%" ),
cex.axis=axis.size )

axis(1,  at = c(0.01)^ff, 
     labels = c("1%" ),
cex.axis=axis.size )

axis(2, las=2, at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%" ),
cex.axis=axis.size )





########################################
##                                    ##
##  PANEL 2: Asexual triple infection ##
##          	                    ## 
########################################


line_seq_x <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff

line_seq_y <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff



plot(x=100, y=100, 
xlim=c(0,0.1^ff), ylim=c(0,0.1^ff), 
xlab="expected prevalence", ylab="observed prevalence",
main="(B) Asexual triple infection",
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


points(x=c(0,1), y=c(0,1), 
type='l', lwd=1.5, lty="dotted")




points(x=prev_seq^ff, y=ASX_triple_vil_mod_fit$fit[,1]^ff, type='l', lwd=1.5, lty="solid")

polygon(x=c(prev_seq^ff, rev(prev_seq^ff)), 
	  y=c( ASX_triple_vil_mod_fit$fit[,2]^ff, rev(ASX_triple_vil_mod_fit$fit[,3]^ff) ),
	  col=rgb(100/256,100/256,100/256,0.1), border=NA)


################################
## Pv, Pm, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=ASX_triple_vil_exp[0*N_village+v]^ff, yy=ASX_triple_vil_obs[0*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pv_col, colour_2=Pm_col, colour_3=Po_col)

	arrows(x0=ASX_triple_vil_exp[0*N_village+v]^ff, y0=ASX_triple_vil_obs_low[0*N_village+v]^ff, 
	     	 x1=ASX_triple_vil_exp[0*N_village+v]^ff, y1=ASX_triple_vil_obs_high[0*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}


################################
## Pf, Pm, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=ASX_triple_vil_exp[1*N_village+v]^ff, yy=ASX_triple_vil_obs[1*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pm_col, colour_3=Po_col)

	arrows(x0=ASX_triple_vil_exp[1*N_village+v]^ff, y0=ASX_triple_vil_obs_low[1*N_village+v]^ff, 
	     	 x1=ASX_triple_vil_exp[1*N_village+v]^ff, y1=ASX_triple_vil_obs_high[1*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}




################################
## Pf, Pv, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=ASX_triple_vil_exp[2*N_village+v]^ff, yy=ASX_triple_vil_obs[2*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pv_col, colour_3=Po_col)

	arrows(x0=ASX_triple_vil_exp[2*N_village+v]^ff, y0=ASX_triple_vil_obs_low[2*N_village+v]^ff, 
	     	 x1=ASX_triple_vil_exp[2*N_village+v]^ff, y1=ASX_triple_vil_obs_high[2*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}




################################
## Pf, Pv, Pm tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=ASX_triple_vil_exp[3*N_village+v]^ff, yy=ASX_triple_vil_obs[3*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pv_col, colour_3=Pm_col)

	arrows(x0=ASX_triple_vil_exp[3*N_village+v]^ff, y0=ASX_triple_vil_obs_low[3*N_village+v]^ff, 
	     	 x1=ASX_triple_vil_exp[3*N_village+v]^ff, y1=ASX_triple_vil_obs_high[3*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}






axis(1,  at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%" ),
cex.axis=axis.size )

axis(1,  at = c(0.01)^ff, 
     labels = c("1%" ),
cex.axis=axis.size )

axis(2, las=2, at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%" ),
cex.axis=axis.size )




########################################
##                                    ##
##  PANEL 3: Gametocyte co-infection  ##
##          	                    ## 
########################################


line_seq_x <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff

line_seq_y <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff



plot(x=100, y=100, 
xlim=c(0,0.5^ff), ylim=c(0,0.5^ff), 
xlab="expected prevalence", ylab="observed prevalence",
main="(C) Gametocyte co-infection",
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


points(x=c(0,1), y=c(0,1), 
type='l', lwd=1.5, lty="dotted")




points(x=prev_seq^ff, y=GAM_double_vil_mod_fit$fit[,1]^ff, type='l', lwd=1.5, lty="solid")

polygon(x=c(prev_seq^ff, rev(prev_seq^ff)), 
	  y=c( GAM_double_vil_mod_fit$fit[,2]^ff, rev(GAM_double_vil_mod_fit$fit[,3]^ff) ),
	  col=rgb(100/256,100/256,100/256,0.1), border=NA)


################################
## Pf single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=GAM_double_exp_vil[v,1,1]^ff, yy=GAM_double_obs_vil[v,1,1]^ff, rr=rr.size1, 
	                     left_colour=Pf_col, right_colour=Pf_col)

	arrows(x0=GAM_double_exp_vil[v,1,1]^ff, y0=GAM_double_obs_low_vil[v,1,1]^ff, 
	     	 x1=GAM_double_exp_vil[v,1,1]^ff, y1=GAM_double_obs_high_vil[v,1,1]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pf_col, lwd=arrow.size)	
}

################################
## Pv single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=GAM_double_exp_vil[v,2,2]^ff, yy=GAM_double_obs_vil[v,2,2]^ff, rr=rr.size1, 
      	               left_colour=Pv_col, right_colour=Pv_col)

	arrows(x0=GAM_double_exp_vil[v,2,2]^ff, y0=GAM_double_obs_low_vil[v,2,2]^ff, 
	     	 x1=GAM_double_exp_vil[v,2,2]^ff, y1=GAM_double_obs_high_vil[v,2,2]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pv_col, lwd=arrow.size)	
}

################################
## Pm single infection

for(v in 1:N_village)
{	
	semi_circle_plotter( xx=GAM_double_exp_vil[v,3,3]^ff, yy=GAM_double_obs_vil[v,3,3]^ff, rr=rr.size1, 
	                     left_colour=Pm_col, right_colour=Pm_col)

	arrows(x0=GAM_double_exp_vil[v,3,3]^ff, y0=GAM_double_obs_low_vil[v,3,3]^ff, 
	     	 x1=GAM_double_exp_vil[v,3,3]^ff, y1=GAM_double_obs_high_vil[v,3,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Pm_col, lwd=arrow.size)	
}

################################
## Po single infection

for(v in 1:N_village)
{
	semi_circle_plotter( xx=GAM_double_exp_vil[v,4,4]^ff, yy=GAM_double_obs_vil[v,4,4]^ff, rr=rr.size1, 
	                     left_colour=Po_col, right_colour=Po_col)

	arrows(x0=GAM_double_exp_vil[v,4,4]^ff, y0=GAM_double_obs_low_vil[v,4,4]^ff, 
	     	 x1=GAM_double_exp_vil[v,4,4]^ff, y1=GAM_double_obs_high_vil[v,4,4]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col=Po_col, lwd=arrow.size)	
}


################################
## Pf + Pv co-infection

for(v in 1:N_village)
{
	arrows(x0=GAM_double_exp_vil[v,1,2]^ff, y0=GAM_double_obs_low_vil[v,1,2]^ff, 
	     	 x1=GAM_double_exp_vil[v,1,2]^ff, y1=GAM_double_obs_high_vil[v,1,2]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,1,2]^ff, yy=GAM_double_obs_vil[v,1,2]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Pv_col)
}

################################
## Pf + Pm co-infection

for(v in 1:N_village)
{	
	arrows(x0=GAM_double_exp_vil[v,1,3]^ff, y0=GAM_double_obs_low_vil[v,1,3]^ff, 
	     	 x1=GAM_double_exp_vil[v,1,3]^ff, y1=GAM_double_obs_high_vil[v,1,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,1,3]^ff, yy=GAM_double_obs_vil[v,1,3]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Pm_col)
}


################################
## Pf + Po co-infection

for(v in 1:N_village)
{	
	arrows(x0=GAM_double_exp_vil[v,1,4]^ff, y0=GAM_double_obs_low_vil[v,1,4]^ff, 
	     	 x1=GAM_double_exp_vil[v,1,4]^ff, y1=GAM_double_obs_high_vil[v,1,4]^ff, 
	    	 length=0.03, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,1,4]^ff, yy=GAM_double_obs_vil[v,1,4]^ff, rr=rr.size2, 
	                     left_colour=Pf_col, right_colour=Po_col)
}

################################
## Pv + Pm co-infection

for(v in 1:N_village)
{
	arrows(x0=GAM_double_exp_vil[v,2,3]^ff, y0=GAM_double_obs_low_vil[v,2,3]^ff, 
	     	 x1=GAM_double_exp_vil[v,2,3]^ff, y1=GAM_double_obs_high_vil[v,2,3]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,2,3]^ff, yy=GAM_double_obs_vil[v,2,3]^ff, rr=rr.size2, 
	                     left_colour=Pv_col, right_colour=Pm_col)
}

################################
## Pv + Po co-infection

for(v in 1:N_village)
{
	arrows(x0=GAM_double_exp_vil[v,2,4]^ff, y0=GAM_double_obs_low_vil[v,2,4]^ff, 
	     	 x1=GAM_double_exp_vil[v,2,4]^ff, y1=GAM_double_obs_high_vil[v,2,4]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,2,4]^ff, yy=GAM_double_obs_vil[v,2,4]^ff, rr=rr.size2, 
	                     left_colour=Pv_col, right_colour=Po_col)
}

################################
## Pm + Po co-infection

for(v in 1:N_village)
{
	arrows(x0=GAM_double_exp_vil[v,3,4]^ff, y0=GAM_double_obs_low_vil[v,3,4]^ff, 
     		 x1=GAM_double_exp_vil[v,3,4]^ff, y1=GAM_double_obs_high_vil[v,3,4]^ff, 
    		 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	

	semi_circle_plotter( xx=GAM_double_exp_vil[v,3,4]^ff, yy=GAM_double_obs_vil[v,3,4]^ff, rr=rr.size2, 
	                     left_colour=Pm_col, right_colour=Po_col)
}


axis(1,  at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%" ),
cex.axis=axis.size )

axis(1,  at = c(0.01)^ff, 
     labels = c("1%" ),
cex.axis=axis.size )

axis(2, las=2, at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%", "20%", "50%" ),
cex.axis=axis.size )





###########################################
##                                       ##
##  PANEL 4: Gametocyte triple infection ##
##          	                       ## 
###########################################


line_seq_x <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff

line_seq_y <- c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)^ff



plot(x=100, y=100, 
xlim=c(0,0.1^ff), ylim=c(0,0.1^ff), 
xlab="expected prevalence", ylab="observed prevalence",
main="(D) Gametocyte triple infection",
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


points(x=c(0,1), y=c(0,1), 
type='l', lwd=1.5, lty="dotted")




points(x=prev_seq^ff, y=GAM_triple_vil_mod_fit$fit[,1]^ff, type='l', lwd=1.5, lty="solid")

polygon(x=c(prev_seq^ff, rev(prev_seq^ff)), 
	  y=c( GAM_triple_vil_mod_fit$fit[,2]^ff, rev(GAM_triple_vil_mod_fit$fit[,3]^ff) ),
	  col=rgb(100/256,100/256,100/256,0.1), border=NA)


################################
## Pv, Pm, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=GAM_triple_vil_exp[0*N_village+v]^ff, yy=GAM_triple_vil_obs[0*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pv_col, colour_2=Pm_col, colour_3=Po_col)

	arrows(x0=GAM_triple_vil_exp[0*N_village+v]^ff, y0=GAM_triple_vil_obs_low[0*N_village+v]^ff, 
	     	 x1=GAM_triple_vil_exp[0*N_village+v]^ff, y1=GAM_triple_vil_obs_high[0*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}


################################
## Pf, Pm, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=GAM_triple_vil_exp[1*N_village+v]^ff, yy=GAM_triple_vil_obs[1*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pm_col, colour_3=Po_col)

	arrows(x0=GAM_triple_vil_exp[1*N_village+v]^ff, y0=GAM_triple_vil_obs_low[1*N_village+v]^ff, 
	     	 x1=GAM_triple_vil_exp[1*N_village+v]^ff, y1=GAM_triple_vil_obs_high[1*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}




################################
## Pf, Pv, Po tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=GAM_triple_vil_exp[2*N_village+v]^ff, yy=GAM_triple_vil_obs[2*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pv_col, colour_3=Po_col)

	arrows(x0=GAM_triple_vil_exp[2*N_village+v]^ff, y0=GAM_triple_vil_obs_low[2*N_village+v]^ff, 
	     	 x1=GAM_triple_vil_exp[2*N_village+v]^ff, y1=GAM_triple_vil_obs_high[2*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}




################################
## Pf, Pv, Pm tripe infection

for(v in 1:N_village)
{
	tri_circle_plotter( xx=GAM_triple_vil_exp[3*N_village+v]^ff, yy=GAM_triple_vil_obs[3*N_village+v]^ff, rr=rr.size2, 
	                     colour_1=Pf_col, colour_2=Pv_col, colour_3=Pm_col)

	arrows(x0=GAM_triple_vil_exp[3*N_village+v]^ff, y0=GAM_triple_vil_obs_low[3*N_village+v]^ff, 
	     	 x1=GAM_triple_vil_exp[3*N_village+v]^ff, y1=GAM_triple_vil_obs_high[3*N_village+v]^ff, 
	    	 length=arrow.edge, angle=90, code=3, col="black", lwd=arrow.size)	
}






axis(1,  at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1)^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%" ),
cex.axis=axis.size )

axis(1,  at = c(0.01)^ff, 
     labels = c("1%" ),
cex.axis=axis.size )

axis(2, las=2, at = c(0, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1 )^ff, 
     labels = c( "0%", "0.01%", "0.1%", "0.5%", "1%", "2%", "5%", "10%" ),
cex.axis=axis.size )











########################################
##                                    ##
##  PANEL 5: Legend                   ##
##          	                    ## 
########################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("P. falciparum", "P. vivax", "P. malariae", "P. ovale"),
       fill = c(Pf_col, Pv_col, Pm_col, Po_col), 
       border = c(Pf_col, Pv_col, Pm_col, Po_col), 
       ncol=4, cex=1.5, bty="n" )



dev.off()















