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
##  1.1. All data                      ##
##                                     ##
#########################################
#########################################

N_f <- length( which(ALB_par$pf_pos==1) )
N_v <- length( which(ALB_par$pv_pos==1) )
N_m <- length( which(ALB_par$pm_pos==1) )
N_o <- length( which(ALB_par$po_pos==1) )

N_fv <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1) )
N_fm <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1) )
N_fo <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1) )
N_vm <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1) )
N_vo <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1) )
N_mo <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1) )

N_fvm <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1) )
N_fvo <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1) )
N_fmo <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )
N_vmo <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )

N_fvmo <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )


#########################################
##                                     ##
##  1.1.1. All data - double           ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double <- rbind( c(N_f,  N_fv, N_fm, N_fo),
                     c(N_fv, N_v,  N_vm, N_vo),
                     c(N_fm, N_vm, N_m,  N_mo),
                     c(N_fo, N_vo, N_mo, N_o ) )

colnames(ASX_double) <- c("Pf", "Pv", "Po", "Pm")
rownames(ASX_double) <- c("Pf", "Pv", "Po", "Pm")

ASX_double_obs <- ASX_double/N_data

ASX_double_exp <- (diag(ASX_double)/N_data) %o% (diag(ASX_double)/N_data)
diag(ASX_double_exp) <- diag(ASX_double_obs)

ASX_double_obs_low <- ASX_double
ASX_double_obs_high <- ASX_double

for(i in 1:4)
{
	for(j in 1:4)
	{
		ASX_double_obs_low[i,j]  <- binom.confint( ASX_double[i,j], N_data, method="wilson")[1,5]
		ASX_double_obs_high[i,j] <- binom.confint( ASX_double[i,j], N_data, method="wilson")[1,6]
	}
}


#########################################
## Regression model

x_exp <- c( ASX_double_exp[1,2:4], ASX_double_exp[2,3:4], ASX_double_exp[3,4] )
y_obs <- c( ASX_double_obs[1,2:4], ASX_double_obs[2,3:4], ASX_double_obs[3,4] )

ASX_double_mod <- lm( y_obs ~ 0 + x_exp )



c( summary(ASX_double_mod)$coef[1],
   summary(ASX_double_mod)$coef[1] - 1.96*summary(ASX_double_mod)$coef[2], 
   summary(ASX_double_mod)$coef[1] + 1.96*summary(ASX_double_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_mod)$coef[1], 
         sigma = summary(ASX_double_mod)$coef[2], 
         nu    = summary(ASX_double_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_mod)$coef[1], 
         sigma = summary(ASX_double_mod)$coef[2], 
         nu    = summary(ASX_double_mod)$df[2])




#########################################
##                                     ##
##  1.1.2. All data - triple           ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

ASX_triple <- c(N_vmo, N_fmo, N_fvo, N_fvm)
names(ASX_triple) <- c("P_vmo", "P_fmo", "P_fvo", "P_fvm")

ASX_triple_obs <- ASX_triple/N_data

ASX_triple_obs_low  <- ASX_triple_obs
ASX_triple_obs_high <- ASX_triple_obs

for(i in 1:4)
{
	ASX_triple_obs_low[i]  <- binom.confint( ASX_triple[i], N_data, method="wilson")[1,5]
	ASX_triple_obs_high[i] <- binom.confint( ASX_triple[i], N_data, method="wilson")[1,6]
}


ASX_triple_exp <- ASX_triple_obs

ASX_triple_exp[1] <- N_v*N_m*N_o/( N_data^3 )
ASX_triple_exp[2] <- N_f*N_m*N_o/( N_data^3 )
ASX_triple_exp[3] <- N_f*N_v*N_m/( N_data^3 )
ASX_triple_exp[4] <- N_f*N_v*N_o/( N_data^3 )



#########################################
## Regression model

ASX_triple_mod <- lm( ASX_triple_obs ~ 0 + ASX_triple_exp )


c( summary(ASX_double_mod)$coef[1],
   summary(ASX_double_mod)$coef[1] - 1.96*summary(ASX_double_mod)$coef[2], 
   summary(ASX_double_mod)$coef[1] + 1.96*summary(ASX_double_mod)$coef[2] )^2

c( summary(ASX_triple_mod)$coef[1],
   summary(ASX_triple_mod)$coef[1] - 1.96*summary(ASX_triple_mod)$coef[2], 
   summary(ASX_triple_mod)$coef[1] + 1.96*summary(ASX_triple_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])

2*pst( q     = summary(ASX_double_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])



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
##                                     ##
##  1.3. All data                      ##
##       Village & treat stratified    ##
##                                     ##
#########################################
#########################################

N_data_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_f_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_v_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_m_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_o_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fv_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fm_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vm_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_mo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fvm_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fvo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fmo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vmo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fvmo_vil_trt <- matrix(NA, nrow=2, ncol=N_village)


for(v in 1:N_village)
{
	N_data_vil_trt[1,v] <- length(which(ALB_par$village == village_names[v] & ALB_par$treat=="PLACEBO"))

	N_f_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_v_vil_trt[1,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_m_vil_trt[1,v] <- length( which(ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_o_vil_trt[1,v] <- length( which(ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fv_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fm_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fo_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vm_vil_trt[1,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vo_vil_trt[1,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_mo_vil_trt[1,v] <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fvm_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fvo_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fmo_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vmo_vil_trt[1,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fvmo_vil_trt[1,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
}


for(v in 1:N_village)
{
	N_data_vil_trt[2,v] <- length(which(ALB_par$village == village_names[v] & ALB_par$treat=="PRMQN"))

	N_f_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_v_vil_trt[2,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_m_vil_trt[2,v] <- length( which(ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_o_vil_trt[2,v] <- length( which(ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fv_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fm_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fo_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vm_vil_trt[2,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vo_vil_trt[2,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_mo_vil_trt[2,v] <- length( which(ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fvm_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fvo_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fmo_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vmo_vil_trt[2,v] <- length( which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fvmo_vil_trt[2,v] <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
}




#########################################
##                                     ##
##  1.3.1. All data - double           ##
##       Village & treat stratified    ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double_vil_trt     <- array(NA, dim=c(2, N_village, 4, 4))
ASX_double_obs_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))
ASX_double_exp_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))

ASX_double_obs_low_vil_trt  <- array(NA, dim=c(2, N_village, 4, 4))
ASX_double_obs_high_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))


for(trt in 1:2)
{
	for(v in 1:N_village)
	{
		ASX_double_vil_trt[trt,v,,] <- rbind( c(N_f_vil_trt[trt,v],  N_fv_vil_trt[trt,v], N_fm_vil_trt[trt,v], N_fo_vil_trt[trt,v]),
	      	 	                          c(N_fv_vil_trt[trt,v], N_v_vil_trt[trt,v],  N_vm_vil_trt[trt,v], N_vo_vil_trt[trt,v]),
	            		                    c(N_fm_vil_trt[trt,v], N_vm_vil_trt[trt,v], N_m_vil_trt[trt,v],  N_mo_vil_trt[trt,v]),
	              	                          c(N_fo_vil_trt[trt,v], N_vo_vil_trt[trt,v], N_mo_vil_trt[trt,v], N_o_vil_trt[trt,v] ) )

		colnames(ASX_double_vil_trt[trt,v,,]) <- c("Pf", "Pv", "Po", "Pm")
		rownames(ASX_double_vil_trt[trt,v,,]) <- c("Pf", "Pv", "Po", "Pm")

		ASX_double_obs_vil_trt[trt,v,,] <- ASX_double_vil_trt[trt,v,,]/N_data_vil_trt[trt,v]

		ASX_double_exp_vil_trt[trt,v,,] <- (diag(ASX_double_vil_trt[trt,v,,])/N_data_vil_trt[trt,v]) %o% (diag(ASX_double_vil_trt[trt,v,,])/N_data_vil_trt[trt,v])
		diag(ASX_double_exp_vil_trt[trt,v,,]) <- diag(ASX_double_obs_vil_trt[trt,v,,])


		ASX_double_obs_low_vil_trt[trt,v,,]  <- ASX_double_vil_trt[trt,v,,]
		ASX_double_obs_high_vil_trt[trt,v,,] <- ASX_double_vil_trt[trt,v,,]

		for(i in 1:4)
		{
			for(j in 1:4)
			{
				ASX_double_obs_low_vil_trt[trt,v,i,j]  <- binom.confint( ASX_double_vil_trt[trt,v,i,j], N_data_vil_trt[trt,v], method="wilson")[1,5]
				ASX_double_obs_high_vil_trt[trt,v,i,j] <- binom.confint( ASX_double_vil_trt[trt,v,i,j], N_data_vil_trt[trt,v], method="wilson")[1,6]
			}
		}

	}
}

#########################################
## Regression model

x_exp <- c( ASX_double_exp_vil_trt[,,1,2:4], ASX_double_exp_vil_trt[,,2,3:4], ASX_double_exp_vil_trt[,,3,4] )
y_obs <- c( ASX_double_obs_vil_trt[,,1,2:4], ASX_double_obs_vil_trt[,,2,3:4], ASX_double_obs_vil_trt[,,3,4] )

ASX_double_vil_trt_mod <- lm( y_obs ~ 0 + x_exp )





c( summary(ASX_double_vil_trt_mod)$coef[1],
   summary(ASX_double_vil_trt_mod)$coef[1] - 1.96*summary(ASX_double_vil_trt_mod)$coef[2], 
   summary(ASX_double_vil_trt_mod)$coef[1] + 1.96*summary(ASX_double_vil_trt_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_trt_mod)$coef[1], 
         sigma = summary(ASX_double_vil_trt_mod)$coef[2], 
         nu    = summary(ASX_double_vil_trt_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_trt_mod)$coef[1], 
         sigma = summary(ASX_double_vil_trt_mod)$coef[2], 
         nu    = summary(ASX_double_vil_trt_mod)$df[2])



#########################################
##                                     ##
##  1.3.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

ASX_triple_vil_trt <- c(N_vmo_vil_trt, N_fmo_vil_trt, N_fvo_vil_trt, N_fvm_vil_trt)
names(ASX_triple_vil_trt) <- c( rep("P_vmo", 2*N_village), rep("P_fmo", 2*N_village), rep("P_fvo", 2*N_village), rep("P_fvm", 2*N_village) )


ASX_triple_vil_trt_obs      <- ASX_triple_vil_trt
ASX_triple_vil_trt_obs_low  <- ASX_triple_vil_trt
ASX_triple_vil_trt_obs_high <- ASX_triple_vil_trt

for(v in 1:N_village)
{
	ASX_triple_vil_trt_obs[0*N_village + v] <- ASX_triple_vil_trt[0*N_village + v]/N_data_vil_trt[1,v]
	ASX_triple_vil_trt_obs[1*N_village + v] <- ASX_triple_vil_trt[1*N_village + v]/N_data_vil_trt[2,v]

	ASX_triple_vil_trt_obs[2*N_village + v] <- ASX_triple_vil_trt[2*N_village + v]/N_data_vil_trt[1,v]
	ASX_triple_vil_trt_obs[3*N_village + v] <- ASX_triple_vil_trt[3*N_village + v]/N_data_vil_trt[2,v]

	ASX_triple_vil_trt_obs[4*N_village + v] <- ASX_triple_vil_trt[4*N_village + v]/N_data_vil_trt[1,v]
	ASX_triple_vil_trt_obs[5*N_village + v] <- ASX_triple_vil_trt[5*N_village + v]/N_data_vil_trt[2,v]

	ASX_triple_vil_trt_obs[6*N_village + v] <- ASX_triple_vil_trt[6*N_village + v]/N_data_vil_trt[1,v]
	ASX_triple_vil_trt_obs[7*N_village + v] <- ASX_triple_vil_trt[7*N_village + v]/N_data_vil_trt[2,v]
}



for(v in 1:N_village)
{
	ASX_triple_vil_trt_obs_low[0*N_village + v] <- binom.confint( ASX_triple_vil_trt[0*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	ASX_triple_vil_trt_obs_low[1*N_village + v] <- binom.confint( ASX_triple_vil_trt[1*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	ASX_triple_vil_trt_obs_low[2*N_village + v] <- binom.confint( ASX_triple_vil_trt[2*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	ASX_triple_vil_trt_obs_low[3*N_village + v] <- binom.confint( ASX_triple_vil_trt[3*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	ASX_triple_vil_trt_obs_low[4*N_village + v] <- binom.confint( ASX_triple_vil_trt[4*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	ASX_triple_vil_trt_obs_low[5*N_village + v] <- binom.confint( ASX_triple_vil_trt[5*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	ASX_triple_vil_trt_obs_low[6*N_village + v] <- binom.confint( ASX_triple_vil_trt[6*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	ASX_triple_vil_trt_obs_low[7*N_village + v] <- binom.confint( ASX_triple_vil_trt[7*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]


	ASX_triple_vil_trt_obs_high[0*N_village + v] <- binom.confint( ASX_triple_vil_trt[0*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	ASX_triple_vil_trt_obs_high[1*N_village + v] <- binom.confint( ASX_triple_vil_trt[1*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	ASX_triple_vil_trt_obs_high[2*N_village + v] <- binom.confint( ASX_triple_vil_trt[2*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	ASX_triple_vil_trt_obs_high[3*N_village + v] <- binom.confint( ASX_triple_vil_trt[3*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	ASX_triple_vil_trt_obs_high[4*N_village + v] <- binom.confint( ASX_triple_vil_trt[4*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	ASX_triple_vil_trt_obs_high[5*N_village + v] <- binom.confint( ASX_triple_vil_trt[5*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	ASX_triple_vil_trt_obs_high[6*N_village + v] <- binom.confint( ASX_triple_vil_trt[6*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	ASX_triple_vil_trt_obs_high[7*N_village + v] <- binom.confint( ASX_triple_vil_trt[7*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]
}





ASX_triple_vil_trt_exp <- ASX_triple_vil_trt_obs

for(v in 1:N_village)
{
	ASX_triple_vil_trt_exp[0*N_village + v] <- N_v_vil_trt[1,v]*N_m_vil_trt[1,v]*N_o_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	ASX_triple_vil_trt_exp[1*N_village + v] <- N_v_vil_trt[2,v]*N_m_vil_trt[2,v]*N_o_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	ASX_triple_vil_trt_exp[2*N_village + v] <- N_f_vil_trt[1,v]*N_m_vil_trt[1,v]*N_o_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	ASX_triple_vil_trt_exp[3*N_village + v] <- N_f_vil_trt[2,v]*N_m_vil_trt[2,v]*N_o_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	ASX_triple_vil_trt_exp[4*N_village + v] <- N_f_vil_trt[1,v]*N_v_vil_trt[1,v]*N_o_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	ASX_triple_vil_trt_exp[5*N_village + v] <- N_f_vil_trt[2,v]*N_v_vil_trt[2,v]*N_o_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	ASX_triple_vil_trt_exp[6*N_village + v] <- N_f_vil_trt[1,v]*N_v_vil_trt[1,v]*N_m_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	ASX_triple_vil_trt_exp[7*N_village + v] <- N_f_vil_trt[2,v]*N_v_vil_trt[2,v]*N_m_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )
}


#########################################
## Regression model

ASX_triple_vil_trt_mod <- lm( ASX_triple_vil_trt_obs ~ 0 + ASX_triple_vil_trt_exp )



c( summary(ASX_double_vil_trt_mod)$coef[1],
   summary(ASX_double_vil_trt_mod)$coef[1] - 1.96*summary(ASX_double_vil_trt_mod)$coef[2], 
   summary(ASX_double_vil_trt_mod)$coef[1] + 1.96*summary(ASX_double_vil_trt_mod)$coef[2] )^2

c( summary(ASX_triple_vil_trt_mod)$coef[1],
   summary(ASX_triple_vil_trt_mod)$coef[1] - 1.96*summary(ASX_triple_vil_trt_mod)$coef[2], 
   summary(ASX_triple_vil_trt_mod)$coef[1] + 1.96*summary(ASX_triple_vil_trt_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_trt_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])






#########################################
#########################################
##                                     ##
##  1.4. All data                      ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################
#########################################

#########################################
##                                     ##
##  1.4.1. All data - double           ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################

#########################################
## Regression model

x_exp <- c( ASX_double_exp_vil_trt[1,,1,2:4], ASX_double_exp_vil_trt[1,,2,3:4], ASX_double_exp_vil_trt[1,,3,4] )
y_obs <- c( ASX_double_obs_vil_trt[1,,1,2:4], ASX_double_obs_vil_trt[1,,2,3:4], ASX_double_obs_vil_trt[1,,3,4] )

ASX_double_vil_plac_mod <- lm( y_obs ~ 0 + x_exp )





c( summary(ASX_double_vil_plac_mod)$coef[1],
   summary(ASX_double_vil_plac_mod)$coef[1] - 1.96*summary(ASX_double_vil_plac_mod)$coef[2], 
   summary(ASX_double_vil_plac_mod)$coef[1] + 1.96*summary(ASX_double_vil_plac_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_plac_mod)$coef[1], 
         sigma = summary(ASX_double_vil_plac_mod)$coef[2], 
         nu    = summary(ASX_double_vil_plac_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_plac_mod)$coef[1], 
         sigma = summary(ASX_double_vil_plac_mod)$coef[2], 
         nu    = summary(ASX_double_vil_plac_mod)$df[2])



#########################################
##                                     ##
##  1.4.2. All data - triple           ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################

#########################################
## Regression model

ASX_triple_vil_plac_obs <- ASX_triple_vil_trt_obs[c(1:5,11:15,21:25,31:35)]
ASX_triple_vil_plac_exp <- ASX_triple_vil_trt_exp[c(1:5,11:15,21:25,31:35)]


ASX_triple_vil_plac_mod <- lm( ASX_triple_vil_plac_obs ~ 0 + ASX_triple_vil_plac_exp )

c( summary(ASX_double_vil_plac_mod)$coef[1],
   summary(ASX_double_vil_plac_mod)$coef[1] - 1.96*summary(ASX_double_vil_plac_mod)$coef[2], 
   summary(ASX_double_vil_plac_mod)$coef[1] + 1.96*summary(ASX_double_vil_plac_mod)$coef[2] )^2

c( summary(ASX_triple_vil_plac_mod)$coef[1],
   summary(ASX_triple_vil_plac_mod)$coef[1] - 1.96*summary(ASX_triple_vil_plac_mod)$coef[2], 
   summary(ASX_triple_vil_plac_mod)$coef[1] + 1.96*summary(ASX_triple_vil_plac_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_plac_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])






#########################################
#########################################
##                                     ##
##  1.5. Baseline data                 ##
##                                     ##
#########################################
#########################################

N_data_base <- length( which(ALB_par$study_day<0) )

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
##  1.5.1. All data - double           ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

ASX_double_base <- rbind( c(N_f_base,  N_fv_base, N_fm_base, N_fo_base),
                          c(N_fv_base, N_v_base,  N_vm_base, N_vo_base),
                          c(N_fm_base, N_vm_base, N_m_base,  N_mo_base),
                          c(N_fo_base, N_vo_base, N_mo_base, N_o_base ) )

colnames(ASX_double_base) <- c("Pf", "Pv", "Po", "Pm")
rownames(ASX_double_base) <- c("Pf", "Pv", "Po", "Pm")

ASX_double_base_obs <- ASX_double_base/N_data_base

ASX_double_base_exp <- (diag(ASX_double_base)/N_data_base) %o% (diag(ASX_double_base)/N_data_base)
diag(ASX_double_base_exp) <- diag(ASX_double_base_obs)

ASX_double_base_obs_low <- ASX_double_base
ASX_double_base_obs_high <- ASX_double_base

for(i in 1:4)
{
	for(j in 1:4)
	{
		ASX_double_base_obs_low[i,j]  <- binom.confint( ASX_double_base[i,j], N_data_base, method="wilson")[1,5]
		ASX_double_base_obs_high[i,j] <- binom.confint( ASX_double_base[i,j], N_data_base, method="wilson")[1,6]
	}
}


#########################################
## Regression model

x_exp <- c( ASX_double_base_exp[1,2:4], ASX_double_base_exp[2,3:4], ASX_double_base_exp[3,4] )
y_obs <- c( ASX_double_base_obs[1,2:4], ASX_double_base_obs[2,3:4], ASX_double_base_obs[3,4] )

ASX_double_base_mod <- lm( y_obs ~ 0 + x_exp )




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
##                                     ##
##  1.5.2. All data - triple           ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

ASX_triple_base <- c(N_vmo_base, N_fmo_base, N_fvo_base, N_fvm_base)
names(ASX_triple_base) <- c("P_vmo", "P_fmo", "P_fvo", "P_fvm")

ASX_triple_base_obs <- ASX_triple_base/N_data_base

ASX_triple_base_obs_low  <- ASX_triple_base
ASX_triple_base_obs_high <- ASX_triple_base

for(i in 1:4)
{
	ASX_triple_base_obs_low[i]  <- binom.confint( ASX_triple_base[i], N_data_base, method="wilson")[1,5]
	ASX_triple_base_obs_high[i] <- binom.confint( ASX_triple_base[i], N_data_base, method="wilson")[1,6]
}


ASX_triple_base_exp <- ASX_triple_base_obs

ASX_triple_base_exp[1] <- N_v_base*N_m_base*N_o_base/( N_data_base^3 )
ASX_triple_base_exp[2] <- N_f_base*N_m_base*N_o_base/( N_data_base^3 )
ASX_triple_base_exp[3] <- N_f_base*N_v_base*N_m_base/( N_data_base^3 )
ASX_triple_base_exp[4] <- N_f_base*N_v_base*N_o_base/( N_data_base^3 )



#########################################
## Regression model

ASX_triple_base_mod <- lm( ASX_triple_base_obs ~ 0 + ASX_triple_base_exp )



c( summary(ASX_double_base_mod)$coef[1],
   summary(ASX_double_base_mod)$coef[1] - 1.96*summary(ASX_double_base_mod)$coef[2], 
   summary(ASX_double_base_mod)$coef[1] + 1.96*summary(ASX_double_base_mod)$coef[2] )^2

c( summary(ASX_triple_base_mod)$coef[1],
   summary(ASX_triple_base_mod)$coef[1] - 1.96*summary(ASX_triple_base_mod)$coef[2], 
   summary(ASX_triple_base_mod)$coef[1] + 1.96*summary(ASX_triple_base_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])

2*pst( q     = summary(ASX_double_base_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])




#########################################
#########################################
##                                     ##
##  1.6. Baseline data                 ##
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
##  1.6.1. All data - double           ##
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



c( summary(ASX_double_vil_base_mod)$coef[1],
   summary(ASX_double_vil_base_mod)$coef[1] - 1.96*summary(ASX_double_vil_base_mod)$coef[2], 
   summary(ASX_double_vil_base_mod)$coef[1] + 1.96*summary(ASX_double_vil_base_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])



#########################################
##                                     ##
##  1.6.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

ASX_triple_vil_base <- c(N_vmo_vil_base, N_fmo_vil_base, N_fvo_vil_base, N_fvm_vil_base)
names(ASX_triple_vil_base) <- c( rep("P_vmo", N_village), rep("P_fmo", N_village), rep("P_fvo", N_village), rep("P_fvm", N_village) )


ASX_triple_vil_base_obs  <- ASX_triple_vil_base
ASX_triple_vil_base_obs_low  <- ASX_triple_vil_base
ASX_triple_vil_base_obs_high <- ASX_triple_vil_base


for(v in 1:N_village)
{
	ASX_triple_vil_base_obs[0*N_village + v] <- ASX_triple_vil_base[0*N_village + v]/N_data_vil_base[v]
	ASX_triple_vil_base_obs[1*N_village + v] <- ASX_triple_vil_base[1*N_village + v]/N_data_vil_base[v]
	ASX_triple_vil_base_obs[2*N_village + v] <- ASX_triple_vil_base[2*N_village + v]/N_data_vil_base[v]
	ASX_triple_vil_base_obs[3*N_village + v] <- ASX_triple_vil_base[3*N_village + v]/N_data_vil_base[v]
}

for(v in 1:N_village)
{
	ASX_triple_vil_base_obs_low[0*N_village + v] <- binom.confint( ASX_triple_vil_base[0*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	ASX_triple_vil_base_obs_low[1*N_village + v] <- binom.confint( ASX_triple_vil_base[1*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	ASX_triple_vil_base_obs_low[2*N_village + v] <- binom.confint( ASX_triple_vil_base[2*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	ASX_triple_vil_base_obs_low[3*N_village + v] <- binom.confint( ASX_triple_vil_base[3*N_village + v], N_data_vil_base[v], method="wilson")[1,5]

	ASX_triple_vil_base_obs_high[0*N_village + v] <- binom.confint( ASX_triple_vil_base[0*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	ASX_triple_vil_base_obs_high[1*N_village + v] <- binom.confint( ASX_triple_vil_base[1*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	ASX_triple_vil_base_obs_high[2*N_village + v] <- binom.confint( ASX_triple_vil_base[2*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	ASX_triple_vil_base_obs_high[3*N_village + v] <- binom.confint( ASX_triple_vil_base[3*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
}





ASX_triple_vil_base_exp <- ASX_triple_vil_base_obs

for(v in 1:N_village)
{
	ASX_triple_vil_base_exp[0*N_village + v] <- N_v_vil_base[v]*N_m_vil_base[v]*N_o_vil_base[v]/( N_data_vil_base[v]^3 )
	ASX_triple_vil_base_exp[1*N_village + v] <- N_f_vil_base[v]*N_m_vil_base[v]*N_o_vil_base[v]/( N_data_vil_base[v]^3 )
	ASX_triple_vil_base_exp[2*N_village + v] <- N_f_vil_base[v]*N_v_vil_base[v]*N_o_vil_base[v]/( N_data_vil_base[v]^3 )
	ASX_triple_vil_base_exp[3*N_village + v] <- N_f_vil_base[v]*N_v_vil_base[v]*N_m_vil_base[v]/( N_data_vil_base[v]^3 )
}


#########################################
## Regression model

ASX_triple_vil_base_mod <- lm( ASX_triple_vil_base_obs ~ 0 + ASX_triple_vil_base_exp )




c( summary(ASX_double_vil_base_mod)$coef[1],
   summary(ASX_double_vil_base_mod)$coef[1] - 1.96*summary(ASX_double_vil_base_mod)$coef[2], 
   summary(ASX_double_vil_base_mod)$coef[1] + 1.96*summary(ASX_double_vil_base_mod)$coef[2] )^2

c( summary(ASX_triple_vil_base_mod)$coef[1],
   summary(ASX_triple_vil_base_mod)$coef[1] - 1.96*summary(ASX_triple_vil_base_mod)$coef[2], 
   summary(ASX_triple_vil_base_mod)$coef[1] + 1.96*summary(ASX_triple_vil_base_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_base_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])







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
##  2.1. All data                      ##
##                                     ##
#########################################
#########################################

N_fg <- length( which(ALB_par$pfg_pos==1) )
N_vg <- length( which(ALB_par$pvg_pos==1) )
N_mg <- length( which(ALB_par$pmg_pos==1) )
N_og <- length( which(ALB_par$pog_pos==1) )

N_fvg <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1) )
N_fmg <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1) )
N_fog <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1) )
N_vmg <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1) )
N_vog <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1) )
N_mog <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )

N_fvmg <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1) )
N_fvog <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1) )
N_fmog <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )
N_vmog <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )

N_fvmog <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )


#########################################
##                                     ##
##  2.1.1. All data - double           ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

GAM_double <- rbind( c(N_fg,  N_fvg, N_fmg, N_fog),
                     c(N_fvg, N_vg,  N_vmg, N_vog),
                     c(N_fmg, N_vmg, N_mg,  N_mog),
                     c(N_fog, N_vog, N_mog, N_og ) )

colnames(GAM_double) <- c("Pfg", "Pvg", "Pog", "Pmg")
rownames(GAM_double) <- c("Pfg", "Pvg", "Pog", "Pmg")

GAM_double_obs <- GAM_double/N_data

GAM_double_exp <- (diag(GAM_double)/N_data) %o% (diag(GAM_double)/N_data)
diag(GAM_double_exp) <- diag(GAM_double_obs)

GAM_double_obs_low <- GAM_double
GAM_double_obs_high <- GAM_double

for(i in 1:4)
{
	for(j in 1:4)
	{
		GAM_double_obs_low[i,j]  <- binom.confint( GAM_double[i,j], N_data, method="wilson")[1,5]
		GAM_double_obs_high[i,j] <- binom.confint( GAM_double[i,j], N_data, method="wilson")[1,6]
	}
}


#########################################
## Regression model

x_exp <- c( GAM_double_exp[1,2:4], GAM_double_exp[2,3:4], GAM_double_exp[3,4] )
y_obs <- c( GAM_double_obs[1,2:4], GAM_double_obs[2,3:4], GAM_double_obs[3,4] )

GAM_double_mod <- lm( y_obs ~ 0 + x_exp )




c( summary(GAM_double_mod)$coef[1],
   summary(GAM_double_mod)$coef[1] - 1.96*summary(GAM_double_mod)$coef[2], 
   summary(GAM_double_mod)$coef[1] + 1.96*summary(GAM_double_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_mod)$coef[1], 
         sigma = summary(GAM_double_mod)$coef[2], 
         nu    = summary(GAM_double_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_mod)$coef[1], 
         sigma = summary(GAM_double_mod)$coef[2], 
         nu    = summary(GAM_double_mod)$df[2])




#########################################
##                                     ##
##  2.1.2. All data - triple           ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

GAM_triple <- c(N_vmog, N_fmog, N_fvog, N_fvmg)
names(GAM_triple) <- c("P_vmo", "P_fmo", "P_fvo", "P_fvm")

GAM_triple_obs <- GAM_triple/N_data

GAM_triple_obs_low  <- GAM_triple_obs
GAM_triple_obs_high <- GAM_triple_obs

for(i in 1:4)
{
	GAM_triple_obs_low[i]  <- binom.confint( GAM_triple[i], N_data, method="wilson")[1,5]
	GAM_triple_obs_high[i] <- binom.confint( GAM_triple[i], N_data, method="wilson")[1,6]
}


GAM_triple_exp <- GAM_triple_obs

GAM_triple_exp[1] <- N_vg*N_mg*N_og/( N_data^3 )
GAM_triple_exp[2] <- N_fg*N_mg*N_og/( N_data^3 )
GAM_triple_exp[3] <- N_fg*N_vg*N_mg/( N_data^3 )
GAM_triple_exp[4] <- N_fg*N_vg*N_og/( N_data^3 )



#########################################
## Regression model

GAM_triple_mod <- lm( GAM_triple_obs ~ 0 + GAM_triple_exp )



c( summary(GAM_double_mod)$coef[1],
   summary(GAM_double_mod)$coef[1] - 1.96*summary(GAM_double_mod)$coef[2], 
   summary(GAM_double_mod)$coef[1] + 1.96*summary(GAM_double_mod)$coef[2] )^2

c( summary(GAM_triple_mod)$coef[1],
   summary(GAM_triple_mod)$coef[1] - 1.96*summary(GAM_triple_mod)$coef[2], 
   summary(GAM_triple_mod)$coef[1] + 1.96*summary(GAM_triple_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])

2*pst( q     = summary(GAM_double_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])



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







#########################################
#########################################
##                                     ##
##  2.3. All data                      ##
##       Village & treat stratified    ##
##                                     ##
#########################################
#########################################

N_data_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_mg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_og_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fvg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fmg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vmg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_mog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fvmg_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fvog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_fmog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)
N_vmog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)

N_fvmog_vil_trt <- matrix(NA, nrow=2, ncol=N_village)


for(v in 1:N_village)
{
	N_data_vil_trt[1,v] <- length(which(ALB_par$village == village_names[v] & ALB_par$treat=="PLACEBO"))

	N_fg_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vg_vil_trt[1,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_mg_vil_trt[1,v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_og_vil_trt[1,v] <- length( which(ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fvg_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fmg_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fog_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vmg_vil_trt[1,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vog_vil_trt[1,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_mog_vil_trt[1,v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fvmg_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fvog_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_fmog_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
	N_vmog_vil_trt[1,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )

	N_fvmog_vil_trt[1,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PLACEBO") )
}


for(v in 1:N_village)
{
	N_data_vil_trt[2,v] <- length(which(ALB_par$village == village_names[v] & ALB_par$treat=="PRMQN"))

	N_fg_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vg_vil_trt[2,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_mg_vil_trt[2,v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_og_vil_trt[2,v] <- length( which(ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fvg_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fmg_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fog_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vmg_vil_trt[2,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vog_vil_trt[2,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_mog_vil_trt[2,v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fvmg_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fvog_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_fmog_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
	N_vmog_vil_trt[2,v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )

	N_fvmog_vil_trt[2,v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$treat=="PRMQN") )
}




#########################################
##                                     ##
##  2.3.1. All data - double           ##
##       Village & treat stratified    ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

GAM_double_vil_trt     <- array(NA, dim=c(2, N_village, 4, 4))
GAM_double_obs_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))
GAM_double_exp_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))

GAM_double_obs_low_vil_trt  <- array(NA, dim=c(2, N_village, 4, 4))
GAM_double_obs_high_vil_trt <- array(NA, dim=c(2, N_village, 4, 4))


for(trt in 1:2)
{
	for(v in 1:N_village)
	{
		GAM_double_vil_trt[trt,v,,] <- rbind( c(N_fg_vil_trt[trt,v],  N_fvg_vil_trt[trt,v], N_fmg_vil_trt[trt,v], N_fog_vil_trt[trt,v]),
	      	 	                          c(N_fvg_vil_trt[trt,v], N_vg_vil_trt[trt,v],  N_vmg_vil_trt[trt,v], N_vog_vil_trt[trt,v]),
	            		                    c(N_fmg_vil_trt[trt,v], N_vmg_vil_trt[trt,v], N_mg_vil_trt[trt,v],  N_mog_vil_trt[trt,v]),
	              	                          c(N_fog_vil_trt[trt,v], N_vog_vil_trt[trt,v], N_mog_vil_trt[trt,v], N_og_vil_trt[trt,v] ) )

		colnames(GAM_double_vil_trt[trt,v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")
		rownames(GAM_double_vil_trt[trt,v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")

		GAM_double_obs_vil_trt[trt,v,,] <- GAM_double_vil_trt[trt,v,,]/N_data_vil_trt[trt,v]

		GAM_double_exp_vil_trt[trt,v,,] <- (diag(GAM_double_vil_trt[trt,v,,])/N_data_vil_trt[trt,v]) %o% (diag(GAM_double_vil_trt[trt,v,,])/N_data_vil_trt[trt,v])
		diag(GAM_double_exp_vil_trt[trt,v,,]) <- diag(GAM_double_obs_vil_trt[trt,v,,])


		GAM_double_obs_low_vil_trt[trt,v,,]  <- GAM_double_vil_trt[trt,v,,]
		GAM_double_obs_high_vil_trt[trt,v,,] <- GAM_double_vil_trt[trt,v,,]

		for(i in 1:4)
		{
			for(j in 1:4)
			{
				GAM_double_obs_low_vil_trt[trt,v,i,j]  <- binom.confint( GAM_double_vil_trt[trt,v,i,j], N_data_vil_trt[trt,v], method="wilson")[1,5]
				GAM_double_obs_high_vil_trt[trt,v,i,j] <- binom.confint( GAM_double_vil_trt[trt,v,i,j], N_data_vil_trt[trt,v], method="wilson")[1,6]
			}
		}

	}
}

#########################################
## Regression model

x_exp <- c( GAM_double_exp_vil_trt[,,1,2:4], GAM_double_exp_vil_trt[,,2,3:4], GAM_double_exp_vil_trt[,,3,4] )
y_obs <- c( GAM_double_obs_vil_trt[,,1,2:4], GAM_double_obs_vil_trt[,,2,3:4], GAM_double_obs_vil_trt[,,3,4] )

GAM_double_vil_trt_mod <- lm( y_obs ~ 0 + x_exp )





c( summary(GAM_double_vil_trt_mod)$coef[1],
   summary(GAM_double_vil_trt_mod)$coef[1] - 1.96*summary(GAM_double_vil_trt_mod)$coef[2], 
   summary(GAM_double_vil_trt_mod)$coef[1] + 1.96*summary(GAM_double_vil_trt_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_trt_mod)$coef[1], 
         sigma = summary(GAM_double_vil_trt_mod)$coef[2], 
         nu    = summary(GAM_double_vil_trt_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_trt_mod)$coef[1], 
         sigma = summary(GAM_double_vil_trt_mod)$coef[2], 
         nu    = summary(GAM_double_vil_trt_mod)$df[2])



#########################################
##                                     ##
##  2.3.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

GAM_triple_vil_trt <- c(N_vmog_vil_trt, N_fmog_vil_trt, N_fvog_vil_trt, N_fvmg_vil_trt)
names(GAM_triple_vil_trt) <- c( rep("P_vmog", 2*N_village), rep("P_fmog", 2*N_village), rep("P_fvog", 2*N_village), rep("P_fvmg", 2*N_village) )


GAM_triple_vil_trt_obs      <- GAM_triple_vil_trt
GAM_triple_vil_trt_obs_low  <- GAM_triple_vil_trt
GAM_triple_vil_trt_obs_high <- GAM_triple_vil_trt


for(v in 1:N_village)
{
	GAM_triple_vil_trt_obs[0*N_village + v] <- GAM_triple_vil_trt[0*N_village + v]/N_data_vil_trt[1,v]
	GAM_triple_vil_trt_obs[1*N_village + v] <- GAM_triple_vil_trt[1*N_village + v]/N_data_vil_trt[2,v]

	GAM_triple_vil_trt_obs[2*N_village + v] <- GAM_triple_vil_trt[2*N_village + v]/N_data_vil_trt[1,v]
	GAM_triple_vil_trt_obs[3*N_village + v] <- GAM_triple_vil_trt[3*N_village + v]/N_data_vil_trt[2,v]

	GAM_triple_vil_trt_obs[4*N_village + v] <- GAM_triple_vil_trt[4*N_village + v]/N_data_vil_trt[1,v]
	GAM_triple_vil_trt_obs[5*N_village + v] <- GAM_triple_vil_trt[5*N_village + v]/N_data_vil_trt[2,v]

	GAM_triple_vil_trt_obs[6*N_village + v] <- GAM_triple_vil_trt[6*N_village + v]/N_data_vil_trt[1,v]
	GAM_triple_vil_trt_obs[7*N_village + v] <- GAM_triple_vil_trt[7*N_village + v]/N_data_vil_trt[2,v]
}

for(v in 1:N_village)
{
	GAM_triple_vil_trt_obs_low[0*N_village + v] <- binom.confint( GAM_triple_vil_trt[0*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	GAM_triple_vil_trt_obs_low[1*N_village + v] <- binom.confint( GAM_triple_vil_trt[1*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	GAM_triple_vil_trt_obs_low[2*N_village + v] <- binom.confint( GAM_triple_vil_trt[2*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	GAM_triple_vil_trt_obs_low[3*N_village + v] <- binom.confint( GAM_triple_vil_trt[3*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	GAM_triple_vil_trt_obs_low[4*N_village + v] <- binom.confint( GAM_triple_vil_trt[4*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	GAM_triple_vil_trt_obs_low[5*N_village + v] <- binom.confint( GAM_triple_vil_trt[5*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]

	GAM_triple_vil_trt_obs_low[6*N_village + v] <- binom.confint( GAM_triple_vil_trt[6*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,5]
	GAM_triple_vil_trt_obs_low[7*N_village + v] <- binom.confint( GAM_triple_vil_trt[7*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,5]


	GAM_triple_vil_trt_obs_high[0*N_village + v] <- binom.confint( GAM_triple_vil_trt[0*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	GAM_triple_vil_trt_obs_high[1*N_village + v] <- binom.confint( GAM_triple_vil_trt[1*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	GAM_triple_vil_trt_obs_high[2*N_village + v] <- binom.confint( GAM_triple_vil_trt[2*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	GAM_triple_vil_trt_obs_high[3*N_village + v] <- binom.confint( GAM_triple_vil_trt[3*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	GAM_triple_vil_trt_obs_high[4*N_village + v] <- binom.confint( GAM_triple_vil_trt[4*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	GAM_triple_vil_trt_obs_high[5*N_village + v] <- binom.confint( GAM_triple_vil_trt[5*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]

	GAM_triple_vil_trt_obs_high[6*N_village + v] <- binom.confint( GAM_triple_vil_trt[6*N_village + v], N_data_vil_trt[1,v], method="wilson")[1,6]
	GAM_triple_vil_trt_obs_high[7*N_village + v] <- binom.confint( GAM_triple_vil_trt[7*N_village + v], N_data_vil_trt[2,v], method="wilson")[1,6]
}





GAM_triple_vil_trt_exp <- GAM_triple_vil_trt_obs

for(v in 1:N_village)
{
	GAM_triple_vil_trt_exp[0*N_village + v] <- N_vg_vil_trt[1,v]*N_mg_vil_trt[1,v]*N_og_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	GAM_triple_vil_trt_exp[1*N_village + v] <- N_vg_vil_trt[2,v]*N_mg_vil_trt[2,v]*N_og_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	GAM_triple_vil_trt_exp[2*N_village + v] <- N_fg_vil_trt[1,v]*N_mg_vil_trt[1,v]*N_og_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	GAM_triple_vil_trt_exp[3*N_village + v] <- N_fg_vil_trt[2,v]*N_mg_vil_trt[2,v]*N_og_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	GAM_triple_vil_trt_exp[4*N_village + v] <- N_fg_vil_trt[1,v]*N_vg_vil_trt[1,v]*N_og_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	GAM_triple_vil_trt_exp[5*N_village + v] <- N_fg_vil_trt[2,v]*N_vg_vil_trt[2,v]*N_og_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )

	GAM_triple_vil_trt_exp[6*N_village + v] <- N_fg_vil_trt[1,v]*N_vg_vil_trt[1,v]*N_mg_vil_trt[1,v]/( N_data_vil_trt[1,v]^3 )
	GAM_triple_vil_trt_exp[7*N_village + v] <- N_fg_vil_trt[2,v]*N_vg_vil_trt[2,v]*N_mg_vil_trt[2,v]/( N_data_vil_trt[2,v]^3 )
}


#########################################
## Regression model

GAM_triple_vil_trt_mod <- lm( GAM_triple_vil_trt_obs ~ 0 + GAM_triple_vil_trt_exp )



c( summary(GAM_double_vil_trt_mod)$coef[1],
   summary(GAM_double_vil_trt_mod)$coef[1] - 1.96*summary(GAM_double_vil_trt_mod)$coef[2], 
   summary(GAM_double_vil_trt_mod)$coef[1] + 1.96*summary(GAM_double_vil_trt_mod)$coef[2] )^2

c( summary(GAM_triple_vil_trt_mod)$coef[1],
   summary(GAM_triple_vil_trt_mod)$coef[1] - 1.96*summary(GAM_triple_vil_trt_mod)$coef[2], 
   summary(GAM_triple_vil_trt_mod)$coef[1] + 1.96*summary(GAM_triple_vil_trt_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_trt_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])






#########################################
#########################################
##                                     ##
##  2.4. All data                      ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################
#########################################

#########################################
##                                     ##
##  2.4.1. All data - double           ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################

#########################################
## Regression model

x_exp <- c( GAM_double_exp_vil_trt[1,,1,2:4], GAM_double_exp_vil_trt[1,,2,3:4], GAM_double_exp_vil_trt[1,,3,4] )
y_obs <- c( GAM_double_obs_vil_trt[1,,1,2:4], GAM_double_obs_vil_trt[1,,2,3:4], GAM_double_obs_vil_trt[1,,3,4] )

GAM_double_vil_plac_mod <- lm( y_obs ~ 0 + x_exp )




c( summary(GAM_double_vil_plac_mod)$coef[1],
   summary(GAM_double_vil_plac_mod)$coef[1] - 1.96*summary(GAM_double_vil_plac_mod)$coef[2], 
   summary(GAM_double_vil_plac_mod)$coef[1] + 1.96*summary(GAM_double_vil_plac_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_plac_mod)$coef[1], 
         sigma = summary(GAM_double_vil_plac_mod)$coef[2], 
         nu    = summary(GAM_double_vil_plac_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_plac_mod)$coef[1], 
         sigma = summary(GAM_double_vil_plac_mod)$coef[2], 
         nu    = summary(GAM_double_vil_plac_mod)$df[2])



#########################################
##                                     ##
##  2.4.2. All data - triple           ##
##       Village & treat stratified    ##
##       Placebo only                  ##
##                                     ##
#########################################

#########################################
## Regression model

GAM_triple_vil_plac_obs <- GAM_triple_vil_trt_obs[c(1:5,11:15,21:25,31:35)]
GAM_triple_vil_plac_exp <- GAM_triple_vil_trt_exp[c(1:5,11:15,21:25,31:35)]


GAM_triple_vil_plac_mod <- lm( GAM_triple_vil_plac_obs ~ 0 + GAM_triple_vil_plac_exp )




c( summary(GAM_double_vil_plac_mod)$coef[1],
   summary(GAM_double_vil_plac_mod)$coef[1] - 1.96*summary(GAM_double_vil_plac_mod)$coef[2], 
   summary(GAM_double_vil_plac_mod)$coef[1] + 1.96*summary(GAM_double_vil_plac_mod)$coef[2] )^2

c( summary(GAM_triple_vil_plac_mod)$coef[1],
   summary(GAM_triple_vil_plac_mod)$coef[1] - 1.96*summary(GAM_triple_vil_plac_mod)$coef[2], 
   summary(GAM_triple_vil_plac_mod)$coef[1] + 1.96*summary(GAM_triple_vil_plac_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_plac_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])






#########################################
#########################################
##                                     ##
##  2.5. Baseline data                 ##
##                                     ##
#########################################
#########################################

N_data_base <- length( which(ALB_par$study_day<0) )

N_fg_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$study_day<0) )
N_vg_base <- length( which(ALB_par$pvg_pos==1 & ALB_par$study_day<0) )
N_mg_base <- length( which(ALB_par$pmg_pos==1 & ALB_par$study_day<0) )
N_og_base <- length( which(ALB_par$pog_pos==1 & ALB_par$study_day<0) )

N_fvg_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$study_day<0) )
N_fmg_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$study_day<0) )
N_fog_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )
N_vmg_base <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$study_day<0) )
N_vog_base <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )
N_mog_base <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )

N_fvmg_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$study_day<0) )
N_fvog_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )
N_fmog_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )
N_vmog_base <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )

N_fvmog_base <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$study_day<0) )


#########################################
##                                     ##
##  2.5.1. All data - double           ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

GAM_double_base <- rbind( c(N_fg_base,  N_fvg_base, N_fmg_base, N_fog_base),
                          c(N_fvg_base, N_vg_base,  N_vmg_base, N_vog_base),
                          c(N_fmg_base, N_vmg_base, N_mg_base,  N_mog_base),
                          c(N_fog_base, N_vog_base, N_mog_base, N_og_base ) )

colnames(GAM_double_base) <- c("Pfg", "Pvg", "Pog", "Pmg")
rownames(GAM_double_base) <- c("Pfg", "Pvg", "Pog", "Pmg")

GAM_double_base_obs <- GAM_double_base/N_data_base

GAM_double_base_exp <- (diag(GAM_double_base)/N_data_base) %o% (diag(GAM_double_base)/N_data_base)
diag(GAM_double_base_exp) <- diag(GAM_double_base_obs)

GAM_double_base_obs_low <- GAM_double_base
GAM_double_base_obs_high <- GAM_double_base

for(i in 1:4)
{
	for(j in 1:4)
	{
		GAM_double_base_obs_low[i,j]  <- binom.confint( GAM_double_base[i,j], N_data_base, method="wilson")[1,5]
		GAM_double_base_obs_high[i,j] <- binom.confint( GAM_double_base[i,j], N_data_base, method="wilson")[1,6]
	}
}


#########################################
## Regression model

x_exp <- c( GAM_double_base_exp[1,2:4], GAM_double_base_exp[2,3:4], GAM_double_base_exp[3,4] )
y_obs <- c( GAM_double_base_obs[1,2:4], GAM_double_base_obs[2,3:4], GAM_double_base_obs[3,4] )

GAM_double_base_mod <- lm( y_obs ~ 0 + x_exp )




c( summary(GAM_double_base_mod)$coef[1],
   summary(GAM_double_base_mod)$coef[1] - 1.96*summary(GAM_double_base_mod)$coef[2], 
   summary(GAM_double_base_mod)$coef[1] + 1.96*summary(GAM_double_base_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_base_mod)$coef[1], 
         sigma = summary(GAM_double_base_mod)$coef[2], 
         nu    = summary(GAM_double_base_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_base_mod)$coef[1], 
         sigma = summary(GAM_double_base_mod)$coef[2], 
         nu    = summary(GAM_double_base_mod)$df[2])




#########################################
##                                     ##
##  2.5.2. All data - triple           ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

GAM_triple_base <- c(N_vmog_base, N_fmog_base, N_fvog_base, N_fvmg_base)
names(GAM_triple_base) <- c("P_vmog", "P_fmog", "P_fvog", "P_fvmg")

GAM_triple_base_obs <- GAM_triple_base/N_data_base

GAM_triple_base_obs_low  <- GAM_triple_base_obs
GAM_triple_base_obs_high <- GAM_triple_base_obs

for(i in 1:4)
{
	GAM_triple_base_obs_low[i]  <- binom.confint( GAM_triple_base[i], N_data_base, method="wilson")[1,5]
	GAM_triple_base_obs_high[i] <- binom.confint( GAM_triple_base[i], N_data_base, method="wilson")[1,6]
}


GAM_triple_base_exp <- GAM_triple_base_obs

GAM_triple_base_exp[1] <- N_vg_base*N_mg_base*N_og_base/( N_data_base^3 )
GAM_triple_base_exp[2] <- N_fg_base*N_mg_base*N_og_base/( N_data_base^3 )
GAM_triple_base_exp[3] <- N_fg_base*N_vg_base*N_mg_base/( N_data_base^3 )
GAM_triple_base_exp[4] <- N_fg_base*N_vg_base*N_og_base/( N_data_base^3 )



#########################################
## Regression model

GAM_triple_base_mod <- lm( GAM_triple_base_obs ~ 0 + GAM_triple_base_exp )




c( summary(GAM_double_base_mod)$coef[1],
   summary(GAM_double_base_mod)$coef[1] - 1.96*summary(GAM_double_base_mod)$coef[2], 
   summary(GAM_double_base_mod)$coef[1] + 1.96*summary(GAM_double_base_mod)$coef[2] )^2

c( summary(GAM_triple_base_mod)$coef[1],
   summary(GAM_triple_base_mod)$coef[1] - 1.96*summary(GAM_triple_base_mod)$coef[2], 
   summary(GAM_triple_base_mod)$coef[1] + 1.96*summary(GAM_triple_base_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])

2*pst( q     = summary(GAM_double_base_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])






#########################################
#########################################
##                                     ##
##  2.6. Baseline data                 ##
##       Village stratified            ##
##                                     ##
#########################################
#########################################

N_data_vil_base <- rep(NA, N_village)

N_fg_vil_base <- rep(NA, N_village)
N_vg_vil_base <- rep(NA, N_village)
N_mg_vil_base <- rep(NA, N_village)
N_og_vil_base <- rep(NA, N_village)

N_fvg_vil_base <- rep(NA, N_village)
N_fmg_vil_base <- rep(NA, N_village)
N_fog_vil_base <- rep(NA, N_village)
N_vmg_vil_base <- rep(NA, N_village)
N_vog_vil_base <- rep(NA, N_village)
N_mog_vil_base <- rep(NA, N_village)

N_fvmg_vil_base <- rep(NA, N_village)
N_fvog_vil_base <- rep(NA, N_village)
N_fmog_vil_base <- rep(NA, N_village)
N_vmog_vil_base <- rep(NA, N_village)

N_fvmog_vil_base <- rep(NA, N_village)

for(v in 1:N_village)
{
	N_data_vil_base[v] <- length(which(ALB_par$village == village_names[v] & ALB_par$study_day<0))

	N_fg_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vg_vil_base[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_mg_vil_base[v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_og_vil_base[v] <- length( which(ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fvg_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fmg_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fog_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vmg_vil_base[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vog_vil_base[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_mog_vil_base[v] <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fvmg_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fvog_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_fmog_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
	N_vmog_vil_base[v] <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )

	N_fvmog_vil_base[v] <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1 & ALB_par$village==village_names[v] & ALB_par$study_day<0) )
}



#########################################
##                                     ##
##  2.6.1. All data - double           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Double infection data                                     

GAM_double_vil_base <- array(NA, dim=c(N_village, 4, 4))
GAM_double_obs_vil_base <- array(NA, dim=c(N_village, 4, 4))
GAM_double_exp_vil_base <- array(NA, dim=c(N_village, 4, 4))

GAM_double_obs_low_vil_base <- array(NA, dim=c(N_village, 4, 4))
GAM_double_obs_high_vil_base <- array(NA, dim=c(N_village, 4, 4))


for(v in 1:N_village)
{
	GAM_double_vil_base[v,,] <- rbind( c(N_fg_vil_base[v],  N_fvg_vil_base[v], N_fmg_vil_base[v], N_fog_vil_base[v]),
      	 	                       c(N_fvg_vil_base[v], N_vg_vil_base[v],  N_vmg_vil_base[v], N_vog_vil_base[v]),
            		                 c(N_fmg_vil_base[v], N_vmg_vil_base[v], N_mg_vil_base[v],  N_mog_vil_base[v]),
              	                       c(N_fog_vil_base[v], N_vog_vil_base[v], N_mog_vil_base[v], N_og_vil_base[v] ) )

	colnames(GAM_double_vil_base[v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")
	rownames(GAM_double_vil_base[v,,]) <- c("Pfg", "Pvg", "Pog", "Pmg")

	GAM_double_obs_vil_base[v,,] <- GAM_double_vil_base[v,,]/N_data_vil_base[v]

	GAM_double_exp_vil_base[v,,] <- (diag(GAM_double_vil_base[v,,])/N_data_vil_base[v]) %o% (diag(GAM_double_vil_base[v,,])/N_data_vil_base[v])
	diag(GAM_double_exp_vil_base[v,,]) <- diag(GAM_double_obs_vil_base[v,,])


	GAM_double_obs_low_vil_base[v,,]  <- GAM_double_vil_base[v,,]
	GAM_double_obs_high_vil_base[v,,] <- GAM_double_vil_base[v,,]

	for(i in 1:4)
	{
		for(j in 1:4)
		{
			GAM_double_obs_low_vil_base[v,i,j]  <- binom.confint( GAM_double_vil_base[v,i,j], N_data_vil_base[v], method="wilson")[1,5]
			GAM_double_obs_high_vil_base[v,i,j] <- binom.confint( GAM_double_vil_base[v,i,j], N_data_vil_base[v], method="wilson")[1,6]
		}
	}

}


#########################################
## Regression model

x_exp <- c( GAM_double_exp_vil_base[,1,2:4], GAM_double_exp_vil_base[,2,3:4], GAM_double_exp_vil_base[,3,4] )
y_obs <- c( GAM_double_obs_vil_base[,1,2:4], GAM_double_obs_vil_base[,2,3:4], GAM_double_obs_vil_base[,3,4] )

GAM_double_vil_base_mod <- lm( y_obs ~ 0 + x_exp )



c( summary(GAM_double_vil_base_mod)$coef[1],
   summary(GAM_double_vil_base_mod)$coef[1] - 1.96*summary(GAM_double_vil_base_mod)$coef[2], 
   summary(GAM_double_vil_base_mod)$coef[1] + 1.96*summary(GAM_double_vil_base_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_base_mod)$coef[1], 
         sigma = summary(GAM_double_vil_base_mod)$coef[2], 
         nu    = summary(GAM_double_vil_base_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_base_mod)$coef[1], 
         sigma = summary(GAM_double_vil_base_mod)$coef[2], 
         nu    = summary(GAM_double_vil_base_mod)$df[2])



#########################################
##                                     ##
##  2.6.2. All data - triple           ##
##       Village stratified            ##
##                                     ##
#########################################

#########################################
## Triple infection data                                     

GAM_triple_vil_base <- c(N_vmog_vil_base, N_fmog_vil_base, N_fvog_vil_base, N_fvmg_vil_base)
names(GAM_triple_vil_base) <- c( rep("P_vmog", N_village), rep("P_fmog", N_village), rep("P_fvog", N_village), rep("P_fvmg", N_village) )

GAM_triple_vil_base_obs      <- GAM_triple_vil_base
GAM_triple_vil_base_obs_low  <- GAM_triple_vil_base
GAM_triple_vil_base_obs_high <- GAM_triple_vil_base

for(v in 1:N_village)
{
	GAM_triple_vil_base_obs[0*N_village + v] <- GAM_triple_vil_base[0*N_village + v]/N_data_vil_base[v]
	GAM_triple_vil_base_obs[1*N_village + v] <- GAM_triple_vil_base[1*N_village + v]/N_data_vil_base[v]
	GAM_triple_vil_base_obs[2*N_village + v] <- GAM_triple_vil_base[2*N_village + v]/N_data_vil_base[v]
	GAM_triple_vil_base_obs[3*N_village + v] <- GAM_triple_vil_base[3*N_village + v]/N_data_vil_base[v]
}



for(v in 1:N_village)
{
	GAM_triple_vil_base_obs_low[0*N_village + v] <- binom.confint( GAM_triple_vil_base[0*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	GAM_triple_vil_base_obs_low[1*N_village + v] <- binom.confint( GAM_triple_vil_base[1*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	GAM_triple_vil_base_obs_low[2*N_village + v] <- binom.confint( GAM_triple_vil_base[2*N_village + v], N_data_vil_base[v], method="wilson")[1,5]
	GAM_triple_vil_base_obs_low[3*N_village + v] <- binom.confint( GAM_triple_vil_base[3*N_village + v], N_data_vil_base[v], method="wilson")[1,5]

	GAM_triple_vil_base_obs_high[0*N_village + v] <- binom.confint( GAM_triple_vil_base[0*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	GAM_triple_vil_base_obs_high[1*N_village + v] <- binom.confint( GAM_triple_vil_base[1*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	GAM_triple_vil_base_obs_high[2*N_village + v] <- binom.confint( GAM_triple_vil_base[2*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
	GAM_triple_vil_base_obs_high[3*N_village + v] <- binom.confint( GAM_triple_vil_base[3*N_village + v], N_data_vil_base[v], method="wilson")[1,6]
}





GAM_triple_vil_base_exp <- GAM_triple_vil_base_obs

for(v in 1:N_village)
{
	GAM_triple_vil_base_exp[0*N_village + v] <- N_vg_vil_base[v]*N_mg_vil_base[v]*N_og_vil_base[v]/( N_data_vil_base[v]^3 )
	GAM_triple_vil_base_exp[1*N_village + v] <- N_fg_vil_base[v]*N_mg_vil_base[v]*N_og_vil_base[v]/( N_data_vil_base[v]^3 )
	GAM_triple_vil_base_exp[2*N_village + v] <- N_fg_vil_base[v]*N_vg_vil_base[v]*N_og_vil_base[v]/( N_data_vil_base[v]^3 )
	GAM_triple_vil_base_exp[3*N_village + v] <- N_fg_vil_base[v]*N_vg_vil_base[v]*N_mg_vil_base[v]/( N_data_vil_base[v]^3 )
}


#########################################
## Regression model

GAM_triple_vil_base_mod <- lm( GAM_triple_vil_base_obs ~ 0 + GAM_triple_vil_base_exp )


c( summary(GAM_double_vil_base_mod)$coef[1],
   summary(GAM_double_vil_base_mod)$coef[1] - 1.96*summary(GAM_double_vil_base_mod)$coef[2], 
   summary(GAM_double_vil_base_mod)$coef[1] + 1.96*summary(GAM_double_vil_base_mod)$coef[2] )^2

c( summary(GAM_triple_vil_base_mod)$coef[1],
   summary(GAM_triple_vil_base_mod)$coef[1] - 1.96*summary(GAM_triple_vil_base_mod)$coef[2], 
   summary(GAM_triple_vil_base_mod)$coef[1] + 1.96*summary(GAM_triple_vil_base_mod)$coef[2] )

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_base_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])



########################################
########################################
##                                    ##
##  ######  ####  #####  ##    #####  ##
##    ##   ##  ## ##  ## ##    ##     ##
##    ##   ###### #####  ##    ####   ##
##    ##   ##  ## ##  ## ##    ##     ##
##    ##   ##  ## #####  ##### #####  ##
##                                    ##
########################################
########################################

########################################
##                                    ##
##  Asexual parasites                 ##
##                                    ##
########################################

########################################
## All samples; double                                   

c( summary(ASX_double_mod)$coef[1],
   summary(ASX_double_mod)$coef[1] - 1.96*summary(ASX_double_mod)$coef[2], 
   summary(ASX_double_mod)$coef[1] + 1.96*summary(ASX_double_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_mod)$coef[1], 
         sigma = summary(ASX_double_mod)$coef[2], 
         nu    = summary(ASX_double_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_mod)$coef[1], 
         sigma = summary(ASX_double_mod)$coef[2], 
         nu    = summary(ASX_double_mod)$df[2])


########################################
## All samples; triple

c( summary(ASX_triple_mod)$coef[1],
   summary(ASX_triple_mod)$coef[1] - 1.96*summary(ASX_triple_mod)$coef[2], 
   summary(ASX_triple_mod)$coef[1] + 1.96*summary(ASX_triple_mod)$coef[2] )


c( summary(ASX_triple_mod)$coef[1],
   summary(ASX_triple_mod)$coef[1] - 1.96*summary(ASX_triple_mod)$coef[2], 
   summary(ASX_triple_mod)$coef[1] + 1.96*summary(ASX_triple_mod)$coef[2] )/summary(ASX_double_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])

2*pst( q     = summary(ASX_double_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_mod)$coef[1], 
       sigma = summary(ASX_triple_mod)$coef[2], 
       nu    = summary(ASX_triple_mod)$df[2])



########################################
## All samples (village strat); double                                   

c( summary(ASX_double_vil_mod)$coef[1],
   summary(ASX_double_vil_mod)$coef[1] - 1.96*summary(ASX_double_vil_mod)$coef[2], 
   summary(ASX_double_vil_mod)$coef[1] + 1.96*summary(ASX_double_vil_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_mod)$coef[1], 
         sigma = summary(ASX_double_vil_mod)$coef[2], 
         nu    = summary(ASX_double_vil_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_mod)$coef[1], 
         sigma = summary(ASX_double_vil_mod)$coef[2], 
         nu    = summary(ASX_double_vil_mod)$df[2])


########################################
## All samples (village strat); triple


c( summary(ASX_triple_vil_mod)$coef[1],
   summary(ASX_triple_vil_mod)$coef[1] - 1.96*summary(ASX_triple_vil_mod)$coef[2], 
   summary(ASX_triple_vil_mod)$coef[1] + 1.96*summary(ASX_triple_vil_mod)$coef[2] )

c( summary(ASX_triple_vil_mod)$coef[1],
   summary(ASX_triple_vil_mod)$coef[1] - 1.96*summary(ASX_triple_vil_mod)$coef[2], 
   summary(ASX_triple_vil_mod)$coef[1] + 1.96*summary(ASX_triple_vil_mod)$coef[2] )/summary(ASX_double_vil_mod)$coef[1]^2


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



########################################
## All samples (village & treat strat); double                                   

c( summary(ASX_double_vil_trt_mod)$coef[1],
   summary(ASX_double_vil_trt_mod)$coef[1] - 1.96*summary(ASX_double_vil_trt_mod)$coef[2], 
   summary(ASX_double_vil_trt_mod)$coef[1] + 1.96*summary(ASX_double_vil_trt_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_trt_mod)$coef[1], 
         sigma = summary(ASX_double_vil_trt_mod)$coef[2], 
         nu    = summary(ASX_double_vil_trt_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_trt_mod)$coef[1], 
         sigma = summary(ASX_double_vil_trt_mod)$coef[2], 
         nu    = summary(ASX_double_vil_trt_mod)$df[2])


########################################
## All samples (village & treat strat); triple

c( summary(ASX_triple_vil_trt_mod)$coef[1],
   summary(ASX_triple_vil_trt_mod)$coef[1] - 1.96*summary(ASX_triple_vil_trt_mod)$coef[2], 
   summary(ASX_triple_vil_trt_mod)$coef[1] + 1.96*summary(ASX_triple_vil_trt_mod)$coef[2] )

c( summary(ASX_triple_vil_trt_mod)$coef[1],
   summary(ASX_triple_vil_trt_mod)$coef[1] - 1.96*summary(ASX_triple_vil_trt_mod)$coef[2], 
   summary(ASX_triple_vil_trt_mod)$coef[1] + 1.96*summary(ASX_triple_vil_trt_mod)$coef[2] )/summary(ASX_double_vil_trt_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_trt_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_trt_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_trt_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_trt_mod)$df[2])


########################################
## All samples (village strat, placeboe); double                                   

c( summary(ASX_double_vil_plac_mod)$coef[1],
   summary(ASX_double_vil_plac_mod)$coef[1] - 1.96*summary(ASX_double_vil_plac_mod)$coef[2], 
   summary(ASX_double_vil_plac_mod)$coef[1] + 1.96*summary(ASX_double_vil_plac_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_plac_mod)$coef[1], 
         sigma = summary(ASX_double_vil_plac_mod)$coef[2], 
         nu    = summary(ASX_double_vil_plac_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_plac_mod)$coef[1], 
         sigma = summary(ASX_double_vil_plac_mod)$coef[2], 
         nu    = summary(ASX_double_vil_plac_mod)$df[2])


########################################
## All samples (village strat, placebo); triple

c( summary(ASX_triple_vil_plac_mod)$coef[1],
   summary(ASX_triple_vil_plac_mod)$coef[1] - 1.96*summary(ASX_triple_vil_plac_mod)$coef[2], 
   summary(ASX_triple_vil_plac_mod)$coef[1] + 1.96*summary(ASX_triple_vil_plac_mod)$coef[2] )

c( summary(ASX_triple_vil_plac_mod)$coef[1],
   summary(ASX_triple_vil_plac_mod)$coef[1] - 1.96*summary(ASX_triple_vil_plac_mod)$coef[2], 
   summary(ASX_triple_vil_plac_mod)$coef[1] + 1.96*summary(ASX_triple_vil_plac_mod)$coef[2] )/summary(ASX_double_vil_plac_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_plac_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_plac_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_plac_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_plac_mod)$df[2])




########################################
## Baseline samples; double                                   

c( summary(ASX_double_base_mod)$coef[1],
   summary(ASX_double_base_mod)$coef[1] - 1.96*summary(ASX_double_base_mod)$coef[2], 
   summary(ASX_double_base_mod)$coef[1] + 1.96*summary(ASX_double_base_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_base_mod)$coef[1], 
         sigma = summary(ASX_double_base_mod)$coef[2], 
         nu    = summary(ASX_double_base_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_base_mod)$coef[1], 
         sigma = summary(ASX_double_base_mod)$coef[2], 
         nu    = summary(ASX_double_base_mod)$df[2])


########################################
## Baseline samples; triple

c( summary(ASX_triple_base_mod)$coef[1],
   summary(ASX_triple_base_mod)$coef[1] - 1.96*summary(ASX_triple_base_mod)$coef[2], 
   summary(ASX_triple_base_mod)$coef[1] + 1.96*summary(ASX_triple_base_mod)$coef[2] )

c( summary(ASX_triple_base_mod)$coef[1],
   summary(ASX_triple_base_mod)$coef[1] - 1.96*summary(ASX_triple_base_mod)$coef[2], 
   summary(ASX_triple_base_mod)$coef[1] + 1.96*summary(ASX_triple_base_mod)$coef[2] )/summary(ASX_double_base_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])

2*pst( q     = summary(ASX_double_base_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_base_mod)$coef[1], 
       sigma = summary(ASX_triple_base_mod)$coef[2], 
       nu    = summary(ASX_triple_base_mod)$df[2])


########################################
## Baseline samples (village strat); double                                   

c( summary(ASX_double_vil_base_mod)$coef[1],
   summary(ASX_double_vil_base_mod)$coef[1] - 1.96*summary(ASX_double_vil_base_mod)$coef[2], 
   summary(ASX_double_vil_base_mod)$coef[1] + 1.96*summary(ASX_double_vil_base_mod)$coef[2] )

2*pst(0, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])

2*pst(1, mu    = summary(ASX_double_vil_base_mod)$coef[1], 
         sigma = summary(ASX_double_vil_base_mod)$coef[2], 
         nu    = summary(ASX_double_vil_base_mod)$df[2])


########################################
## Baseline samples (village strat); triple


c( summary(ASX_triple_vil_base_mod)$coef[1],
   summary(ASX_triple_vil_base_mod)$coef[1] - 1.96*summary(ASX_triple_vil_base_mod)$coef[2], 
   summary(ASX_triple_vil_base_mod)$coef[1] + 1.96*summary(ASX_triple_vil_base_mod)$coef[2] )

c( summary(ASX_triple_vil_base_mod)$coef[1],
   summary(ASX_triple_vil_base_mod)$coef[1] - 1.96*summary(ASX_triple_vil_base_mod)$coef[2], 
   summary(ASX_triple_vil_base_mod)$coef[1] + 1.96*summary(ASX_triple_vil_base_mod)$coef[2] )/summary(ASX_double_vil_base_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])

2*pst( q     = summary(ASX_double_vil_base_mod)$coef[1]^2, 
       mu    = summary(ASX_triple_vil_base_mod)$coef[1], 
       sigma = summary(ASX_triple_vil_base_mod)$coef[2], 
       nu    = summary(ASX_triple_vil_base_mod)$df[2])



########################################
##                                    ##
##  Gametocytes                       ##
##                                    ##
########################################

########################################
## All samples; double                                   

c( summary(GAM_double_mod)$coef[1],
   summary(GAM_double_mod)$coef[1] - 1.96*summary(GAM_double_mod)$coef[2], 
   summary(GAM_double_mod)$coef[1] + 1.96*summary(GAM_double_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_mod)$coef[1], 
         sigma = summary(GAM_double_mod)$coef[2], 
         nu    = summary(GAM_double_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_mod)$coef[1], 
         sigma = summary(GAM_double_mod)$coef[2], 
         nu    = summary(GAM_double_mod)$df[2])




########################################
## All samples; triple

c( summary(GAM_triple_mod)$coef[1],
   summary(GAM_triple_mod)$coef[1] - 1.96*summary(GAM_triple_mod)$coef[2], 
   summary(GAM_triple_mod)$coef[1] + 1.96*summary(GAM_triple_mod)$coef[2] )

c( summary(GAM_triple_mod)$coef[1],
   summary(GAM_triple_mod)$coef[1] - 1.96*summary(GAM_triple_mod)$coef[2], 
   summary(GAM_triple_mod)$coef[1] + 1.96*summary(GAM_triple_mod)$coef[2] )/summary(GAM_double_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])

2*pst( q     = summary(GAM_double_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_mod)$coef[1], 
       sigma = summary(GAM_triple_mod)$coef[2], 
       nu    = summary(GAM_triple_mod)$df[2])



########################################
## All samples (village strat); double                                   

c( summary(GAM_double_vil_mod)$coef[1],
   summary(GAM_double_vil_mod)$coef[1] - 1.96*summary(GAM_double_vil_mod)$coef[2], 
   summary(GAM_double_vil_mod)$coef[1] + 1.96*summary(GAM_double_vil_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_mod)$coef[1], 
         sigma = summary(GAM_double_vil_mod)$coef[2], 
         nu    = summary(GAM_double_vil_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_mod)$coef[1], 
         sigma = summary(GAM_double_vil_mod)$coef[2], 
         nu    = summary(GAM_double_vil_mod)$df[2])


########################################
## All samples (village strat); triple

c( summary(GAM_triple_vil_mod)$coef[1],
   summary(GAM_triple_vil_mod)$coef[1] - 1.96*summary(GAM_triple_vil_mod)$coef[2], 
   summary(GAM_triple_vil_mod)$coef[1] + 1.96*summary(GAM_triple_vil_mod)$coef[2] )

c( summary(GAM_triple_vil_mod)$coef[1],
   summary(GAM_triple_vil_mod)$coef[1] - 1.96*summary(GAM_triple_vil_mod)$coef[2], 
   summary(GAM_triple_vil_mod)$coef[1] + 1.96*summary(GAM_triple_vil_mod)$coef[2] )/summary(GAM_double_vil_mod)$coef[1]^2

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



########################################
## All samples (village & treat strat); double                                   

c( summary(GAM_double_vil_trt_mod)$coef[1],
   summary(GAM_double_vil_trt_mod)$coef[1] - 1.96*summary(GAM_double_vil_trt_mod)$coef[2], 
   summary(GAM_double_vil_trt_mod)$coef[1] + 1.96*summary(GAM_double_vil_trt_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_trt_mod)$coef[1], 
         sigma = summary(GAM_double_vil_trt_mod)$coef[2], 
         nu    = summary(GAM_double_vil_trt_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_trt_mod)$coef[1], 
         sigma = summary(GAM_double_vil_trt_mod)$coef[2], 
         nu    = summary(GAM_double_vil_trt_mod)$df[2])


########################################
## All samples (village & treat strat); triple


c( summary(GAM_triple_vil_trt_mod)$coef[1],
   summary(GAM_triple_vil_trt_mod)$coef[1] - 1.96*summary(GAM_triple_vil_trt_mod)$coef[2], 
   summary(GAM_triple_vil_trt_mod)$coef[1] + 1.96*summary(GAM_triple_vil_trt_mod)$coef[2] )

c( summary(GAM_triple_vil_trt_mod)$coef[1],
   summary(GAM_triple_vil_trt_mod)$coef[1] - 1.96*summary(GAM_triple_vil_trt_mod)$coef[2], 
   summary(GAM_triple_vil_trt_mod)$coef[1] + 1.96*summary(GAM_triple_vil_trt_mod)$coef[2] )/summary(GAM_double_vil_trt_mod)$coef[1]^2


2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_trt_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_trt_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_trt_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_trt_mod)$df[2])


########################################
## All samples (village strat, placeboe); double                                   

c( summary(GAM_double_vil_plac_mod)$coef[1],
   summary(GAM_double_vil_plac_mod)$coef[1] - 1.96*summary(GAM_double_vil_plac_mod)$coef[2], 
   summary(GAM_double_vil_plac_mod)$coef[1] + 1.96*summary(GAM_double_vil_plac_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_plac_mod)$coef[1], 
         sigma = summary(GAM_double_vil_plac_mod)$coef[2], 
         nu    = summary(GAM_double_vil_plac_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_plac_mod)$coef[1], 
         sigma = summary(GAM_double_vil_plac_mod)$coef[2], 
         nu    = summary(GAM_double_vil_plac_mod)$df[2])


########################################
## All samples (village strat, placebo); triple

c( summary(GAM_triple_vil_plac_mod)$coef[1],
   summary(GAM_triple_vil_plac_mod)$coef[1] - 1.96*summary(GAM_triple_vil_plac_mod)$coef[2], 
   summary(GAM_triple_vil_plac_mod)$coef[1] + 1.96*summary(GAM_triple_vil_plac_mod)$coef[2] )


c( summary(GAM_triple_vil_plac_mod)$coef[1],
   summary(GAM_triple_vil_plac_mod)$coef[1] - 1.96*summary(GAM_triple_vil_plac_mod)$coef[2], 
   summary(GAM_triple_vil_plac_mod)$coef[1] + 1.96*summary(GAM_triple_vil_plac_mod)$coef[2] )/summary(GAM_double_vil_plac_mod)$coef[1]^2

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_plac_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_plac_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_plac_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_plac_mod)$df[2])




########################################
## Baseline samples; double                                   

c( summary(GAM_double_base_mod)$coef[1],
   summary(GAM_double_base_mod)$coef[1] - 1.96*summary(GAM_double_base_mod)$coef[2], 
   summary(GAM_double_base_mod)$coef[1] + 1.96*summary(GAM_double_base_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_base_mod)$coef[1], 
         sigma = summary(GAM_double_base_mod)$coef[2], 
         nu    = summary(GAM_double_base_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_base_mod)$coef[1], 
         sigma = summary(GAM_double_base_mod)$coef[2], 
         nu    = summary(GAM_double_base_mod)$df[2])


########################################
## Baseline samples; triple

c( summary(GAM_triple_base_mod)$coef[1],
   summary(GAM_triple_base_mod)$coef[1] - 1.96*summary(GAM_triple_base_mod)$coef[2], 
   summary(GAM_triple_base_mod)$coef[1] + 1.96*summary(GAM_triple_base_mod)$coef[2] )

c( summary(GAM_triple_base_mod)$coef[1],
   summary(GAM_triple_base_mod)$coef[1] - 1.96*summary(GAM_triple_base_mod)$coef[2], 
   summary(GAM_triple_base_mod)$coef[1] + 1.96*summary(GAM_triple_base_mod)$coef[2] )/summary(GAM_double_base_mod)$coef[1]^2

2*pst( q     = 0, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])

2*pst( q     = summary(GAM_double_base_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_base_mod)$coef[1], 
       sigma = summary(GAM_triple_base_mod)$coef[2], 
       nu    = summary(GAM_triple_base_mod)$df[2])


########################################
## Baseline samples (village strat); double                                   

c( summary(GAM_double_vil_base_mod)$coef[1],
   summary(GAM_double_vil_base_mod)$coef[1] - 1.96*summary(GAM_double_vil_base_mod)$coef[2], 
   summary(GAM_double_vil_base_mod)$coef[1] + 1.96*summary(GAM_double_vil_base_mod)$coef[2] )

2*pst(0, mu    = summary(GAM_double_vil_base_mod)$coef[1], 
         sigma = summary(GAM_double_vil_base_mod)$coef[2], 
         nu    = summary(GAM_double_vil_base_mod)$df[2])

2*pst(1, mu    = summary(GAM_double_vil_base_mod)$coef[1], 
         sigma = summary(GAM_double_vil_base_mod)$coef[2], 
         nu    = summary(GAM_double_vil_base_mod)$df[2])


########################################
## Baseline samples (village strat); triple

c( summary(GAM_triple_vil_base_mod)$coef[1],
   summary(GAM_triple_vil_base_mod)$coef[1] - 1.96*summary(GAM_triple_vil_base_mod)$coef[2], 
   summary(GAM_triple_vil_base_mod)$coef[1] + 1.96*summary(GAM_triple_vil_base_mod)$coef[2] )

c( summary(GAM_triple_vil_base_mod)$coef[1],
   summary(GAM_triple_vil_base_mod)$coef[1] - 1.96*summary(GAM_triple_vil_base_mod)$coef[2], 
   summary(GAM_triple_vil_base_mod)$coef[1] + 1.96*summary(GAM_triple_vil_base_mod)$coef[2] )/summary(GAM_double_vil_base_mod)$coef[1]^2

2*pst( q     = 0, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])

2*pst( q     = 1, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])

2*pst( q     = summary(GAM_double_vil_base_mod)$coef[1]^2, 
       mu    = summary(GAM_triple_vil_base_mod)$coef[1], 
       sigma = summary(GAM_triple_vil_base_mod)$coef[2], 
       nu    = summary(GAM_triple_vil_base_mod)$df[2])





