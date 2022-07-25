library(VennDiagram)
library(binom)
library(fields)



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
	ALB_par_village[i] <- as.character(ALB_epi$villagegroup_corrected[which(ALB_epi$studyid == ALB_par$studyid[i])][1])
}

ALB_par$village <- ALB_par_village



ALB_par <- ALB_par[which(ALB_par$study_day<0),]

N_data <- nrow(ALB_par)


paras_mod <- function( x_ref, x_compare )
{
	x_ref <- x_ref[which(x_ref>0)]
	x_compare <- x_compare[which(x_compare>0)]


	paras_dens <- log(c( x_ref, x_compare ))

	paras_factor <- as.factor(c( rep(0, length(x_ref)), rep(1, length(x_compare)) ))
	paras_factor <- relevel(paras_factor, ref="0")

	paras_lm <- summary(lm( paras_dens ~ paras_factor ))		

	paras_effect <- exp( c( paras_lm$coef[2,1],
	                        paras_lm$coef[2,1] - 1.96*paras_lm$coef[2,2],
	                        paras_lm$coef[2,1] + 1.96*paras_lm$coef[2,2] ) )

	paras_pval <- paras_lm$coef[2,4]

	list( paras_effect, paras_pval )
}



COINF <- array(NA, dim=c(8,8,3) )
colnames(COINF) <- c("Pf", "Pv", "Pm", "Po", "Pfg", "Pvg", "Pmg", "Pog")
rownames(COINF) <- c("Pf", "Pv", "Pm", "Po", "Pfg", "Pvg", "Pmg", "Pog")

####################################
## P. falciparum asexuals

COINF[1,2,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[1,3,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[1,4,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[1,5,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[1,6,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[1,7,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[1,8,] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. vivax asexuals

COINF[2,1,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[2,3,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[2,4,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[2,5,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[2,6,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[2,7,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[2,8,] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. malariae asexuals

COINF[3,1,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[3,2,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[3,4,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[3,5,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[3,6,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[3,7,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[3,8,] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. ovale asexuals

COINF[4,1,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[4,2,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[4,3,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[4,5,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[4,6,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[4,7,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[4,8,] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. falciparum gametocytes

COINF[5,1,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[5,2,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[5,3,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[5,4,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[5,6,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[5,7,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[5,8,] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. vivax gametocytes

COINF[6,1,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[6,2,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[6,3,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[6,4,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[6,5,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[6,7,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1)] )[[1]]

COINF[6,8,] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. malariae gametocytes

COINF[7,1,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[7,2,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[7,3,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[7,4,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[7,5,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[7,6,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[7,8,] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1)] )[[1]]


####################################
## P. ovale gametocytes

COINF[8,1,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pf_pos==1)] )[[1]]

COINF[8,2,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pv_pos==1)] )[[1]]

COINF[8,3,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pm_pos==1)] )[[1]]

COINF[8,4,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$po_pos==1)] )[[1]]

COINF[8,5,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pfg_pos==1)] )[[1]]

COINF[8,6,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pvg_pos==1)] )[[1]]

COINF[8,7,] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pmg_pos==1)] )[[1]]




COINF_p <- matrix(NA, nrow=8, ncol=8)
colnames(COINF_p) <- c("Pf", "Pv", "Pm", "Po", "Pfg", "Pvg", "Pmg", "Pog")
rownames(COINF_p) <- c("Pf", "Pv", "Pm", "Po", "Pfg", "Pvg", "Pmg", "Pog")

####################################
## P. falciparum asexuals

COINF_p[1,2] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[1,3] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[1,4] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[1,5] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[1,6] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[1,7] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[1,8] <- paras_mod( x_ref     = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pf_copyno[which(ALB_par$pf_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. vivax asexuals

COINF_p[2,1] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[2,3] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[2,4] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[2,5] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[2,6] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[2,7] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[2,8] <- paras_mod( x_ref     = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pv_copyno[which(ALB_par$pv_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. malariae asexuals

COINF_p[3,1] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[3,2] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[3,4] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[3,5] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[3,6] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[3,7] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[3,8] <- paras_mod( x_ref     = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pm_copyno[which(ALB_par$pm_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. ovale asexuals

COINF_p[4,1] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[4,2] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[4,3] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[4,5] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[4,6] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[4,7] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[4,8] <- paras_mod( x_ref     = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$po_copyno[which(ALB_par$po_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. falciparum gametocytes

COINF_p[5,1] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[5,2] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[5,3] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[5,4] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[5,6] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[5,7] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[5,8] <- paras_mod( x_ref     = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. vivax gametocytes

COINF_p[6,1] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[6,2] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[6,3] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[6,4] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[6,5] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pfg_pos==0)],
                         x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[6,7] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1)] )[[2]]

COINF_p[6,8] <- paras_mod( x_ref     = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. malariae gametocytes

COINF_p[7,1] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[7,2] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[7,3] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[7,4] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[7,5] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[7,6] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[7,8] <- paras_mod( x_ref     = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==0)],
                          x_compare = ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1)] )[[2]]


####################################
## P. ovale gametocytes

COINF_p[8,1] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pf_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pf_pos==1)] )[[2]]

COINF_p[8,2] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pv_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pv_pos==1)] )[[2]]

COINF_p[8,3] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pm_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pm_pos==1)] )[[2]]

COINF_p[8,4] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$po_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$po_pos==1)] )[[2]]

COINF_p[8,5] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pfg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pfg_pos==1)] )[[2]]

COINF_p[8,6] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pvg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pvg_pos==1)] )[[2]]

COINF_p[8,7] <- paras_mod( x_ref     = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pmg_pos==0)],
                          x_compare = ALB_par$pog_copyno[which(ALB_par$pog_pos==1 & ALB_par$pmg_pos==1)] )[[2]]



###############################################
## Benjamini-Hochberg correction of P values

COINF_p_adjust <- matrix( p.adjust( COINF_p, method="BH" ), nrow=8, ncol=8)


COINF[,,1] <- t(COINF[,,1])
COINF[,,2] <- t(COINF[,,2])
COINF[,,3] <- t(COINF[,,3])

COINF_p <- t(COINF_p)
COINF_p_adjust <- t(COINF_p_adjust)



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
Pm_col <- "forestgreen"
Po_col <- "dodgerblue"

edge = 0.15

main.size = 1.1
axis.size = 1.1
lab.size  = 1.2

point.size = 0.5
line.size  = 1
rr.size1    = 0.005
rr.size2    = 0.0075
arrow.edge = 0.015


effect_bands <- c( 0.0, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 30 )
N_bands <- length(effect_bands)-1

##colour_bands <- rainbow( N_bands, start=4/6, end=2/6)

colour_bands <- c("navy", "mediumblue", "dodgerblue", "lightskyblue1",
                  "darkseagreen1", "greenyellow", "limegreen", "forestgreen" )

tiff( file="Figure4_coinfection_density_Baseline.tif", width=20, height=18, units="cm", res=500)



lay.mat <- rbind( c(1), 
                  c(2) )
layout(lay.mat, heights=c(10,1.75) )
layout.show(2)


par(mar=c(6.5,7,0.5,0.5))
par(mgp=c(0.75, 0.25, 0))


plot(x=100, y=100,
xlim=c(0,8), ylim=c(0,8),
xlab="", ylab="",
main="",
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n',
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:8)
{
	for(j in 1:8)
	{
		polygon( x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                     border=NA, col="grey" )

		if( is.na(COINF[i,j,1])==FALSE )
		{
			for( k in 1:N_bands )
			{

				if( (COINF[i,j,1] > effect_bands[k]) && (COINF[i,j,1] < effect_bands[k+1]) )
				{
					polygon( x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                        	         border=NA, col=colour_bands[k] )
				}
			}
		}

		if( is.na(COINF_p[i,j])==FALSE )
		{
			if( COINF_p[i,j] < 0.05 )
			{
				polygon( x=c(i-1,i,i,i-1), y=c(j-1,j-1,j-1+edge,j-1+edge),
            	               border=NA, col="orange" )

				polygon( x=c(i-1,i,i,i-1), y=c(j-edge,j-edge,j,j),
            	               border=NA, col="orange" )

				polygon( x=c(i-1,i-1+edge,i-1+edge,i-1), y=c(j-1,j-1,j,j),
            	               border=NA, col="orange" )

				polygon( x=c(i,i-edge,i-edge,i), y=c(j-1,j-1,j,j),
            	               border=NA, col="orange" )
			}
		}

		
		if( is.na(COINF_p_adjust[i,j])==FALSE )
		{
			if( COINF_p_adjust[i,j] < 0.05 )
			{
				polygon( x=c(i-1,i,i,i-1), y=c(j-1,j-1,j-1+edge,j-1+edge),
            	               border=NA, col="red" )

				polygon( x=c(i-1,i,i,i-1), y=c(j-edge,j-edge,j,j),
            	               border=NA, col="red" )

				polygon( x=c(i-1,i-1+edge,i-1+edge,i-1), y=c(j-1,j-1,j,j),
            	               border=NA, col="red" )

				polygon( x=c(i,i-edge,i-edge,i), y=c(j-1,j-1,j,j),
            	               border=NA, col="red" )
			}
		}

	}
}


axis(1, at=seq(from=0.5, by=1, length=8), 
        label=c("Pf asexual \n co-infection", "Pv asexual \n co-infection", "Pm asexual \n co-infection", "Po asexual \n co-infection",
                "Pf gametocyte \n co-infection", "Pv gametocyte \n co-infection", "Pm gametocyte \n co-infection", "Po gametocyte \n co-infection"), 
     las=2, tick=FALSE, cex.axis=axis.size)

axis(2, at=seq(from=0.5, by=1, length=8), 
        label=c("Pf asexual \n density", "Pv asexual \n density", "Pm asexual \n density", "Po asexual \n density",
                "Pf gametocyte \n density", "Pv gametocyte \n density", "Pm gametocyte \n density", "Po gametocyte \n density"),
        las=2, tick=FALSE, cex.axis=axis.size)



#############
## LEGEND


par(mar=c(3,3,1,1))
par(mgp=c(1.9, 0.75, 0))

plot(x=100, y=100,
xlim=c(0,100), ylim=c(0,1),
xlab="fold change in parasite density", ylab="",
main="",
cex.lab=1.5,
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n')

for(k in 1:N_bands)
{
	polygon(y=c(0,1,1,0), x=c(k-1,k-1,k,k)*100/N_bands,
              border=NA, col=colour_bands[k])	
}


axis(1, at=seq(from=0, to=100, length=N_bands+1), 
     label=effect_bands, 
     cex.axis=1.25)


dev.off()



Pfpos_Pmpos <- ALB_par$pf_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 ) ]

Pfpos_Pmneg <- ALB_par$pf_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pm_pos==0 ) ]


exp(mean(log(Pfpos_Pmpos)))

exp(mean(log(Pfpos_Pmneg)))




Pmpos_Pfpos <- ALB_par$pm_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pm_pos==1 ) ]

Pmpos_Pfneg <- ALB_par$pm_copyno[ which(ALB_par$pf_pos==0 & ALB_par$pm_pos==1 ) ]


exp(mean(log(Pmpos_Pfpos)))

exp(mean(log(Pmpos_Pfneg)))




Pvpos_Pvgpos <- ALB_par$pv_copyno[ which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==1 ) ]

Pvpos_Pvgneg <- ALB_par$pv_copyno[ which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==0 ) ]


exp(mean(log(Pvpos_Pvgpos)))

exp(mean(log(Pvpos_Pvgneg)))




Pvgpos_Pvpos <- ALB_par$pvg_copyno[ which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==1 ) ]

Pvgpos_Pvneg <- ALB_par$pvg_copyno[ which(ALB_par$pvg_pos==1 & ALB_par$pv_pos==0 ) ]


exp(mean(log(Pvgpos_Pvpos)))

exp(mean(log(Pvgpos_Pvneg)))








Pfpos_single <- ALB_par$pf_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pfpos_single, prob=c(0.5, 0.025, 0.975))

length(Pfpos_single)


Pfpos_co <- ALB_par$pf_copyno[ which(ALB_par$pf_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pfpos_co, prob=c(0.5, 0.025, 0.975))

length(Pfpos_co)






Pvpos_single <- ALB_par$pv_copyno[ which(ALB_par$pv_pos==1 & ALB_par$pf_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pvpos_single, prob=c(0.5, 0.025, 0.975))

length(Pvpos_single)


Pvpos_co <- ALB_par$pv_copyno[ which(ALB_par$pv_pos==1 & 
                                      (ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pvpos_co, prob=c(0.5, 0.025, 0.975))

length(Pvpos_co)




Pmpos_single <- ALB_par$pm_copyno[ which(ALB_par$pm_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pf_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pmpos_single, prob=c(0.5, 0.025, 0.975))

length(Pmpos_single)


Pmpos_co <- ALB_par$pm_copyno[ which(ALB_par$pm_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pf_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pmpos_co, prob=c(0.5, 0.025, 0.975))

length(Pmpos_co)




Popos_single <- ALB_par$po_copyno[ which(ALB_par$po_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$pf_pos==0 ) ]


quantile(Popos_single, prob=c(0.5, 0.025, 0.975))

length(Popos_single)


Popos_co <- ALB_par$po_copyno[ which(ALB_par$po_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$pf_pos==1) ) ]

quantile(Popos_co, prob=c(0.5, 0.025, 0.975))

length(Popos_co)






ALB_par$pog_pos[which( is.na(ALB_par$pog_pos) == TRUE )] <- 0
ALB_par$pmg_pos[which( is.na(ALB_par$pmg_pos) == TRUE )] <- 0


Pfgpos_single <- ALB_par$pfg_copyno[ which(ALB_par$pfg_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pfgpos_single, prob=c(0.5, 0.025, 0.975))

length(Pfgpos_single)


Pfgpos_co <- ALB_par$pfg_copyno[ which(ALB_par$pfg_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pfgpos_co, prob=c(0.5, 0.025, 0.975))

length(Pfgpos_co)






Pvgpos_single <- ALB_par$pvg_copyno[ which(ALB_par$pvg_pos==1 & ALB_par$pf_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pvgpos_single, prob=c(0.5, 0.025, 0.975))

length(Pvgpos_single)


Pvgpos_co <- ALB_par$pvg_copyno[ which(ALB_par$pvg_pos==1 & 
                                      (ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pvgpos_co, prob=c(0.5, 0.025, 0.975))

length(Pvgpos_co)




Pmgpos_single <- ALB_par$pmg_copyno[ which(ALB_par$pmg_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pf_pos==0 & ALB_par$po_pos==0 ) ]


quantile(Pmgpos_single, prob=c(0.5, 0.025, 0.975))

length(Pmgpos_single)


Pmgpos_co <- ALB_par$pmg_copyno[ which(ALB_par$pmg_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pf_pos==1 | ALB_par$po_pos==1) ) ]

quantile(Pmgpos_co, prob=c(0.5, 0.025, 0.975))

length(Pmgpos_co)




Pogpos_single <- ALB_par$pog_copyno[ which(ALB_par$pog_pos==1 & ALB_par$pv_pos==0 & 
                                         ALB_par$pm_pos==0 & ALB_par$pf_pos==0 ) ]


quantile(Pogpos_single, prob=c(0.5, 0.025, 0.975))

length(Pogpos_single)


Pogpos_co <- ALB_par$pog_copyno[ which(ALB_par$pog_pos==1 & 
                                      (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$pf_pos==1) ) ]

quantile(Pogpos_co, prob=c(0.5, 0.025, 0.975))

length(Pogpos_co)






length( which(ALB_par$pf_pos==1) )

length( which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==1) )

length( which(ALB_par$pf_pos==1 & ALB_par$pfg_pos==0) )


length(which( ALB_par$pf_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$po_pos==0 ))

length(which( ALB_par$pfg_pos==1 & ALB_par$pf_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$po_pos==0 ))


length(which( ALB_par$pf_pos==1 & 
       (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ))


length(which( ALB_par$pf_pos==1 & ALB_par$pfg_pos==1 &
       (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ))









length( which(ALB_par$pv_pos==1) )

length( which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==1) )

length( which(ALB_par$pv_pos==1 & ALB_par$pvg_pos==0) )


length(which( ALB_par$pv_pos==1 & ALB_par$pf_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$po_pos==0 ))

length(which( ALB_par$pvg_pos==1 & ALB_par$pv_pos==1 & ALB_par$pf_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$po_pos==0 ))


length(which( ALB_par$pv_pos==1 & 
       (ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ))


length(which( ALB_par$pv_pos==1 & ALB_par$pvg_pos==1 &
       (ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ))










length( which(ALB_par$pm_pos==1) )

length( which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==1) )

length( which(ALB_par$pm_pos==1 & ALB_par$pmg_pos==0) )


length(which( ALB_par$pm_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pf_pos==0 & ALB_par$po_pos==0 ))

length(which( ALB_par$pmg_pos==1 & ALB_par$pm_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pf_pos==0 & ALB_par$po_pos==0 ))


length(which( ALB_par$pm_pos==1 & 
       (ALB_par$pv_pos==1 | ALB_par$pf_pos==1 | ALB_par$po_pos==1) ))


length(which( ALB_par$pm_pos==1 & ALB_par$pmg_pos==1 &
       (ALB_par$pv_pos==1 | ALB_par$pf_pos==1 | ALB_par$po_pos==1) ))







length( which(ALB_par$po_pos==1) )

length( which(ALB_par$po_pos==1 & ALB_par$pog_pos==1) )

length( which(ALB_par$po_pos==1 & ALB_par$pog_pos==0) )


length(which( ALB_par$po_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$pf_pos==0 ))

length(which( ALB_par$pog_pos==1 & ALB_par$po_pos==1 & ALB_par$pv_pos==0 & 
              ALB_par$pm_pos==0 & ALB_par$pf_pos==0 ))


length(which( ALB_par$po_pos==1 & 
       (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$pf_pos==1) ))


length(which( ALB_par$po_pos==1 & ALB_par$pog_pos==1 &
       (ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$pf_pos==1) ))




