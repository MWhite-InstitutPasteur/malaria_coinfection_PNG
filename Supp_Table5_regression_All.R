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

###################################
## Read in data

ALB_par <- read.csv("Albinama_data_PAR.csv", sep=";")

ALB_epi <- read.csv("Albimana_data_EPI.csv")


###################################
## Process data frame

N_sam <- nrow(ALB_par)

ALB_data <- matrix(NA, nrow=N_sam, ncol=27)
colnames(ALB_data) <- c("studyid", "obs_id", "study_day", 
                        "age", "sex", "village", "bednet", "treat",
                        "qmal_pos",
                        "Pf_pos", "Pv_pos", "Po_pos", "Pm_pos",
                        "Pfg_pos", "Pvg_pos", "Pog_pos", "Pmg_pos",
                        "N_P_pos", "N_Pg_pos",
                        "sing_inf", "doub_inf", "trip_inf", "quad_inf",
                        "sing_inf_gam", "doub_inf_gam", "trip_inf_gam", "quad_inf_gam" )

ALB_data <- as.data.frame(ALB_data)

ALB_data[,1:8] <- ALB_epi[,c(1,2,3,67,61,69,64,65)]

for(i in 1:N_sam)
{
	index_i <- which( ALB_par[,2] == ALB_data[i,2] )

	ALB_data[i,9:17] <- ALB_par[index_i,c(9, 11, 13, 23, 19, 17, 15, 25, 21)]
}

drop <- c()

for(i in 1:N_sam)
{
	if( is.na(ALB_data[i,9]) == FALSE )
	{
		if( ALB_data[i,9] == 0 )
		{
			ALB_data[i,10:17] <- 0
		}
	}else{
		drop <- c(drop, i)
	}
}


ALB_data <- ALB_data[-drop,]


N_sam <- nrow(ALB_data)


for(i in 1:N_sam)
{
	na_index <- which(is.na(ALB_data[i,10:17])==TRUE)
	
	if( length(na_index) > 0 )
	{
		ALB_data[i,9+na_index] <- 0
	}
}

ALB_data[,20:27] <- 0

for(i in 1:N_sam)
{
	ALB_data[i,18] <- sum( ALB_data[i,10:13] )
	ALB_data[i,19] <- sum( ALB_data[i,14:17] )

	if( ALB_data[i,18] > 0 ){ ALB_data[i,20] <- 1 }
	if( ALB_data[i,18] > 1 ){ ALB_data[i,21] <- 1 }
	if( ALB_data[i,18] > 2 ){ ALB_data[i,22] <- 1 }
	if( ALB_data[i,18] > 3 ){ ALB_data[i,23] <- 1 }

	if( ALB_data[i,19] > 0 ){ ALB_data[i,24] <- 1 }
	if( ALB_data[i,19] > 1 ){ ALB_data[i,25] <- 1 }
	if( ALB_data[i,19] > 2 ){ ALB_data[i,26] <- 1 }
	if( ALB_data[i,19] > 3 ){ ALB_data[i,27] <- 1 }
}


ALB_data$village[which(ALB_data$village == "Balif ")] <- "Balif"


################################################
## Set the reference to be a female in Albinama
## with no bed net who did not receive PQ

ALB_data$village <- relevel( ALB_data$village, ref="Albinama")

ALB_data$sex <- relevel( ALB_data$sex, ref="female")

ALB_data$bednet <- relevel( ALB_data$bednet, ref="no")

ALB_data$treat <- relevel( ALB_data$treat, ref="PLACEBO")

ALB_data$age <- relevel( ALB_data$age, ref=10)


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
## P. falciaprum

Pf_glm <- glm( Pf_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pf_glm)

cbind(
round( cbind(  exp( summary(Pf_glm)$coef[,1] ),
               exp( summary(Pf_glm)$coef[,1] - 1.96*summary(Pf_glm)$coef[,2] ),
               exp( summary(Pf_glm)$coef[,1] + 1.96*summary(Pf_glm)$coef[,2] ) ), 3),

round( summary(Pf_glm)$coef[,4], 6 )
)


#########################################
## P. vivax

Pv_glm <- glm( Pv_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pv_glm)

cbind(
round( cbind(  exp( summary(Pv_glm)$coef[,1] ),
               exp( summary(Pv_glm)$coef[,1] - 1.96*summary(Pv_glm)$coef[,2] ),
               exp( summary(Pv_glm)$coef[,1] + 1.96*summary(Pv_glm)$coef[,2] ) ), 3),

round( summary(Pv_glm)$coef[,4], 6 )
)


#########################################
## P. malariae

Pm_glm <- glm( Pm_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pm_glm)

cbind(
round( cbind(  exp( summary(Pm_glm)$coef[,1] ),
               exp( summary(Pm_glm)$coef[,1] - 1.96*summary(Pm_glm)$coef[,2] ),
               exp( summary(Pm_glm)$coef[,1] + 1.96*summary(Pm_glm)$coef[,2] ) ), 3),

round( summary(Pm_glm)$coef[,4], 6 )
)


#########################################
## P. ovale

Po_glm <- glm( Po_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Po_glm)

cbind(
round( cbind(  exp( summary(Po_glm)$coef[,1] ),
               exp( summary(Po_glm)$coef[,1] - 1.96*summary(Po_glm)$coef[,2] ),
               exp( summary(Po_glm)$coef[,1] + 1.96*summary(Po_glm)$coef[,2] ) ), 3),

round( summary(Po_glm)$coef[,4], 6 )
)


#########################################
## P. falciaprum gametocytes

Pfg_glm <- glm( Pfg_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pfg_glm)

cbind(
round( cbind(  exp( summary(Pfg_glm)$coef[,1] ),
               exp( summary(Pfg_glm)$coef[,1] - 1.96*summary(Pfg_glm)$coef[,2] ),
               exp( summary(Pfg_glm)$coef[,1] + 1.96*summary(Pfg_glm)$coef[,2] ) ), 3),

round( summary(Pfg_glm)$coef[,4], 6 )
)


#########################################
## P. vivax gametocytes

Pvg_glm <- glm( Pvg_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pvg_glm)

cbind(
round( cbind(  exp( summary(Pvg_glm)$coef[,1] ),
               exp( summary(Pvg_glm)$coef[,1] - 1.96*summary(Pvg_glm)$coef[,2] ),
               exp( summary(Pvg_glm)$coef[,1] + 1.96*summary(Pvg_glm)$coef[,2] ) ), 3),

round( summary(Pvg_glm)$coef[,4], 6 )
)


#########################################
## P. malariae gametocytes

Pmg_glm <- glm( Pmg_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pmg_glm)

cbind(
round( cbind(  exp( summary(Pmg_glm)$coef[,1] ),
               exp( summary(Pmg_glm)$coef[,1] - 1.96*summary(Pmg_glm)$coef[,2] ),
               exp( summary(Pmg_glm)$coef[,1] + 1.96*summary(Pmg_glm)$coef[,2] ) ), 3),

round( summary(Pmg_glm)$coef[,4], 6 )
)


#########################################
## P. ovale gametocytes

Pog_glm <- glm( Pog_pos ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(Pog_glm)

cbind(
round( cbind(  exp( summary(Pog_glm)$coef[,1] ),
               exp( summary(Pog_glm)$coef[,1] - 1.96*summary(Pog_glm)$coef[,2] ),
               exp( summary(Pog_glm)$coef[,1] + 1.96*summary(Pog_glm)$coef[,2] ) ), 3),

round( summary(Pog_glm)$coef[,4], 6 )
)




#########################################
## Single infection (asexuals)

sing_inf_glm <- glm( sing_inf ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(sing_inf_glm)

cbind(
round( cbind(  exp( summary(sing_inf_glm)$coef[,1] ),
               exp( summary(sing_inf_glm)$coef[,1] - 1.96*summary(sing_inf_glm)$coef[,2] ),
               exp( summary(sing_inf_glm)$coef[,1] + 1.96*summary(sing_inf_glm)$coef[,2] ) ), 3),

round( summary(sing_inf_glm)$coef[,4], 6 )
)


#########################################
## Double infection (asexuals)

doub_inf_glm <- glm( doub_inf ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(doub_inf_glm)

cbind(
round( cbind(  exp( summary(doub_inf_glm)$coef[,1] ),
               exp( summary(doub_inf_glm)$coef[,1] - 1.96*summary(doub_inf_glm)$coef[,2] ),
               exp( summary(doub_inf_glm)$coef[,1] + 1.96*summary(doub_inf_glm)$coef[,2] ) ), 3),

round( summary(doub_inf_glm)$coef[,4], 6 )
)


#########################################
## Triple infection (asexuals)

trip_inf_glm <- glm( trip_inf ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(trip_inf_glm)

cbind(
round( cbind(  exp( summary(trip_inf_glm)$coef[,1] ),
               exp( summary(trip_inf_glm)$coef[,1] - 1.96*summary(trip_inf_glm)$coef[,2] ),
               exp( summary(trip_inf_glm)$coef[,1] + 1.96*summary(trip_inf_glm)$coef[,2] ) ), 3),

round( summary(trip_inf_glm)$coef[,4], 6 )
)


#########################################
## Quadruple infection (asexuals)

quad_inf_glm <- glm( quad_inf ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(quad_inf_glm)



cbind(
round( cbind(  exp( summary(quad_inf_glm)$coef[,1] ),
               exp( summary(quad_inf_glm)$coef[,1] - 1.96*summary(quad_inf_glm)$coef[,2] ),
               exp( summary(quad_inf_glm)$coef[,1] + 1.96*summary(quad_inf_glm)$coef[,2] ) ), 3),

round( summary(quad_inf_glm)$coef[,4], 6 )
)



#########################################
## Single infection (gametocytes)

sing_inf_gam_glm <- glm( sing_inf_gam ~ age + sex + village + bednet + treat, 
                    data = ALB_data, family = "binomial")

summary(sing_inf_gam_glm)

cbind(
round( cbind(  exp( summary(sing_inf_gam_glm)$coef[,1] ),
               exp( summary(sing_inf_gam_glm)$coef[,1] - 1.96*summary(sing_inf_gam_glm)$coef[,2] ),
               exp( summary(sing_inf_gam_glm)$coef[,1] + 1.96*summary(sing_inf_gam_glm)$coef[,2] ) ), 3),

round( summary(sing_inf_gam_glm)$coef[,4], 6 )
)


#########################################
## Double infection (gametocytes

doub_inf_gam_glm <- glm( doub_inf_gam ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(doub_inf_gam_glm)

cbind(
round( cbind(  exp( summary(doub_inf_gam_glm)$coef[,1] ),
               exp( summary(doub_inf_gam_glm)$coef[,1] - 1.96*summary(doub_inf_gam_glm)$coef[,2] ),
               exp( summary(doub_inf_gam_glm)$coef[,1] + 1.96*summary(doub_inf_gam_glm)$coef[,2] ) ), 3),

round( summary(doub_inf_gam_glm)$coef[,4], 6 )
)


#########################################
## Triple infection (gametocytes

trip_inf_gam_glm <- glm( trip_inf_gam ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(trip_inf_gam_glm)

cbind(
round( cbind(  exp( summary(trip_inf_gam_glm)$coef[,1] ),
               exp( summary(trip_inf_gam_glm)$coef[,1] - 1.96*summary(trip_inf_gam_glm)$coef[,2] ),
               exp( summary(trip_inf_gam_glm)$coef[,1] + 1.96*summary(trip_inf_gam_glm)$coef[,2] ) ), 3),

round( summary(trip_inf_gam_glm)$coef[,4], 6 )
)


#########################################
## Quadruple infection (gametocytes

quad_inf_gam_glm <- glm( quad_inf_gam ~ age + sex + village + bednet + treat, 
                     data = ALB_data, family = "binomial")

summary(quad_inf_gam_glm)

cbind(
round( cbind(  exp( summary(quad_inf_gam_glm)$coef[,1] ),
               exp( summary(quad_inf_gam_glm)$coef[,1] - 1.96*summary(quad_inf_gam_glm)$coef[,2] ),
               exp( summary(quad_inf_gam_glm)$coef[,1] + 1.96*summary(quad_inf_gam_glm)$coef[,2] ) ), 3),

round( summary(quad_inf_gam_glm)$coef[,4], 6 )
)



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
##  Epi overview                      ##
##                                    ##
########################################

studyid_part <- unique(ALB_data$studyid)
N_part <- length(studyid_part)

age_part <- rep(NA, N_part)
bednet_part <- rep(NA, N_part)
treat_part <- rep(NA, N_part)
sex_part <- rep(NA, N_part)


for(i in 1:N_part)
{
	age_part[i] <- ALB_data$age[which(ALB_data$studyid == studyid_part[i])][1]	

	bednet_part[i] <- as.character(ALB_data$bednet[which(ALB_data$studyid == studyid_part[i])][1])	

	treat_part[i] <- as.character(ALB_data$treat[which(ALB_data$studyid == studyid_part[i])][1])	

	sex_part[i] <- as.character(ALB_data$sex[which(ALB_data$studyid == studyid_part[i])][1])	
}

quantile( age_part, prob=c(0.5, 0.025, 0.975) )

table( bednet_part )

table( treat_part )

table( sex_part )


Pf_part <- rep(NA, N_part)
Pv_part <- rep(NA, N_part)
Pm_part <- rep(NA, N_part)
Po_part <- rep(NA, N_part)


Pfg_part <- rep(NA, N_part)
Pvg_part <- rep(NA, N_part)
Pmg_part <- rep(NA, N_part)
Pog_part <- rep(NA, N_part)


for(i in 1:N_part)
{
	if( sum(ALB_data$Pf_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pf_part[i] <- 1	
	}else{
		Pf_part[i] <- 0	
	}

	if( sum(ALB_data$Pv_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pv_part[i] <- 1	
	}else{
		Pv_part[i] <- 0	
	}

	if( sum(ALB_data$Pm_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pm_part[i] <- 1	
	}else{
		Pm_part[i] <- 0	
	}

	if( sum(ALB_data$Po_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Po_part[i] <- 1	
	}else{
		Po_part[i] <- 0	
	}

	if( sum(ALB_data$Pfg_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pfg_part[i] <- 1	
	}else{
		Pfg_part[i] <- 0	
	}

	if( sum(ALB_data$Pvg_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pvg_part[i] <- 1	
	}else{
		Pvg_part[i] <- 0	
	}

	if( sum(ALB_data$Pmg_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pmg_part[i] <- 1	
	}else{
		Pmg_part[i] <- 0	
	}

	if( sum(ALB_data$Pog_pos[which(ALB_data$studyid == studyid_part[i])]) > 0 )
	{
		Pog_part[i] <- 1	
	}else{
		Pog_part[i] <- 0	
	}
}


table( Pf_part )
table( Pv_part )
table( Pm_part )
table( Po_part )



table( Pfg_part )
table( Pvg_part )
table( Pmg_part )
table( Pog_part )

                                    



########################################
## Asexual infection table

Tab2 <- matrix(NA, ncol=16, nrow=9)
rownames(Tab2) <- c("age", "gender", "bednet", "treatment",
                    "Albinama", "Amahup", "Balanga", "Balif", "Bolumita")
colnames(Tab2) <- c("med_Pf", "low_Pf", "high_Pf", "P_Pf",
                    "med_Pv", "low_Pv", "high_Pv", "P_Pv",
                    "med_Pm", "low_Pm", "high_Pm", "P_Pm",
                    "med_Po", "low_Po", "high_Po", "P_Po")

index <- c(2, 3, 8, 9, 1, 4, 5, 6, 7)



Tab2[,1:4] <- c( exp( summary(Pf_glm)$coef[index,1] ),
	           exp( summary(Pf_glm)$coef[index,1] - 1.96*summary(Pf_glm)$coef[index,2] ),
	           exp( summary(Pf_glm)$coef[index,1] + 1.96*summary(Pf_glm)$coef[index,2] ),
	           summary(Pf_glm)$coef[index,4] )


Tab2[,5:8] <- c( exp( summary(Pv_glm)$coef[index,1] ),
	           exp( summary(Pv_glm)$coef[index,1] - 1.96*summary(Pv_glm)$coef[index,2] ),
	           exp( summary(Pv_glm)$coef[index,1] + 1.96*summary(Pv_glm)$coef[index,2] ),
	           summary(Pv_glm)$coef[index,4] )


Tab2[,9:12] <- c( exp( summary(Pm_glm)$coef[index,1] ),
	            exp( summary(Pm_glm)$coef[index,1] - 1.96*summary(Pm_glm)$coef[index,2] ),
	            exp( summary(Pm_glm)$coef[index,1] + 1.96*summary(Pm_glm)$coef[index,2] ),
	            summary(Pm_glm)$coef[index,4] )


Tab2[,13:16] <- c( exp( summary(Po_glm)$coef[index,1] ),
	             exp( summary(Po_glm)$coef[index,1] - 1.96*summary(Po_glm)$coef[index,2] ),
	             exp( summary(Po_glm)$coef[index,1] + 1.96*summary(Po_glm)$coef[index,2] ),
	             summary(Po_glm)$coef[index,4] )


Tab2[,c(1:3,5:7,9:11,13:15)] <- round( Tab2[,c(1:3,5:7,9:11,13:15)], 3)


########################################
## Gametocyte infection table

Tab3 <- matrix(NA, ncol=16, nrow=9)
rownames(Tab3) <- c("age", "gender", "bednet", "treatment",
                    "Albinama", "Amahup", "Balanga", "Balif", "Bolumita")
colnames(Tab3) <- c("med_Pfg", "low_Pfg", "high_Pfg", "P_Pfg",
                    "med_Pvg", "low_Pvg", "high_Pvg", "P_Pvg",
                    "med_Pmg", "low_Pmg", "high_Pmg", "P_Pmg",
                    "med_Pog", "low_Pog", "high_Pog", "P_Pog")

index <- c(2, 3, 8, 9, 1, 4, 5, 6, 7)



Tab3[,1:4] <- c( exp( summary(Pfg_glm)$coef[index,1] ),
	           exp( summary(Pfg_glm)$coef[index,1] - 1.96*summary(Pfg_glm)$coef[index,2] ),
	           exp( summary(Pfg_glm)$coef[index,1] + 1.96*summary(Pfg_glm)$coef[index,2] ),
	           summary(Pfg_glm)$coef[index,4] )


Tab3[,5:8] <- c( exp( summary(Pvg_glm)$coef[index,1] ),
	           exp( summary(Pvg_glm)$coef[index,1] - 1.96*summary(Pvg_glm)$coef[index,2] ),
	           exp( summary(Pvg_glm)$coef[index,1] + 1.96*summary(Pvg_glm)$coef[index,2] ),
	           summary(Pvg_glm)$coef[index,4] )


Tab3[,9:12] <- c( exp( summary(Pmg_glm)$coef[index,1] ),
	            exp( summary(Pmg_glm)$coef[index,1] - 1.96*summary(Pmg_glm)$coef[index,2] ),
	            exp( summary(Pmg_glm)$coef[index,1] + 1.96*summary(Pmg_glm)$coef[index,2] ),
	            summary(Pmg_glm)$coef[index,4] )


Tab3[,13:16] <- c( exp( summary(Pog_glm)$coef[index,1] ),
	             exp( summary(Pog_glm)$coef[index,1] - 1.96*summary(Pog_glm)$coef[index,2] ),
	             exp( summary(Pog_glm)$coef[index,1] + 1.96*summary(Pog_glm)$coef[index,2] ),
	             summary(Pog_glm)$coef[index,4] )


Tab3[,c(1:3,5:7,9:11,13:15)] <- round( Tab3[,c(1:3,5:7,9:11,13:15)], 3)





########################################
## Asexual coinfection table

Tab4 <- matrix(NA, ncol=16, nrow=9)
rownames(Tab4) <- c("age", "gender", "bednet", "treatment",
                    "Albinama", "Amahup", "Balanga", "Balif", "Bolumita")
colnames(Tab4) <- c("med_sing", "low_sing", "high_sing", "P_sing",
                    "med_doub", "low_doub", "high_doub", "P_doub",
                    "med_trip", "low_trip", "high_trip", "P_trip",
                    "med_quad", "low_quad", "high_quad", "P_quad")

index <- c(2, 3, 8, 9, 1, 4, 5, 6, 7)
	


Tab4[,1:4] <- c( exp( summary(sing_inf_glm)$coef[index,1] ),
	           exp( summary(sing_inf_glm)$coef[index,1] - 1.96*summary(sing_inf_glm)$coef[index,2] ),
	           exp( summary(sing_inf_glm)$coef[index,1] + 1.96*summary(sing_inf_glm)$coef[index,2] ),
	           summary(sing_inf_glm)$coef[index,4] )


Tab4[,5:8] <- c( exp( summary(doub_inf_glm)$coef[index,1] ),
	           exp( summary(doub_inf_glm)$coef[index,1] - 1.96*summary(doub_inf_glm)$coef[index,2] ),
	           exp( summary(doub_inf_glm)$coef[index,1] + 1.96*summary(doub_inf_glm)$coef[index,2] ),
	           summary(doub_inf_glm)$coef[index,4] )


Tab4[,9:12] <- c( exp( summary(trip_inf_glm)$coef[index,1] ),
	            exp( summary(trip_inf_glm)$coef[index,1] - 1.96*summary(trip_inf_glm)$coef[index,2] ),
	            exp( summary(trip_inf_glm)$coef[index,1] + 1.96*summary(trip_inf_glm)$coef[index,2] ),
	            summary(trip_inf_glm)$coef[index,4] )


Tab4[,13:16] <- c( exp( summary(quad_inf_glm)$coef[index,1] ),
	             exp( summary(quad_inf_glm)$coef[index,1] - 1.96*summary(quad_inf_glm)$coef[index,2] ),
	             exp( summary(quad_inf_glm)$coef[index,1] + 1.96*summary(quad_inf_glm)$coef[index,2] ),
	             summary(quad_inf_glm)$coef[index,4] )


Tab4[,c(1:3,5:7,9:11,13:15)] <- round( Tab4[,c(1:3,5:7,9:11,13:15)], 3)


########################################
## Gametocyte coinfection table

Tab5 <- matrix(NA, ncol=16, nrow=9)
rownames(Tab5) <- c("age", "gender", "bednet", "treatment",
                        "Albinama", "Amahup", "Balanga", "Balif", "Bolumita")
colnames(Tab5) <- c("med_sing", "low_sing", "high_sing", "P_sing",
                    "med_doub", "low_doub", "high_doub", "P_doub",
                    "med_trip", "low_trip", "high_trip", "P_trip",
                    "med_quad", "low_quad", "high_quad", "P_quad")


index <- c(2, 3, 8, 9, 1, 4, 5, 6, 7)


Tab5[,1:4] <- c( exp( summary(sing_inf_gam_glm)$coef[index,1] ),
	           exp( summary(sing_inf_gam_glm)$coef[index,1] - 1.96*summary(sing_inf_gam_glm)$coef[index,2] ),
	           exp( summary(sing_inf_gam_glm)$coef[index,1] + 1.96*summary(sing_inf_gam_glm)$coef[index,2] ),
	           summary(sing_inf_gam_glm)$coef[index,4] )


Tab5[,5:8] <- c( exp( summary(doub_inf_gam_glm)$coef[index,1] ),
	           exp( summary(doub_inf_gam_glm)$coef[index,1] - 1.96*summary(doub_inf_gam_glm)$coef[index,2] ),
	           exp( summary(doub_inf_gam_glm)$coef[index,1] + 1.96*summary(doub_inf_gam_glm)$coef[index,2] ),
	           summary(doub_inf_gam_glm)$coef[index,4] )


Tab5[,9:12] <- c( exp( summary(trip_inf_gam_glm)$coef[index,1] ),
	            exp( summary(trip_inf_gam_glm)$coef[index,1] - 1.96*summary(trip_inf_gam_glm)$coef[index,2] ),
	            exp( summary(trip_inf_gam_glm)$coef[index,1] + 1.96*summary(trip_inf_gam_glm)$coef[index,2] ),
	            summary(trip_inf_gam_glm)$coef[index,4] )


Tab5[,13:16] <- c( exp( summary(quad_inf_gam_glm)$coef[index,1] ),
	             exp( summary(quad_inf_gam_glm)$coef[index,1] - 1.96*summary(quad_inf_gam_glm)$coef[index,2] ),
	             exp( summary(quad_inf_gam_glm)$coef[index,1] + 1.96*summary(quad_inf_gam_glm)$coef[index,2] ),
	             summary(quad_inf_gam_glm)$coef[index,4] )


Tab5[,c(1:3,5:7,9:11,13:15)] <- round( Tab5[,c(1:3,5:7,9:11,13:15)], 3)





