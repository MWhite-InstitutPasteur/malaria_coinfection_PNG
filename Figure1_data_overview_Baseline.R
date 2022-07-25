library(VennDiagram)
library(binom)

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


###################################
##                               ##
##  PART 1: Asexual prevalence   ##
##          	               ## 
###################################

Pf_prev_1 <- length( which(ALB_par$pf_pos==1) )/N_data
Pv_prev_1 <- length( which(ALB_par$pv_pos==1) )/N_data
Pm_prev_1 <- length( which(ALB_par$pm_pos==1) )/N_data
Po_prev_1 <- length( which(ALB_par$po_pos==1) )/N_data


Pf_prev_2 <- length( intersect( which(ALB_par$pf_pos==1), 
                                which(ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) )/N_data
Pv_prev_2 <- length( intersect( which(ALB_par$pv_pos==1), 
                                which(ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) )/N_data
Pm_prev_2 <- length( intersect( which(ALB_par$pm_pos==1), 
                                which(ALB_par$pf_pos==1 | ALB_par$pv_pos==1 | ALB_par$po_pos==1) ) )/N_data
Po_prev_2 <- length( intersect( which(ALB_par$po_pos==1), 
                                which(ALB_par$pf_pos==1 | ALB_par$pv_pos==1 | ALB_par$pm_pos==1) ) )/N_data

Pf_prev_3 <- length( intersect( which(ALB_par$pf_pos==1), 
                                which(rowSums(ALB_par[,c(13,19,23)])>=2) ) )/N_data
Pv_prev_3 <- length( intersect( which(ALB_par$pv_pos==1), 
                                which(rowSums(ALB_par[,c(11,19,23)])>=2) ) )/N_data
Pm_prev_3 <- length( intersect( which(ALB_par$pm_pos==1), 
                                which(rowSums(ALB_par[,c(11,13,23)])>=2) ) )/N_data
Po_prev_3 <- length( intersect( which(ALB_par$po_pos==1), 
                                which(rowSums(ALB_par[,c(11,13,19)])>=2) ) )/N_data

Pf_prev_4 <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )/N_data
Pv_prev_4 <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )/N_data
Pm_prev_4 <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )/N_data
Po_prev_4 <- length( which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) )/N_data



Pf_prev_1_CI <- as.numeric(binom.confint( Pf_prev_1*N_data, N_data, method="wilson")[1,5:6])
Pv_prev_1_CI <- as.numeric(binom.confint( Pv_prev_1*N_data, N_data, method="wilson")[1,5:6])
Pm_prev_1_CI <- as.numeric(binom.confint( Pm_prev_1*N_data, N_data, method="wilson")[1,5:6])
Po_prev_1_CI <- as.numeric(binom.confint( Po_prev_1*N_data, N_data, method="wilson")[1,5:6])


Pf_prev_2_CI <- as.numeric(binom.confint( Pf_prev_2*N_data, N_data, method="wilson")[1,5:6])
Pv_prev_2_CI <- as.numeric(binom.confint( Pv_prev_2*N_data, N_data, method="wilson")[1,5:6])
Pm_prev_2_CI <- as.numeric(binom.confint( Pm_prev_2*N_data, N_data, method="wilson")[1,5:6])
Po_prev_2_CI <- as.numeric(binom.confint( Po_prev_2*N_data, N_data, method="wilson")[1,5:6])


Pf_prev_3_CI <- as.numeric(binom.confint( Pf_prev_3*N_data, N_data, method="wilson")[1,5:6])
Pv_prev_3_CI <- as.numeric(binom.confint( Pv_prev_3*N_data, N_data, method="wilson")[1,5:6])
Pm_prev_3_CI <- as.numeric(binom.confint( Pm_prev_3*N_data, N_data, method="wilson")[1,5:6])
Po_prev_3_CI <- as.numeric(binom.confint( Po_prev_3*N_data, N_data, method="wilson")[1,5:6])


Pf_prev_4_CI <- as.numeric(binom.confint( Pf_prev_4*N_data, N_data, method="wilson")[1,5:6])
Pv_prev_4_CI <- as.numeric(binom.confint( Pv_prev_4*N_data, N_data, method="wilson")[1,5:6])
Pm_prev_4_CI <- as.numeric(binom.confint( Pm_prev_4*N_data, N_data, method="wilson")[1,5:6])
Po_prev_4_CI <- as.numeric(binom.confint( Po_prev_4*N_data, N_data, method="wilson")[1,5:6])


###################################
##                               ##
##  PART 2: Asexual density      ##
##          	               ## 
###################################

Pf_copy_1 <- ALB_par$pf_copyno[which(ALB_par$pf_pos==1)]
Pv_copy_1 <- ALB_par$pv_copyno[which(ALB_par$pv_pos==1)]
Pm_copy_1 <- ALB_par$pm_copyno[which(ALB_par$pm_pos==1)]
Po_copy_1 <- ALB_par$po_copyno[which(ALB_par$po_pos==1)]


Pf_copy_2 <- ALB_par$pf_copyno[ intersect( which(ALB_par$pf_pos==1), 
                                           which(ALB_par$pv_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]
Pv_copy_2 <- ALB_par$pv_copyno[ intersect( which(ALB_par$pv_pos==1), 
                                           which(ALB_par$pf_pos==1 | ALB_par$pm_pos==1 | ALB_par$po_pos==1) ) ]
Pm_copy_2 <- ALB_par$pm_copyno[ intersect( which(ALB_par$pm_pos==1), 
                                           which(ALB_par$pf_pos==1 | ALB_par$pv_pos==1 | ALB_par$po_pos==1) ) ]
Po_copy_2 <- ALB_par$po_copyno[ intersect( which(ALB_par$po_pos==1), 
                                           which(ALB_par$pf_pos==1 | ALB_par$pv_pos==1 | ALB_par$pm_pos==1) ) ]


Pf_copy_3 <- ALB_par$pf_copyno[ intersect( which(ALB_par$pf_pos==1), 
                                           which(rowSums(ALB_par[,c(13,19,23)])>=2) ) ]
Pv_copy_3 <- ALB_par$pv_copyno[ intersect( which(ALB_par$pv_pos==1), 
                                           which(rowSums(ALB_par[,c(11,19,23)])>=2) ) ]
Pm_copy_3 <- ALB_par$pm_copyno[ intersect( which(ALB_par$pm_pos==1), 
                                           which(rowSums(ALB_par[,c(11,13,23)])>=2) ) ]
Po_copy_3 <- ALB_par$po_copyno[ intersect( which(ALB_par$po_pos==1), 
                                           which(rowSums(ALB_par[,c(11,13,19)])>=2) ) ]

Pf_copy_4 <- ALB_par$pf_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) ]
Pv_copy_4 <- ALB_par$pv_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) ]
Pm_copy_4 <- ALB_par$pm_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) ]
Po_copy_4 <- ALB_par$po_copyno[ which(ALB_par$pf_pos==1 & ALB_par$pv_pos==1 & ALB_par$pm_pos==1 & ALB_par$po_pos==1) ]


boxplot_stats <- function( xx )
{
	qq <- quantile( xx, prob=c(0.025, 0.25, 0.5, 0.75, 0.975) )

	outliers <- xx[ which( xx < qq[1] | xx > qq[5]) ]

	list( qq, outliers )
}


Pf_copy_1_box <- boxplot_stats( Pf_copy_1 )
Pv_copy_1_box <- boxplot_stats( Pv_copy_1 )
Pm_copy_1_box <- boxplot_stats( Pm_copy_1 )
Po_copy_1_box <- boxplot_stats( Po_copy_1 )

Pf_copy_2_box <- boxplot_stats( Pf_copy_2 )
Pv_copy_2_box <- boxplot_stats( Pv_copy_2 )
Pm_copy_2_box <- boxplot_stats( Pm_copy_2 )
Po_copy_2_box <- boxplot_stats( Po_copy_2 )

Pf_copy_3_box <- boxplot_stats( Pf_copy_3 )
Pv_copy_3_box <- boxplot_stats( Pv_copy_3 )
Pm_copy_3_box <- boxplot_stats( Pm_copy_3 )
Po_copy_3_box <- boxplot_stats( Po_copy_3 )

Pf_copy_4_box <- boxplot_stats( Pf_copy_4 )
Pv_copy_4_box <- boxplot_stats( Pv_copy_4 )
Pm_copy_4_box <- boxplot_stats( Pm_copy_4 )
Po_copy_4_box <- boxplot_stats( Po_copy_4 )



#####################################
##                                 ##
##  PART 3: Asexual Venn diagram   ##
##          	                 ## 
#####################################

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


######################################
##                                  ##
##  PART 4: Gametocyte prevalence   ##
##          	                  ## 
######################################

Pfg_prev_1 <- length( which(ALB_par$pfg_pos==1) )/N_data
Pvg_prev_1 <- length( which(ALB_par$pvg_pos==1) )/N_data
Pmg_prev_1 <- length( which(ALB_par$pmg_pos==1) )/N_data
Pog_prev_1 <- length( which(ALB_par$pog_pos==1) )/N_data


Pfg_prev_2 <- length( intersect( which(ALB_par$pfg_pos==1), 
                                 which(ALB_par$pvg_pos==1 | ALB_par$pmg_pos==1 | ALB_par$pog_pos==1) ) )/N_data
Pvg_prev_2 <- length( intersect( which(ALB_par$pvg_pos==1), 
                                 which(ALB_par$pfg_pos==1 | ALB_par$pmg_pos==1 | ALB_par$pog_pos==1) ) )/N_data
Pmg_prev_2 <- length( intersect( which(ALB_par$pmg_pos==1), 
                                 which(ALB_par$pfg_pos==1 | ALB_par$pvg_pos==1 | ALB_par$pog_pos==1) ) )/N_data
Pog_prev_2 <- length( intersect( which(ALB_par$pog_pos==1), 
                                 which(ALB_par$pfg_pos==1 | ALB_par$pvg_pos==1 | ALB_par$pmg_pos==1) ) )/N_data

Pfg_prev_3 <- length( intersect( which(ALB_par$pfg_pos==1), 
                                 which(rowSums(ALB_par[,c(15,21,25)])>=2) ) )/N_data
Pvg_prev_3 <- length( intersect( which(ALB_par$pvg_pos==1), 
                                 which(rowSums(ALB_par[,c(17,21,25)])>=2) ) )/N_data
Pmg_prev_3 <- length( intersect( which(ALB_par$pmg_pos==1), 
                                 which(rowSums(ALB_par[,c(17,15,25)])>=2) ) )/N_data
Pog_prev_3 <- length( intersect( which(ALB_par$pog_pos==1), 
                                 which(rowSums(ALB_par[,c(17,15,21)])>=2) ) )/N_data

Pfg_prev_4 <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )/N_data
Pvg_prev_4 <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )/N_data
Pmg_prev_4 <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )/N_data
Pog_prev_4 <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )/N_data



Pfg_prev_1_CI <- as.numeric(binom.confint( Pfg_prev_1*N_data, N_data, method="wilson")[1,5:6])
Pvg_prev_1_CI <- as.numeric(binom.confint( Pvg_prev_1*N_data, N_data, method="wilson")[1,5:6])
Pmg_prev_1_CI <- as.numeric(binom.confint( Pmg_prev_1*N_data, N_data, method="wilson")[1,5:6])
Pog_prev_1_CI <- as.numeric(binom.confint( Pog_prev_1*N_data, N_data, method="wilson")[1,5:6])


Pfg_prev_2_CI <- as.numeric(binom.confint( Pfg_prev_2*N_data, N_data, method="wilson")[1,5:6])
Pvg_prev_2_CI <- as.numeric(binom.confint( Pvg_prev_2*N_data, N_data, method="wilson")[1,5:6])
Pmg_prev_2_CI <- as.numeric(binom.confint( Pmg_prev_2*N_data, N_data, method="wilson")[1,5:6])
Pog_prev_2_CI <- as.numeric(binom.confint( Pog_prev_2*N_data, N_data, method="wilson")[1,5:6])


Pfg_prev_3_CI <- as.numeric(binom.confint( Pfg_prev_3*N_data, N_data, method="wilson")[1,5:6])
Pvg_prev_3_CI <- as.numeric(binom.confint( Pvg_prev_3*N_data, N_data, method="wilson")[1,5:6])
Pmg_prev_3_CI <- as.numeric(binom.confint( Pmg_prev_3*N_data, N_data, method="wilson")[1,5:6])
Pog_prev_3_CI <- as.numeric(binom.confint( Pog_prev_3*N_data, N_data, method="wilson")[1,5:6])


Pfg_prev_4_CI <- as.numeric(binom.confint( Pfg_prev_4*N_data, N_data, method="wilson")[1,5:6])
Pvg_prev_4_CI <- as.numeric(binom.confint( Pvg_prev_4*N_data, N_data, method="wilson")[1,5:6])
Pmg_prev_4_CI <- as.numeric(binom.confint( Pmg_prev_4*N_data, N_data, method="wilson")[1,5:6])
Pog_prev_4_CI <- as.numeric(binom.confint( Pog_prev_4*N_data, N_data, method="wilson")[1,5:6])


###################################
##                               ##
##  PART 5: Gametocyte density   ##
##          	               ## 
###################################

Pfg_copy_1 <- ALB_par$pfg_copyno[which(ALB_par$pfg_pos==1)]
Pvg_copy_1 <- ALB_par$pvg_copyno[which(ALB_par$pvg_pos==1)]
Pmg_copy_1 <- ALB_par$pmg_copyno[which(ALB_par$pmg_pos==1)]
Pog_copy_1 <- ALB_par$pog_copyno[which(ALB_par$pog_pos==1)]


Pfg_copy_2 <- ALB_par$pfg_copyno[ intersect( which(ALB_par$pfg_pos==1), 
                                             which(ALB_par$pvg_pos==1 | ALB_par$pmg_pos==1 | ALB_par$pog_pos==1) ) ]
Pvg_copy_2 <- ALB_par$pvg_copyno[ intersect( which(ALB_par$pvg_pos==1), 
                                             which(ALB_par$pfg_pos==1 | ALB_par$pmg_pos==1 | ALB_par$pog_pos==1) ) ]
Pmg_copy_2 <- ALB_par$pmg_copyno[ intersect( which(ALB_par$pmg_pos==1), 
                                             which(ALB_par$pfg_pos==1 | ALB_par$pvg_pos==1 | ALB_par$pog_pos==1) ) ]
Pog_copy_2 <- ALB_par$pog_copyno[ intersect( which(ALB_par$pog_pos==1), 
                                             which(ALB_par$pfg_pos==1 | ALB_par$pvg_pos==1 | ALB_par$pmg_pos==1) ) ]


Pfg_copy_3 <- ALB_par$pfg_copyno[ intersect( which(ALB_par$pfg_pos==1), 
                                             which(rowSums(ALB_par[,c(15,21,25)])>=2) ) ]
Pvg_copy_3 <- ALB_par$pvg_copyno[ intersect( which(ALB_par$pvg_pos==1), 
                                             which(rowSums(ALB_par[,c(17,21,25)])>=2) ) ]
Pmg_copy_3 <- ALB_par$pmg_copyno[ intersect( which(ALB_par$pmg_pos==1), 
                                             which(rowSums(ALB_par[,c(17,15,25)])>=2) ) ]
Pog_copy_3 <- ALB_par$pog_copyno[ intersect( which(ALB_par$pog_pos==1), 
                                             which(rowSums(ALB_par[,c(17,15,21)])>=2) ) ]

Pfg_copy_4 <- ALB_par$pfg_copyno[ which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) ]
Pvg_copy_4 <- ALB_par$pvg_copyno[ which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) ]
Pmg_copy_4 <- ALB_par$pmg_copyno[ which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) ]
Pog_copy_4 <- ALB_par$pog_copyno[ which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) ]


Pfg_copy_1_box <- boxplot_stats( Pfg_copy_1 )
Pvg_copy_1_box <- boxplot_stats( Pvg_copy_1 )
Pmg_copy_1_box <- boxplot_stats( Pmg_copy_1 )
Pog_copy_1_box <- boxplot_stats( Pog_copy_1 )

Pfg_copy_2_box <- boxplot_stats( Pfg_copy_2 )
Pvg_copy_2_box <- boxplot_stats( Pvg_copy_2 )
Pmg_copy_2_box <- boxplot_stats( Pmg_copy_2 )
Pog_copy_2_box <- boxplot_stats( Pog_copy_2 )

Pfg_copy_3_box <- boxplot_stats( Pfg_copy_3 )
Pvg_copy_3_box <- boxplot_stats( Pvg_copy_3 )
Pmg_copy_3_box <- boxplot_stats( Pmg_copy_3 )
Pog_copy_3_box <- boxplot_stats( Pog_copy_3 )

Pfg_copy_4_box <- boxplot_stats( Pfg_copy_4 )
Pvg_copy_4_box <- boxplot_stats( Pvg_copy_4 )
Pmg_copy_4_box <- boxplot_stats( Pmg_copy_4 )
Pog_copy_4_box <- boxplot_stats( Pog_copy_4 )


########################################
##                                    ##
##  PART 6: Gametocyte Venn diagram   ##
##          	                    ## 
########################################

N_g_f <- length( which(ALB_par$pfg_pos==1) )
N_g_v <- length( which(ALB_par$pvg_pos==1) )
N_g_m <- length( which(ALB_par$pmg_pos==1) )
N_g_o <- length( which(ALB_par$pog_pos==1) )

N_g_fv <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1) )
N_g_fm <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1) )
N_g_fo <- length( which(ALB_par$pfg_pos==1 & ALB_par$pog_pos==1) )
N_g_vm <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1) )
N_g_vo <- length( which(ALB_par$pvg_pos==1 & ALB_par$pog_pos==1) )
N_g_mo <- length( which(ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )

N_g_fvm <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1) )
N_g_fvo <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pog_pos==1) )
N_g_fmo <- length( which(ALB_par$pfg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )
N_g_vmo <- length( which(ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )

N_g_fvmo <- length( which(ALB_par$pfg_pos==1 & ALB_par$pvg_pos==1 & ALB_par$pmg_pos==1 & ALB_par$pog_pos==1) )



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

main.size = 2
axis.size = 1
lab.size  = 2
line.size = 2.0
point.size = 0.5


Pf_col <- "red"
Pv_col <- "orange"
Pm_col <- "forestgreen"
Po_col <- "dodgerblue"



tiff( file="Figure1_data_overview_Baseline.tif", width=30, height=20, units="cm", res=500)


lay.mat = rbind( c( 1, 2, 3 ),
                 c( 4, 5, 6 ),
                 c( 7, 7, 7 ) )

layout(lay.mat, heights=c(10,10,1))
layout.show(7)

par(mar = c(3,4.5,2.5,1))

par(mgp = c(3, 0.75, 0))


###################################
##                               ##
##  PANEL 1: Asexual prevalence  ##
##          	               ## 
###################################


line_seq_x <- c(0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.6)



plot(x=100, y=100, 
xlim=c(-0.5,18),  ylim=c(0,0.6), 
xlab="", ylab="prevalence",
main="(A) Asexual parasite prevalence",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)



#############################
## Single

polygon( x=c(0,1,1,0), y=c(0,0,Pf_prev_1,Pf_prev_1),
         col=Pf_col, border=NA )

arrows( x0=0.5, y0=Pf_prev_1_CI[1], x1=0.5, y1=Pf_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=1 + c(0,1,1,0), y=c(0,0,Pv_prev_1,Pv_prev_1),
         col=Pv_col, border=NA )

arrows( x0=1.5, y0=Pv_prev_1_CI[1], x1=1.5, y1=Pv_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=2 + c(0,1,1,0), y=c(0,0,Pm_prev_1,Pm_prev_1),
         col=Pm_col, border=NA )

arrows( x0=2.5, y0=Pm_prev_1_CI[1], x1=2.5, y1=Pm_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=3 + c(0,1,1,0), y=c(0,0,Po_prev_1,Po_prev_1),
         col=Po_col, border=NA )

arrows( x0=3.5, y0=Po_prev_1_CI[1], x1=3.5, y1=Po_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )



#############################
## Double

polygon( x=4.5 + c(0,1,1,0), y=c(0,0,Pf_prev_2,Pf_prev_2),
         col=Pf_col, border=NA )

arrows( x0=5, y0=Pf_prev_2_CI[1], x1=5, y1=Pf_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=5.5 + c(0,1,1,0), y=c(0,0,Pv_prev_2,Pv_prev_2),
         col=Pv_col, border=NA )

arrows( x0=6, y0=Pv_prev_2_CI[1], x1=6, y1=Pv_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=6.5 + c(0,1,1,0), y=c(0,0,Pm_prev_2,Pm_prev_2),
         col=Pm_col, border=NA )

arrows( x0=7, y0=Pm_prev_2_CI[1], x1=7, y1=Pm_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=7.5 + c(0,1,1,0), y=c(0,0,Po_prev_2,Po_prev_2),
         col=Po_col, border=NA )

arrows( x0=8, y0=Po_prev_2_CI[1], x1=8, y1=Po_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )



#############################
## Triple

polygon( x=9 + c(0,1,1,0), y=c(0,0,Pf_prev_3,Pf_prev_3),
         col=Pf_col, border=NA )

arrows( x0=9.5, y0=Pf_prev_3_CI[1], x1=9.5, y1=Pf_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=10 + c(0,1,1,0), y=c(0,0,Pv_prev_3,Pv_prev_3),
         col=Pv_col, border=NA )

arrows( x0=10.5, y0=Pv_prev_3_CI[1], x1=10.5, y1=Pv_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=11 + c(0,1,1,0), y=c(0,0,Pm_prev_3,Pm_prev_3),
         col=Pm_col, border=NA )

arrows( x0=11.5, y0=Pm_prev_3_CI[1], x1=11.5, y1=Pm_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=12 + c(0,1,1,0), y=c(0,0,Po_prev_3,Po_prev_3),
         col=Po_col, border=NA )

arrows( x0=12.5, y0=Po_prev_3_CI[1], x1=12.5, y1=Po_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


#############################
## Quadruple

polygon( x=13.5 + c(0,1,1,0), y=c(0,0,Pf_prev_4,Pf_prev_4),
         col=Pf_col, border=NA  )

arrows( x0=14, y0=Pf_prev_4_CI[1], x1=14, y1=Pf_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=14.5 + c(0,1,1,0), y=c(0,0,Pv_prev_4,Pv_prev_4),
         col=Pv_col, border=NA  )

arrows( x0=15, y0=Pv_prev_4_CI[1], x1=15, y1=Pv_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=15.5 + c(0,1,1,0), y=c(0,0,Pm_prev_4,Pm_prev_4),
         col=Pm_col, border=NA  )

arrows( x0=16, y0=Pm_prev_4_CI[1], x1=16, y1=Pm_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=16.5 + c(0,1,1,0), y=c(0,0,Po_prev_4,Po_prev_4),
         col=Po_col, border=NA  )

arrows( x0=17, y0=Po_prev_4_CI[1], x1=17, y1=Po_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )




for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

axis(1,  at = c(-0.5, 2, 6.5, 11, 15.5, 18), 
     labels = c("", "single", "double", "triple", "quadruple", "" ),
cex.axis=1.5*axis.size )


axis(2, las=2, 
     at = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60), 
     labels = c( "0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%", "45%", "50%", "55%", "60%" ),
cex.axis=axis.size )



###################################
##                               ##
##  PANEL 2: Asexual density     ##
##          	               ## 
###################################


line_seq_x <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6)



plot(x=100, y=100, 
xlim=c(-0.5,18),  ylim=c(0.001, 1e6), log="y",
xlab="", ylab="copy number",
main="(B) Asexual parasite density",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)



#############################
## Single

arrows( x0=0.5, y0=Pf_copy_1_box[[1]][1], x1=0.5, y1=Pf_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=c(0,1,1,0), y=Pf_copy_1_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=c(0,1), y=Pf_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(0.5, length(Pf_copy_1_box[[2]])), y=Pf_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=1.5, y0=Pv_copy_1_box[[1]][1], x1=1.5, y1=Pv_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=1+c(0,1,1,0), y=Pv_copy_1_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=1+c(0,1), y=Pv_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(1.5, length(Pv_copy_1_box[[2]])), y=Pv_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=2.5, y0=Pm_copy_1_box[[1]][1], x1=2.5, y1=Pm_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=2+c(0,1,1,0), y=Pm_copy_1_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=2+c(0,1), y=Pm_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(2.5, length(Pm_copy_1_box[[2]])), y=Pm_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=3.5, y0=Po_copy_1_box[[1]][1], x1=3.5, y1=Po_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=3+c(0,1,1,0), y=Po_copy_1_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=3+c(0,1), y=Po_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(3.5, length(Po_copy_1_box[[2]])), y=Po_copy_1_box[[2]],
        pch=19, cex=point.size )



#############################
## Double

arrows( x0=5, y0=Pf_copy_2_box[[1]][1], x1=5, y1=Pf_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=4.5 + c(0,1,1,0), y=Pf_copy_2_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=4.5 + c(0,1), y=Pf_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(5, length(Pf_copy_2_box[[2]])), y=Pf_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=6, y0=Pv_copy_2_box[[1]][1], x1=6, y1=Pv_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=5.5+c(0,1,1,0), y=Pv_copy_2_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=5.5+c(0,1), y=Pv_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(6, length(Pv_copy_2_box[[2]])), y=Pv_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=7, y0=Pm_copy_2_box[[1]][1], x1=7, y1=Pm_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=6.5+c(0,1,1,0), y=Pm_copy_2_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=6.5+c(0,1), y=Pm_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(7, length(Pm_copy_2_box[[2]])), y=Pm_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=8, y0=Po_copy_2_box[[1]][1], x1=8, y1=Po_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=7.5+c(0,1,1,0), y=Po_copy_2_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=7.5+c(0,1), y=Po_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(8, length(Po_copy_2_box[[2]])), y=Po_copy_2_box[[2]],
        pch=19, cex=point.size )



#############################
## Triple

arrows( x0=9.5, y0=Pf_copy_3_box[[1]][1], x1=9.5, y1=Pf_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=9 + c(0,1,1,0), y=Pf_copy_3_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=9 + c(0,1), y=Pf_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(9.5, length(Pf_copy_3_box[[2]])), y=Pf_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=10.5, y0=Pv_copy_3_box[[1]][1], x1=10.5, y1=Pv_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=10+c(0,1,1,0), y=Pv_copy_3_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=10+c(0,1), y=Pv_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(10.5, length(Pv_copy_3_box[[2]])), y=Pv_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=11.5, y0=Pm_copy_3_box[[1]][1], x1=11.5, y1=Pm_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=11+c(0,1,1,0), y=Pm_copy_3_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=11+c(0,1), y=Pm_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(11.5, length(Pm_copy_3_box[[2]])), y=Pm_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=12.5, y0=Po_copy_3_box[[1]][1], x1=12.5, y1=Po_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=12+c(0,1,1,0), y=Po_copy_3_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=12+c(0,1), y=Po_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(12.5, length(Po_copy_3_box[[2]])), y=Po_copy_3_box[[2]],
        pch=19, cex=point.size )




#############################
## Quadruple

arrows( x0=14, y0=Pf_copy_4_box[[1]][1], x1=14, y1=Pf_copy_4_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=13.5 + c(0,1,1,0), y=Pf_copy_4_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=13.5 + c(0,1), y=Pf_copy_4_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(14, length(Pf_copy_4_box[[2]])), y=Pf_copy_4_box[[2]],
        pch=19, cex=point.size )



arrows( x0=15, y0=Pv_copy_4_box[[1]][1], x1=15, y1=Pv_copy_4_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=14.5+c(0,1,1,0), y=Pv_copy_4_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=14.5+c(0,1), y=Pv_copy_4_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(15, length(Pv_copy_4_box[[2]])), y=Pv_copy_4_box[[2]],
        pch=19, cex=point.size )



arrows( x0=16, y0=Pm_copy_4_box[[1]][1], x1=16, y1=Pm_copy_4_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=15.5+c(0,1,1,0), y=Pm_copy_4_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=15.5+c(0,1), y=Pm_copy_4_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(16, length(Pm_copy_4_box[[2]])), y=Pm_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=17, y0=Po_copy_4_box[[1]][1], x1=17, y1=Po_copy_4_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=16.5+c(0,1,1,0), y=Po_copy_4_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=16.5+c(0,1), y=Po_copy_4_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(17, length(Po_copy_4_box[[2]])), y=Po_copy_3_box[[2]],
        pch=19, cex=point.size )



for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c(-0.5, 2, 6.5, 11, 15.5, 18), 
     labels = c("", "single", "double", "triple", "quadruple", "" ),
cex.axis=1.5*axis.size )


axis(2, las=2, 
     at     = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e8), 
     labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e8),
cex.axis=axis.size )




#####################################
##                                 ##
##  PANEL 3: Asexual Venn diagram  ##
##          	                 ## 
#####################################

##plot.new()

##draw.quad.venn( area1=N_f, area2=N_v, area3=N_m, area4=N_o, 
##                n12=N_fv, n13=N_fm, n14=N_fo, n23=N_vm n24=N_vo, n34=N_mo, 
##                n123=N_fvm, n124=N_fvo, n134=N_fmo, n234=N_vmo, n1234=N_fvmo,
##                col=c(Pf_col, Pv_col, Pm_col, Po_col),
##                fill=c(Pf_col, Pv_col, Pm_col, Po_col) ) 



plot(x=100, y=100, 
xlim=c(0,1),
xlab="", ylab="",
main="(C) Asexual parasite co-infection",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)





######################################
##                                  ##
##  PANEL 4: Gametocyte prevalence  ##
##          	                  ## 
######################################


line_seq_x <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)



plot(x=100, y=100, 
xlim=c(-0.5,18),  ylim=c(0,0.4), 
xlab="", ylab="prevalence",
main="(D) Gametocyte prevalence",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)



#############################
## Single

polygon( x=c(0,1,1,0), y=c(0,0,Pfg_prev_1,Pfg_prev_1),
         col=Pf_col, border=NA )

arrows( x0=0.5, y0=Pfg_prev_1_CI[1], x1=0.5, y1=Pfg_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=1 + c(0,1,1,0), y=c(0,0,Pvg_prev_1,Pvg_prev_1),
         col=Pv_col, border=NA )

arrows( x0=1.5, y0=Pvg_prev_1_CI[1], x1=1.5, y1=Pvg_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=2 + c(0,1,1,0), y=c(0,0,Pmg_prev_1,Pmg_prev_1),
         col=Pm_col, border=NA )

arrows( x0=2.5, y0=Pmg_prev_1_CI[1], x1=2.5, y1=Pmg_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=3 + c(0,1,1,0), y=c(0,0,Pog_prev_1,Pog_prev_1),
         col=Po_col, border=NA )

arrows( x0=3.5, y0=Pog_prev_1_CI[1], x1=3.5, y1=Pog_prev_1_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )



#############################
## Double

polygon( x=4.5 + c(0,1,1,0), y=c(0,0,Pfg_prev_2,Pfg_prev_2),
         col=Pf_col, border=NA )

arrows( x0=5, y0=Pfg_prev_2_CI[1], x1=5, y1=Pfg_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=5.5 + c(0,1,1,0), y=c(0,0,Pvg_prev_2,Pvg_prev_2),
         col=Pv_col, border=NA )

arrows( x0=6, y0=Pvg_prev_2_CI[1], x1=6, y1=Pvg_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=6.5 + c(0,1,1,0), y=c(0,0,Pmg_prev_2,Pmg_prev_2),
         col=Pm_col, border=NA )

arrows( x0=7, y0=Pmg_prev_2_CI[1], x1=7, y1=Pmg_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=7.5 + c(0,1,1,0), y=c(0,0,Pog_prev_2,Pog_prev_2),
         col=Po_col, border=NA )

arrows( x0=8, y0=Pog_prev_2_CI[1], x1=8, y1=Pog_prev_2_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )



#############################
## Triple

polygon( x=9 + c(0,1,1,0), y=c(0,0,Pfg_prev_3,Pfg_prev_3),
         col=Pf_col, border=NA )

arrows( x0=9.5, y0=Pfg_prev_3_CI[1], x1=9.5, y1=Pfg_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=10 + c(0,1,1,0), y=c(0,0,Pvg_prev_3,Pvg_prev_3),
         col=Pv_col, border=NA )

arrows( x0=10.5, y0=Pvg_prev_3_CI[1], x1=10.5, y1=Pvg_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=11 + c(0,1,1,0), y=c(0,0,Pmg_prev_3,Pmg_prev_3),
         col=Pm_col, border=NA )

arrows( x0=11.5, y0=Pmg_prev_3_CI[1], x1=11.5, y1=Pmg_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=12 + c(0,1,1,0), y=c(0,0,Pog_prev_3,Pog_prev_3),
         col=Po_col, border=NA )

arrows( x0=12.5, y0=Pog_prev_3_CI[1], x1=12.5, y1=Pog_prev_3_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


#############################
## Quadruple

polygon( x=13.5 + c(0,1,1,0), y=c(0,0,Pfg_prev_4,Pfg_prev_4),
         col=Pf_col, border=NA  )

arrows( x0=14, y0=Pfg_prev_4_CI[1], x1=14, y1=Pfg_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=14.5 + c(0,1,1,0), y=c(0,0,Pvg_prev_4,Pvg_prev_4),
         col=Pv_col, border=NA  )

arrows( x0=15, y0=Pvg_prev_4_CI[1], x1=15, y1=Pvg_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=15.5 + c(0,1,1,0), y=c(0,0,Pmg_prev_4,Pmg_prev_4),
         col=Pm_col, border=NA  )

arrows( x0=16, y0=Pmg_prev_4_CI[1], x1=16, y1=Pmg_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )


polygon( x=16.5 + c(0,1,1,0), y=c(0,0,Pog_prev_4,Pog_prev_4),
         col=Po_col, border=NA  )

arrows( x0=17, y0=Pog_prev_4_CI[1], x1=17, y1=Pog_prev_4_CI[2],
        angle=90, code=3, length=0.04, lty="solid"  )






for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}



axis(1,  at = c(-0.5, 2, 6.5, 11, 15.5, 18), 
     labels = c("", "single", "double", "triple", "quadruple", "" ),
cex.axis=1.5*axis.size )


axis(2, las=2, 
     at = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40), 
     labels = c( "0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%" ),
cex.axis=axis.size )



###################################
##                               ##
##  PANEL 5: Gametocyte density  ##
##          	               ## 
###################################


line_seq_x <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4)


plot(x=100, y=100, 
xlim=c(-0.5,18),  ylim=c(0.001, 10000), log="y",
xlab="", ylab="copy number",
main="(E) Gametocyte density",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)



#############################
## Single

arrows( x0=0.5, y0=Pfg_copy_1_box[[1]][1], x1=0.5, y1=Pfg_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=c(0,1,1,0), y=Pfg_copy_1_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=c(0,1), y=Pfg_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(0.5, length(Pfg_copy_1_box[[2]])), y=Pfg_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=1.5, y0=Pvg_copy_1_box[[1]][1], x1=1.5, y1=Pvg_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=1+c(0,1,1,0), y=Pvg_copy_1_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=1+c(0,1), y=Pvg_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(1.5, length(Pvg_copy_1_box[[2]])), y=Pvg_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=2.5, y0=Pmg_copy_1_box[[1]][1], x1=2.5, y1=Pmg_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=2+c(0,1,1,0), y=Pmg_copy_1_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=2+c(0,1), y=Pmg_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(2.5, length(Pmg_copy_1_box[[2]])), y=Pmg_copy_1_box[[2]],
        pch=19, cex=point.size )



arrows( x0=3.5, y0=Pog_copy_1_box[[1]][1], x1=3.5, y1=Pog_copy_1_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=3+c(0,1,1,0), y=Pog_copy_1_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=3+c(0,1), y=Pog_copy_1_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(3.5, length(Pog_copy_1_box[[2]])), y=Pog_copy_1_box[[2]],
        pch=19, cex=point.size )



#############################
## Double

arrows( x0=5, y0=Pfg_copy_2_box[[1]][1], x1=5, y1=Pfg_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=4.5 + c(0,1,1,0), y=Pfg_copy_2_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=4.5 + c(0,1), y=Pfg_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(5, length(Pfg_copy_2_box[[2]])), y=Pfg_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=6, y0=Pvg_copy_2_box[[1]][1], x1=6, y1=Pvg_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=5.5+c(0,1,1,0), y=Pvg_copy_2_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=5.5+c(0,1), y=Pvg_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(6, length(Pvg_copy_2_box[[2]])), y=Pvg_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=7, y0=Pmg_copy_2_box[[1]][1], x1=7, y1=Pmg_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=6.5+c(0,1,1,0), y=Pmg_copy_2_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=6.5+c(0,1), y=Pmg_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(7, length(Pmg_copy_2_box[[2]])), y=Pmg_copy_2_box[[2]],
        pch=19, cex=point.size )



arrows( x0=8, y0=Pog_copy_2_box[[1]][1], x1=8, y1=Pog_copy_2_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=7.5+c(0,1,1,0), y=Pog_copy_2_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=7.5+c(0,1), y=Pog_copy_2_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(8, length(Pog_copy_2_box[[2]])), y=Pog_copy_2_box[[2]],
        pch=19, cex=point.size )



#############################
## Triple

arrows( x0=9.5, y0=Pfg_copy_3_box[[1]][1], x1=9.5, y1=Pfg_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=9 + c(0,1,1,0), y=Pfg_copy_3_box[[1]][c(2,2,4,4)],
         col=Pf_col, border=NA )

points( x=9 + c(0,1), y=Pfg_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(9.5, length(Pfg_copy_3_box[[2]])), y=Pfg_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=10.5, y0=Pvg_copy_3_box[[1]][1], x1=10.5, y1=Pvg_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=10+c(0,1,1,0), y=Pvg_copy_3_box[[1]][c(2,2,4,4)],
         col=Pv_col, border=NA )

points( x=10+c(0,1), y=Pvg_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(10.5, length(Pvg_copy_3_box[[2]])), y=Pvg_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=11.5, y0=Pmg_copy_3_box[[1]][1], x1=11.5, y1=Pmg_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=11+c(0,1,1,0), y=Pmg_copy_3_box[[1]][c(2,2,4,4)],
         col=Pm_col, border=NA )

points( x=11+c(0,1), y=Pmg_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(11.5, length(Pmg_copy_3_box[[2]])), y=Pmg_copy_3_box[[2]],
        pch=19, cex=point.size )



arrows( x0=12.5, y0=Pog_copy_3_box[[1]][1], x1=12.5, y1=Pog_copy_3_box[[1]][5],
        angle=90, code=3, length=0.04, lty="solid"  )

polygon( x=12+c(0,1,1,0), y=Pog_copy_3_box[[1]][c(2,2,4,4)],
         col=Po_col, border=NA )

points( x=12+c(0,1), y=Pog_copy_3_box[[1]][c(3,3)], 
        type='l', lwd=line.size )

points( x=rep(12.5, length(Pog_copy_3_box[[2]])), y=Pog_copy_3_box[[2]],
        pch=19, cex=point.size )




#############################
## Quadruple

##arrows( x0=14, y0=Pfg_copy_4_box[[1]][1], x1=14, y1=Pfg_copy_4_box[[1]][5],
##        angle=90, code=3, length=0.04, lty="solid"  )
##
##polygon( x=13.5 + c(0,1,1,0), y=Pfg_copy_4_box[[1]][c(2,2,4,4)],
##         col=Pf_col, border=NA )
##
##points( x=13.5 + c(0,1), y=Pfg_copy_4_box[[1]][c(3,3)], 
##        type='l', lwd=line.size )
##
##points( x=rep(14, length(Pfg_copy_4_box[[2]])), y=Pfg_copy_4_box[[2]],
        pch=19, cex=point.size )

points( x=14, y=Pfg_copy_4_box[[1]][1],
        pch=19, cex=point.size )



##arrows( x0=15, y0=Pvg_copy_4_box[[1]][1], x1=15, y1=Pvg_copy_4_box[[1]][5],
##        angle=90, code=3, length=0.04, lty="solid"  )
##
##polygon( x=14.5+c(0,1,1,0), y=Pvg_copy_4_box[[1]][c(2,2,4,4)],
##         col=Pv_col, border=NA )
##
##points( x=14.5+c(0,1), y=Pvg_copy_4_box[[1]][c(3,3)], 
##        type='l', lwd=line.size )
##
##points( x=rep(15, length(Pvg_copy_4_box[[2]])), y=Pvg_copy_4_box[[2]],
##        pch=19, cex=point.size )

points( x=15, y=Pvg_copy_4_box[[1]][1],
        pch=19, cex=point.size )



##arrows( x0=16, y0=Pmg_copy_4_box[[1]][1], x1=16, y1=Pmg_copy_4_box[[1]][5],
##        angle=90, code=3, length=0.04, lty="solid"  )
##
##polygon( x=15.5+c(0,1,1,0), y=Pmg_copy_4_box[[1]][c(2,2,4,4)],
##         col=Pm_col, border=NA )
##
##points( x=15.5+c(0,1), y=Pmg_copy_4_box[[1]][c(3,3)], 
##        type='l', lwd=line.size )
##
##points( x=rep(16, length(Pmg_copy_4_box[[2]])), y=Pmg_copy_4_box[[2]],
##        pch=19, cex=point.size )

points( x=16, y=Pmg_copy_4_box[[1]][1],
        pch=19, cex=point.size )



##arrows( x0=17, y0=Pog_copy_4_box[[1]][1], x1=17, y1=Pog_copy_4_box[[1]][5],
##        angle=90, code=3, length=0.04, lty="solid"  )
##
##polygon( x=16.5+c(0,1,1,0), y=Pog_copy_4_box[[1]][c(2,2,4,4)],
##         col=Po_col, border=NA )
##
##points( x=16.5+c(0,1), y=Pog_copy_4_box[[1]][c(3,3)], 
##        type='l', lwd=line.size )
##
##points( x=rep(17, length(Pog_copy_4_box[[2]])), y=Pog_copy_4_box[[2]],
##        pch=19, cex=point.size )

points( x=17, y=Pog_copy_4_box[[1]][1],
        pch=19, cex=point.size )






for(i in 1:length(line_seq_x))
{
	points(x=c(-1e10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}



axis(1,  at = c(-0.5, 2, 6.5, 11, 15.5, 18), 
     labels = c("", "single", "double", "triple", "quadruple", "" ),
cex.axis=1.5*axis.size )


axis(2, las=2, 
     at     = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4), 
     labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4),
cex.axis=axis.size )




########################################
##                                    ##
##  PANEL 6: Gametocyte Venn diagram  ##
##          	                    ## 
########################################



plot(x=1000, y=1000, 
xlim=c(0,1),
xlab="", ylab="",
main="(F) Gametocyte co-infection",
xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty='n',
cex.main=main.size, cex.lab=lab.size)




##draw.quad.venn( area1=N_g_f, area2=N_g_v, area3=N_g_m, area4=N_g_o, 
##                n12=N_g_fv, n13=N_g_fm, n14=N_g_fo, n23=N_g_vm n24=N_g_vo, n34=N_g_mo, 
##                n123=N_g_fvm, n124=N_g_fvo, n134=N_g_fmo, n234=N_g_vmo, n1234=N_g_fvmo,
##                col=c(Pf_col, Pv_col, Pm_col, Po_col),
##                fill=c(Pf_col, Pv_col, Pm_col, Po_col) ) 



########################################
##                                    ##
##  PANEL 7: Legend                   ##
##          	                    ## 
########################################

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("P. falciparum", "P. vivax", "P. malariae", "P. ovale"),
       fill = c(Pf_col, Pv_col, Pm_col, Po_col), 
       border = c(Pf_col, Pv_col, Pm_col, Po_col), 
       ncol=4, cex=2.0, bty="n" )



dev.off()






###################################
###################################
##                               ##
##  ##   ## ##### #   ## #   ##  ## 
##  ##   ## ##    ##  ## ##  ##  ## 
##   ## ##  ####  ### ## ### ##  ##
##    ###   ##    ## ### ## ###  ##
##     #    ##### ##  ## ##  ##  ##
##                               ##
###################################
###################################



tiff( file="Figure1_Asexual_Venn_Baseline.tif", width=20, height=20, units="cm", res=500)



par(mfrow=c(1,1))

draw.quad.venn( area1=N_f, area2=N_v, area3=N_m, area4=N_o, 
                n12=N_fv, n13=N_fm, n14=N_fo, n23=N_vm, n24=N_vo, n34=N_mo, 
                n123=N_fvm, n124=N_fvo, n134=N_fmo, n234=N_vmo, n1234=N_fvmo,
                col=c(Pf_col, Pv_col, Pm_col, Po_col),
                fill=c(Pf_col, Pv_col, Pm_col, Po_col),
                lwd=4, alpha=0.3,
                cex=2, fontfamily="sans" ) 




dev.off()




tiff( file="Figure1_Gametocyte_Venn_Baseline.tif", width=20, height=20, units="cm", res=500)



par(mfrow=c(1,1))

draw.quad.venn( area1=N_g_f, area2=N_g_v, area3=N_g_m, area4=N_g_o, 
                n12=N_g_fv, n13=N_g_fm, n14=N_g_fo, n23=N_g_vm, n24=N_g_vo, n34=N_g_mo, 
                n123=N_g_fvm, n124=N_g_fvo, n134=N_g_fmo, n234=N_g_vmo, n1234=N_g_fvmo,
                col=c(Pf_col, Pv_col, Pm_col, Po_col),
                fill=c(Pf_col, Pv_col, Pm_col, Po_col),
                lwd=4, alpha=0.3,
                cex=2, fontfamily="sans" ) 




dev.off()




