#----------------------------------------------------------------------
#Supplementary material
#Genome-Wide Association Studies in Multi-Environment Trials
#______________________________________________________________________

#----------------------------------------------------------------------
#M1: Single-environment model involving pedigree information
#______________________________________________________________________

# To load library
library(sommer)

#Database
DT<-read.table("DT.txt", sep=" ", header = T, dec=".")

# Subset by environment
Y1<- DT[1:1797,]
Y2<- DT[1798:3594,]
Y3<- DT[3595:5391,]
Y4<- DT[5392:7188,]

# Molecular Marker Matrix
M<- DT[1:599,3:1281]
rownames(M) <- unique(as.factor(DT$Name))
M

# Pedigree Matrix (Kinship)
K<- read.table("pedigree.txt", sep="\t", header = T, dec=",") 
rownames(K)<-unique(as.factor(DT$Name))
colnames(K)<-unique(as.factor(DT$Name))
K <- as.matrix(K)

#-------------------------------------------------------------------------
#Note:
#For the M1 model in the remaining environments, only change 
#the data in the GWAS function. That is: Y2, Y3, or Y4.
#-------------------------------------------------------------------------


# Model M1 (Environment 1) 

M11 <- GWAS(Yield~1,
            random= ~ vs(Name, Gu=K),method = "NR",
            rcov= ~ units, M=M, gTerm = "u:Name",
            data=Y1)

#Variance components
M11$sigma
# Genotype's variance
Vg_M11<- round(M11$sigma[[1]],3)
# Residual variance
Ve_M11<- round(M11$sigma[[2]],3)
#AIC 
M11$AIC
#BIC
M11$BIC
#Variance component's Standard Error
SE_M11<- round(sqrt(diag(M11$sigmaSE)),4)
#Genotype's BLUPs
BLUP_M11<- randef(M11) 
#Scores
scores_M11<- as.data.frame(round(M11$scores[2,1:1109],4))
Markers<-row.names(scores_M11)
scores_M11 <- data.frame(Markers,scores_M11)
colnames(scores_M11)<- c("Markers","Scores_M11")
scores_M11
write.table(scores_M11,"Scores_M11.txt",row.names=FALSE,sep="\t")

#------------------------------------------------------------------------
#M2: Single-environment model including molecular similarity
#________________________________________________________________________

# Additive relationship matrix
A <- A.mat(M)
A

#------------------------------------------------------------------------
#Note:
#For the M2 model in the remaining environments, only change 
#the data in the GWAS function. That is: Y2, Y3, or Y4.
#------------------------------------------------------------------------


# Model M2 (Environment 1) 

M21 <- GWAS(Yield~1,
            random= ~ vs(Name, Gu=A),method = "NR",
            rcov= ~ units, M=M, gTerm = "u:Name",
            data=Y1)

#Variance components
M21$sigma
# Genotype's variance
Vg_M21<- round(M21$sigma[[1]],3)
# Residual variance
Ve_M21<- round(M21$sigma[[2]],3)
#AIC 
M21$AIC
#BIC
M21$BIC
#Variance component's Standard Error
SE_M21<- round(sqrt(diag(M21$sigmaSE)),4)
#Genotype's BLUPs
BLUP_M21<- randef(M21) 
#Scores
scores_M21<- as.data.frame(round(M21$scores[2,1:1109],4))
Markers<-row.names(scores_M21)
scores_M21 <- data.frame(Markers,scores_M21)
colnames(scores_M21)<- c("Markers","scores_M21")
scores_M21
write.table(scores_M21,"Scores_M21.txt",row.names=FALSE,sep="\t")


#----------------------------------------------------------------------
#M3: Multi-environment model involving pedigree information
#______________________________________________________________________


# Number of genotypes
nind <- length(unique(DT$Name))
# Number of environments
nenv <- length(unique(DT$Env))

# Correlation matrix among environments 
Rho<- diag(nenv)
colnames(Rho) <- rownames(Rho) <- unique(DT$Env)
Rho

# Kronecker matrix using the Rho and K matrices
Rho_K<- kronecker(Rho,K,make.dimnames=TRUE)
Rho_K

# Model M3

M3 <- GWAS(Yield~Env,
           random= ~ vs(Name, Gu=K) + vs(Env:Name, Gu=Rho_K),
           rcov= ~ units, M=M, gTerm = "u:Name",method = "NR",
           data=DT)

#Variance components
M3$sigma
# Genotype's variance
Vg_M3<- round(M3$sigma[[1]],3)
# Genotype-environment's variance
Vge_M3<- round(M3$sigma[[2]],3)
# Residual variance
Ve_M3<- round(M3$sigma[[3]],3)
#AIC 
M3$AIC
#BIC
M3$BIC
#Variance component's Standard Error
SE_M3<- round(sqrt(diag(M3$sigmaSE)),4)
#Genotype's BLUPs
BLUP_G_M3<- randef(M3)[[1]] 
# GxE Interaction's BLUPs 
BLUP_GE_M3<- randef(M3)[[2]]
#Scores
scores_M3<- as.data.frame(round(M3$scores[2,1:1109],4))
Markers<-row.names(scores_M3)
scores_M3 <- data.frame(Markers,scores_M3)
colnames(scores_M3)<- c("Markers","Scores_M3")
scores_M3
write.table(scores_M3,"Scores_M3.txt",row.names=FALSE,sep="\t")

#----------------------------------------------------------------------
#M4: Multi-environment model including molecular similarity
#______________________________________________________________________

# Kronecker matrix using the Rho and A matrices
Rho_A<- kronecker(Rho,A,make.dimnames=TRUE)
Rho_A

# # Model M4

M4 <- GWAS(Yield~Env,
           random= ~ vs(Name, Gu=A) + vs(Env:Name, Gu=Rho_A),
           rcov= ~ units, M=M, gTerm = "u:Name",method = "NR",
           data=DT)

#Variance components
M4$sigma
round(M4$sigma,3)
# Genotype's variance
Vg_M4<- round(M4$sigma[[1]],3)
# Genotype-environment's variance
Vge_M4<- round(M4$sigma[[2]],3)
# Residual variance
Ve_M4<- round(M4$sigma[[3]],3)
#AIC 
M4$AIC
#BIC
M4$BIC
#Variance component's Standard Error
SE_M4<- round(sqrt(diag(M4$sigmaSE)),4)
#Genotype's BLUPs
BLUP_G_M4<- randef(M4)[[1]] 
# GxE Interaction's BLUPs 
BLUP_GE_M4<- randef(M4)[[2]]
#Scores
scores_M4<- as.data.frame(round(M4$scores[2,1:1109],4))
Markers<-row.names(scores_M4)
scores_M4 <- data.frame(Markers,scores_M4)
colnames(scores_M4)<- c("Markers","Scores_M4")
scores_M4
write.table(scores_M4,"Scores_M4.txt",row.names=FALSE,sep="\t")