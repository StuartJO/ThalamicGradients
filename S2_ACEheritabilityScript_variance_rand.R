#!/Library/Frameworks/R.framework/Resources/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

require(OpenMx)   # not needed for the simulation, but will be needed when we come to model specification
require(MASS) 
require(R.matlab) 

library(R.matlab)
library(MASS)
library(OpenMx)

set.seed(200)

# load the data
setwd("/projects/kg98/stuarto/SeedReg/Heritability")
data_covar = readMat(args[1])
data = readMat(args[2])


# replace years with centuries (recommended by openMx administrators)
data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

# replace missing covariates with a very different number setting them to "pseudo-missing" values 
# as openMx doesn't allow missing covariate values (recommended by openMx administrators)
data_covar$MZ.ID[is.na(data_covar$MZ.ID)] <- -999
data_covar$MZ.age[is.na(data_covar$MZ.age)] <- -999
data_covar$MZ.sex[is.na(data_covar$MZ.sex)] <- -999

data_covar$DZ.ID[is.na(data_covar$DZ.ID)] <- -999
data_covar$DZ.age[is.na(data_covar$DZ.age)] <- -999
data_covar$DZ.sex[is.na(data_covar$DZ.sex)] <- -999

# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityT <- numeric(numEdges)
heritabilityE <- numeric(numEdges)
heritabilityS <- numeric(numEdges)
heritabilitySp <- numeric(numEdges)
heritabilityVARA <- numeric(numEdges)
heritabilityVARC <- numeric(numEdges)
heritabilityVART <- numeric(numEdges)
heritabilityVARE <- numeric(numEdges)
heritabilitySubRemDZ <- numeric(numEdges)
heritabilitySubRemMZ <- numeric(numEdges)

# (recommended by openMx administrators)
mxOption(NULL,"Default optimizer","SLSQP")

for (edge in c(1:numEdges)){

# replaces outlier with NA
outliersValueMZ<- boxplot.stats(data$Output.MZ[,,edge])$out
data$Output.MZ[,,edge][data$Output.MZ[,,edge] %in% outliersValueMZ]=NaN;
  
outliersValueDZ<- boxplot.stats(data$Output.DZ[,,edge])$out
data$Output.DZ[,,edge][data$Output.DZ[,,edge] %in% outliersValueDZ]=NaN;

heritabilitySubRemDZ[edge] = length(outliersValueDZ);
heritabilitySubRemMZ[edge] = length(outliersValueMZ);
  
# each connection (phenotype of interest) is denoted as edge
myDataMZ<-data.frame(data$Output.MZ[,,edge], data_covar$MZ.age[,1], data_covar$MZ.age[,2],data_covar$MZ.sex[,1],data_covar$MZ.sex[,2])
myDataDZ<-data.frame(data$Output.DZ[,,edge], data_covar$DZ.age[,1], data_covar$DZ.age[,2],data_covar$DZ.sex[,1],data_covar$DZ.sex[,2])

myDataMZ_measure<-data$Output.MZ[,,edge]
myDataDZ_measure<-data$Output.DZ[,,edge]

colnames(myDataMZ) <- c('twin1', 'twin2','ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(myDataDZ) <- c('twin1', 'twin2','ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')
selVars <- c('twin1','twin2')

# "complete" specifies use only cases with data in all columns, otherwise subjects would be excluded if they're missing an additional sibling.
CovMZ = cov(data$Output.MZ[,,edge],use="complete")
CovDZ = cov(data$Output.DZ[,,edge],use="complete")

# mean across MZ twins
MeanMZ = colMeans(data$Output.MZ[,1:2,edge],na.rm=TRUE)
MeanMZ = mean(MeanMZ)
# mean across DZ twins
MeanDZ = colMeans(data$Output.DZ[,1:2,edge],na.rm=TRUE)
MeanDZ = mean(MeanDZ)

mylabels=c("twin1","twin2")
MZsat <- mxModel("MZsat",
                 
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanMZ,MeanMZ), labels =c("b0_mz1","b0_mz2"), name="Intercepts" ),
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type = "Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholMZ" ),
                 mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                 mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel("DZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanDZ,MeanDZ), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type = "Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholDZ" ),
                 mxAlgebra( CholDZ %*% t(CholDZ),name="expCovDZ"),
                 mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

# Model specification starts here
SatModel <- mxModel("twinSat", MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
#---------------------------------------------------------------------------------------------------------------------------
SatModelFit <- mxTryHard(SatModel)
summary(SatModelFit)

# -----------------------------------------------------------------------
#Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
twinACE <- mxModel("twinACE",
                   # Matrices X, Y, and Z to store a, c, and e path coefficients
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="a", name="X" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="c", name="Y" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="e", name="Z" ),
                   # Matrices A, C, and E compute variance components
                   mxAlgebra( expression=X %*% t(X), name="A" ),
                   mxAlgebra( expression=Y %*% t(Y), name="C" ),
                   mxAlgebra( expression=Z %*% t(Z), name="E" ),

                   # Algebra for expected variance/covariance matrix in MZ+sibling
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , A+C),
                                         cbind(A+C  , A+C+E)),
                     name="expCovMZ"),

                   # Algebra for expected variance/covariance matrix in DZ+siblig
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , 0.5%x%A+C),
                                         cbind(0.5%x%A+C , A+C+E)),
                     name="expCovDZ"),

                   mxModel("MZ", mxData( observed=myDataMZ, type="raw" ),
                           # Algebra for making the means a function of the definition variables age and sex
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, c(MeanMZ,MeanMZ), labels =c("b0_mz1","b0_mz2"), name="Intercepts"),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                           mxExpectationNormal( covariance="twinACE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
                           mxFitFunctionML()),

                   mxModel("DZ", mxData( observed=myDataDZ, type="raw" ),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, c(MeanDZ,MeanDZ), labels=c("b0_dz1","b0_dz2"), name="Intercepts"),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                           mxExpectationNormal( covariance="twinACE.expCovDZ", means="expMeanDZ", dimnames=selVars ),
                           mxFitFunctionML()),

                   mxFitFunctionMultigroup( c("MZ.fitfunction",  "DZ.fitfunction"))
)
twinACEFit<-mxTryHard(twinACE)

#Run ACE model
# -----------------------------------------------------------------------
#estCovMZ  <- mxEval(twinACE.expCovMZ, twinACEFit)      # expected covariance matrix for MZ's
#estCovDZ  <- mxEval(twinACE.expCovDZ, twinACEFit)      # expected covariance matrix for DZ's
estVA1     <- mxEval(a*a, twinACEFit)              # additive genetic variance, a^2
estVC1     <- mxEval(c*c, twinACEFit)              # shared enviromnemtal variance, c^2
estVE1     <- mxEval(e*e, twinACEFit)              # unique environmental variance, e^2
estVP1     <- (estVA1+estVC1+estVE1)                  # total variance

heritmodels <- matrix(nrow = 4, ncol = 9)

heritmodels[1,1] <- estVA1/estVP1                          # standardized additive genetic variance
heritmodels[1,2] <- estVC1/estVP1
heritmodels[1,3] <- 0
heritmodels[1,4] <- estVE1/estVP1
heritmodels[1,5] <- twinACEFit@output$status$code
heritmodels[1,6] <- estVA1
heritmodels[1,7] <- estVC1
heritmodels[1,8] <- 0
heritmodels[1,9] <- estVE1

# Generate AE Model - C=0
twinAE      <- twinACE
twinAE      <- mxRename(twinAE, "twinAE")
twinAE      <- omxSetParameters(twinAE, labels="c", free=FALSE, values=0 )
twinAEFit   <- mxTryHard(twinAE)
AESumm      <- summary(twinAEFit)

estVA3    <- mxEval(a*a, twinAEFit)              # additive genetic variance, a^2
estVC3    <- mxEval(c*c, twinAEFit)              # shared enviromnemtal variance, c^2
estVE3    <- mxEval(e*e, twinAEFit)              # unique environmental variance, e^2
estVP3    <- (estVA3+estVC3+estVE3)                  # total variance

heritmodels[2,1] <- estVA3/estVP3                          # standardized additive genetic variance
heritmodels[2,2] <- estVC3/estVP3
heritmodels[2,3] <- 0
heritmodels[2,4] <- estVE3/estVP3
heritmodels[2,5] <- twinAEFit@output$status$code
heritmodels[2,6] <- estVA3
heritmodels[2,7] <- estVC3
heritmodels[2,8] <- 0
heritmodels[2,9] <- estVE3

#Generate CE Model - A=0
twinCE      <- twinACE
twinCE      <- mxRename(twinCE, "twinCE")
twinCE      <- omxSetParameters(twinCE, labels="a", free=FALSE, values=0 )
twinCEFit   <- mxTryHard(twinCE)

estVA4    <- mxEval(a*a, twinCEFit)              # additive genetic variance, a^2
estVC4    <- mxEval(c*c, twinCEFit)              # shared enviromnemtal variance, c^2
estVE4    <- mxEval(e*e, twinCEFit)              # unique environmental variance, e^2
estVP4    <- (estVA4+estVC4+estVE4)              # total variance

heritmodels[3,1] <- estVA4/estVP4                          # standardized additive genetic variance
heritmodels[3,2] <- estVC4/estVP4
heritmodels[3,3] <- 0
heritmodels[3,4] <- estVE4/estVP4
heritmodels[3,5] <- twinCEFit@output$status$code
heritmodels[3,6] <- estVA4
heritmodels[3,7] <- estVC4
heritmodels[3,8] <- 0
heritmodels[3,9] <- estVE4

#Generate E Model - A & C=0
twinE      <- twinCE
twinE      <- mxRename(twinE, "twinE")
twinE      <- omxSetParameters(twinE, labels="c", free=FALSE, values=0 )
twinEFit   <- mxTryHard(twinE)

estVA5    <- mxEval(a*a, twinEFit)              # additive genetic variance, a^2
estVC5    <- mxEval(c*c, twinEFit)   
estVE5    <- mxEval(e*e, twinEFit)              # unique environmental variance, e^2
estVP5    <- (estVA5+estVC5+estVE5)              # total variance

heritmodels[4,1] <- estVA5/estVP5                          # standardized additive genetic variance
heritmodels[4,2] <- estVC5/estVP5
heritmodels[4,3] <- 0
heritmodels[4,4] <- estVE5/estVP5
heritmodels[4,5] <- twinEFit@output$status$code
heritmodels[4,6] <- estVA5
heritmodels[4,7] <- estVC5
heritmodels[4,8] <- 0
heritmodels[4,9] <- estVE5


# model comparison
options('digits' = 5)
# compare saturated to other models
compValuesACE = mxCompare(SatModelFit, c(twinACEFit,twinAEFit,twinCEFit,twinEFit))
#compValuesACE = mxCompare(twinACEFit, c(twinAEFit,twinCEFit,twinEFit))

# find model with the lowest AIC
INDmin = which.min(compValuesACE$AIC[1:5])
AICmin = compValuesACE$AIC[INDmin]

# check if this model is not significantly different from saturated
pmin = compValuesACE$p[INDmin]

if (INDmin==1) {
  INDmin2 = which.min(compValuesACE$AIC[2:6])
  heritabilityA[edge] <- heritmodels[INDmin2,1]
  heritabilityC[edge] <- heritmodels[INDmin2,2]
  heritabilityT[edge] <- heritmodels[INDmin2,2]
  heritabilityE[edge] <- heritmodels[INDmin2,4]
  heritabilityS[edge] <- 1; # 1 if saturated model gave lower AIC value
  heritabilitySp[edge] <- compValuesACE$p[INDmin2+1] # selecting from a table where first row is excluded
  heritabilityVARA[edge] <- heritmodels[INDmin2,6]
  heritabilityVARC[edge] <- heritmodels[INDmin2,7]
  heritabilityVART[edge] <- heritmodels[INDmin2,8]
  heritabilityVARE[edge] <- heritmodels[INDmin2,9]

} else {
  heritabilityA[edge] <- heritmodels[INDmin-1,1]
  heritabilityC[edge] <- heritmodels[INDmin-1,2]
  heritabilityT[edge] <- heritmodels[INDmin-1,3]
  heritabilityE[edge] <- heritmodels[INDmin-1,4]
  heritabilityS[edge] <- heritmodels[INDmin-1,5]
  heritabilitySp[edge] <- compValuesACE$p[INDmin] # selecting from a table where values in first row are NaNs
  heritabilityVARA[edge] <- heritmodels[INDmin-1,6]
  heritabilityVARC[edge] <- heritmodels[INDmin-1,7]
  heritabilityVART[edge] <- heritmodels[INDmin-1,8]
  heritabilityVARE[edge] <- heritmodels[INDmin-1,9]
}
}

heritabilityACE <- data.frame(heritabilityA,heritabilityC,heritabilityT,heritabilityE,heritabilityS,heritabilitySp,heritabilityVARA,heritabilityVARC,heritabilityVART,heritabilityVARE)
setwd("/projects/kg98/stuarto/SeedReg/Heritability")
write.csv(heritabilityACE,sprintf("heritabilityACE_%s.txt",args[2]),row.names=FALSE)


