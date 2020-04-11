# ***
# 0) load packages
#

library(foreign)
library(mice)
library(miceadds)
library(mitml)

# ***
# 1) prepare data
#

# student data set (level 1)
stu.dat <- read.spss("PiSA2012_SPSS_student_DEU.sav", to.data.frame=TRUE)

# school data set (level 2)
sch.dat <- read.spss("PiSA2012_SPSS_school_DEU.sav", to.data.frame=TRUE)

# create combined data sets (level 1 and 2)
ind <- match(stu.dat$SCHOOLID, sch.dat$SCHOOLID)
stu.vars <- c("SCHOOLID","ST04Q01","ESCS","CLSMAN","STUDREL","PV1MATH")
sch.vars <- c("SC11Q01","SC11Q02")

# combine
dat <- cbind(stu.dat[,stu.vars], sch.dat[ind,sch.vars])

# reformat and rename
dat <- within(dat, SCHOOLID <- as.numeric(SCHOOLID))
colnames(dat) <- c("SCHOOLID","GENDER","ESCS","CLSMAN","STUDREL","PVMATH","NSTU","NCMP")

# fix seed for random number generator(s)
seed <- 1234

# ***
# 2.i) joint modeling (JM)
#

fml <- list( GENDER + ESCS + CLSMAN + STUDREL + PVMATH ~ 1 + (1|SCHOOLID) ,
             NSTU + NCMP ~ 1 )

imp.jm <- mitml::jomoImpute(data=dat, formula=fml, n.burn=5000, n.iter=500, m=20,
                            seed=seed)

# completed data sets
for(i in 1:20){
  out <- within(mitml::mitmlComplete(imp.jm,i),{ GENDER <- as.numeric(GENDER)-1 })
  write.table(out, file=paste0("jm",i,".dat"), row.names=FALSE, col.names=FALSE)
  cat(paste0("jm",i,".dat\n"), file="jm_list.dat", append=TRUE)
}

# ***
# 2.ii) single-level FCS (FCS-SL)
#

predMatrix <- matrix(0, ncol(dat), ncol(dat))
rownames(predMatrix) <- colnames(predMatrix) <- colnames(dat)

predMatrix["ESCS",] <- c(0,1,0,1,1,1,1,1)      # regression (flat file)
predMatrix["CLSMAN",] <- c(0,1,1,0,1,1,1,1)    # regression (flat file)
predMatrix["STUDREL",] <- c(0,1,1,1,0,1,1,1)   # regression (flat file)
predMatrix["NSTU",] <- c(0,1,1,1,1,1,0,1)      # regression (flat file)
predMatrix["NCMP",] <- c(0,1,1,1,1,1,1,0)      # regression (flat file)

impMethod <- c("","","norm","norm","norm","","norm","norm")

imp.fcssl <- mice::mice(data=dat, predictorMatrix=predMatrix, method=impMethod,
                        m=20, maxit=50, seed=seed)

# completed data sets
for(i in 1:20){
  out <- within(mice::complete(imp.fcssl,i),{ GENDER <- as.numeric(GENDER)-1
    NSTU <- mitml::clusterMeans(NSTU,SCHOOLID)
    NCMP <- mitml::clusterMeans(NCMP,SCHOOLID)
  })
  write.table(out, file=paste0("fcssl",i,".dat"), row.names=FALSE, col.names=FALSE)
  cat(paste0("fcssl",i,".dat\n"), file="fcssl_list.dat", append=TRUE)
}

# ***
# 2.iii) two-level FCS with manifest cluster means (FCS-MAN)
#

predMatrix <- matrix(0, ncol(dat), ncol(dat))
rownames(predMatrix) <- colnames(predMatrix) <- colnames(dat)

predMatrix["ESCS",] <- c(-2,1,0,3,3,3,1,1)      # random intercepts model
predMatrix["CLSMAN",] <- c(-2,1,3,0,3,3,1,1)    # random intercepts model
predMatrix["STUDREL",] <- c(-2,1,3,3,0,3,1,1)   # random intercepts model
predMatrix["NSTU",] <- c(-2,1,1,1,1,1,0,1)      # regression at level 2
predMatrix["NCMP",] <- c(-2,1,1,1,1,1,1,0)      # regression at level 2

impMethod <- c("","","2l.pan","2l.pan","2l.pan","","2lonly.norm","2lonly.norm")

imp.fcsman <- mice::mice(data=dat, predictorMatrix=predMatrix, method=impMethod,
                         m=20, maxit=50, seed=seed)

# completed data sets
for(i in 1:20){
  out <- within(mice::complete(imp.fcsman,i),{ GENDER <- as.numeric(GENDER)-1 })
  write.table(out, file=paste0("fcsman",i,".dat"), row.names=FALSE, col.names=FALSE)
  cat(paste0("fcsman",i,".dat\n"), file="fcsman_list.dat", append=TRUE)
}

# ***
# 2.iv) two-level FCS with manifest and latent cluster means (FCS-LAT)
#

lm.dat <- within(dat, LM.CLSMAN <- LM.STUDREL <- LM.PVMATH <- NA)

predMatrix <- matrix(0, ncol(lm.dat), ncol(lm.dat))
rownames(predMatrix) <- colnames(predMatrix) <- colnames(lm.dat)

predMatrix["LM.CLSMAN",] <- c(-2,1,1,2,0,0,1,1,0,1,1)    # latent group means
predMatrix["LM.STUDREL",] <- c(-2,1,1,0,2,0,1,1,1,0,1)   # latent group means
predMatrix["LM.PVMATH",] <- c(-2,1,1,0,0,2,1,1,1,1,0)    # latent group means
predMatrix["ESCS",] <- c(-2,1,0,1,1,1,1,1,1,1,1)         # random intercepts model
predMatrix["CLSMAN",] <- c(-2,1,3,0,1,1,1,1,0,1,1)       # random intercepts model
predMatrix["STUDREL",] <- c(-2,1,3,1,0,1,1,1,1,0,1)      # random intercepts model
predMatrix["NSTU",] <- c(-2,1,1,0,0,0,0,1,1,1,1)         # regression at level 2
predMatrix["NCMP",] <- c(-2,1,1,0,0,0,1,0,1,1,1)         # regression at level 2

visitSeq <- c(11,4,9,5,10,3,7,8)

impMethod <- c("","","2l.pan","2l.pan","2l.pan","","2lonly.norm","2lonly.norm",
               "2l.latentgroupmean.mcmc","2l.latentgroupmean.mcmc",
               "2l.latentgroupmean.mcmc")

imp.fcslat <- mice::mice(data=lm.dat, predictorMatrix=predMatrix, method=impMethod,
                         visitSequence=visitSeq, allow.na=TRUE, m=20, maxit=50,
                         seed=seed)

# completed data sets
for(i in 1:20){
  out <- within(mice::complete(imp.fcslat,i),{ GENDER <- as.numeric(GENDER)-1 })
  write.table(out, file=paste0("fcslat",i,".dat"), row.names=FALSE, col.names=FALSE)
  cat(paste0("fcslat",i,".dat\n"), file="fcslat_list.dat", append=TRUE)
}
