# meta-analysis-trainning
This is a biref introduction of all the meta-analysis universe in Rstudio.

First let's install R, and then: Rstudio 
1)	First install R: https://cran.r-project.org/mirrors.html
2)	Now you can install Rstudio: https://www.rstudio.com/products/rstudio/download/
One aditional video explain in detail all the steps for the instalation: https://www.youtube.com/watch?v=NZxSA80lF1I

The packages that we will use are:
install.packages("metafor")
install.packages("meta")
install.packages("pimeta")
install.packages("bamdit")
install.packages("metaDigitise")
install.packages("revtools")
#in case some other package is missing
require("metafor")
require("meta")
require("pimeta")
require("bamdit")
require("metaDigitise")
require("revtools")



#EXPLORING DIffferent  MODELS in meta
####metabin
####
# Use subset of Olkin (1995) to conduct meta-analysis based on
# inverse variance method (with risk ratio as summary measure)
#
library(meta)
data(Olkin1995)
Olkin1995.sub<- Olkin1995[c(41, 47, 51, 59),]
m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
              data=Olkin1995.sub,
              method = "Inverse")
m3 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
              data = Olkin1995, subset = year < 1970,
              sm = "OR", method = "Inverse", studlab = author)
m4 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
              data = Olkin1995, subset = year < 1970,
              sm = "OR",  method.tau = "REML",
              hakn = TRUE, studlab = author)
m5 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
              data = Olkin1995, subset = year < 1970,
              sm = "OR",  method = "GLMM",model.glmm ="CM.AL",
              hakn = TRUE, studlab = author)

forest(m3, layout = "RevMan5")
forest(m1, layout = "RevMan5")
forest(m4, layout = "RevMan5")
forest(m5, layout = "RevMan5")



#####PIMETA + meta
library(meta)  #metagen
library(pimeta) #pima
library(dmetar)
data(ThirdWave)
?dmetar::ThirdWave
head(ThirdWave)
#  generic inverse variance
#

m.gen <- metagen(TE = TE,
                 seTE = seTE,
                 studlab = Author,
                 data = ThirdWave,
                 sm = "SMD",
                 comb.fixed = FALSE,
                 comb.random = TRUE,
                 prediction = TRUE,
                 method.tau = "REML",
                 hakn = TRUE,
                 title = "Third Wave Psychotherapies")

pidat <- pima(ThirdWave$TE, ThirdWave$seTE, parallel = 6, seed = 123456)
m.gen$lower.predict <- pidat$lpi
m.gen$upper.predict <- pidat$upi
summary(m.gen)
forest(m.gen)


#Ok Let's do it simpler with metafor




#DTA MA in METAFOR 
# Epxloring the Profile likelihood

rm(list=ls())
library(metafor)
library(metadat)
dat <- dat.kearon1998

#id	numeric	study id
#author	character	study author(s)
#year	numeric	publication year
#patients	character	patient group (either symptomatic or asymptomatic patients)
#tp	numeric	number of true positives
#np	numeric	number of positive patients (cases)
#tn	numeric	number of true negatives
#nn	numeric	number of negative patients (non-cases)

dat <- to.long(measure="OR", ai=tp, n1i=np, ci=tn, n2i=nn,
               data=dat.kearon1998)
levels(dat$group) <- c("sensitivity", "specificity")
### calculate logit-transformed sensitivities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="sensitivity")
### calculate logit-transformed specificities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="specificity")
### bivariate random-effects model for logit sensitivity and specificity
res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | study, struct="UN", data=dat)
res
forest(res)
profile.rma.mv(res)




library(bamdit)
data("glas")
glas.t <- glas[glas$marker == "Telomerase", 1:4]


glas.t$tn<-glas.t$n2-glas.t$fp
dat <- to.long(measure="OR", ai=tp, n1i=n1, ci=tn, n2i=n2, data=glas.t)
levels(dat$group) <- c("sensitivity", "specificity")
### calculate logit-transformed sensitivities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="sensitivity")
### calculate logit-transformed specificities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="specificity")
### bivariate random-effects model for logit sensitivity and specificity
res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | study, struct="UN", data=dat)
res
forest(res)
profile.rma.mv(res)




# Example: data from Glas et al. (2003).....................................
library(bamdit)
data("glas")
glas.t <- glas[glas$marker == "Telomerase", 1:4]

# Simple visualization ...
plotdata(glas.t,                # Data frame
         two.by.two = FALSE     # Data is given as: (tp, n1, fp, n2)
)

glas.m1 <- metadiag(glas.t,                # Data frame
                    two.by.two = FALSE,    # Data is given as: (tp, n1, fp, n2)
                    re = "normal",         # Random effects distribution
                    re.model = "DS",       # Random effects on D and S
                    link = "logit",        # Link function
                    sd.Fisher.rho   = 1.7, # Prior standard deviation of correlation
                    nr.burnin = 1000,      # Iterations for burnin
                    nr.iterations = 10000, # Total iterations
                    nr.chains = 2,         # Number of chains
                    r2jags = TRUE)         # Use r2jags as interface to jags

summary(glas.m1, digit=3)
bamdit::plotsesp(glas.m1)

plot(glas.m1,                    # Fitted model
     level = c(0.5), # Credibility levels
     parametric.smooth = TRUE)   # Parametric curve
plot(glas.m1,                    # Fitted model
     level = c(0.5, 0.75, 0.95), # Credibility levels
     parametric.smooth = TRUE)   # Parametric curve




library(metaDigitise)
library(revtools)
