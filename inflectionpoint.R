library(readODS)
library(tidyverse)
library(tibble)
library(reshape2)
library(inflection)
library(glue)

#table with AF and ANI values
df<- read_ods("/home/mcoin/Petrotoga/ANI/aniaf_flexpoint_all_better.ods", row_names=FALSE)
x_val<-df$ANIavg
y_val<-df$Afavg
#fit a polynomial of 4th order to the data
polynome<-lm(y_val ~ x_val + I(x_val^2)+ I(x_val^3)+ I(x_val^4))
summary(polynome)
#create a function with the fitted parameters
cfs<-polynome$coefficients
ff=function(x){cfs[1]+cfs[2]*x+cfs[3]*x^2+cfs[4]*x^3+cfs[5]*x^4}
fs<-paste0(cfs[1],'+0',cfs[2],'*x+0',cfs[3],'*x^2+0',cfs[4],'*x^3+0',cfs[5],'*x^4+0', collapse='')
f<-parse(text = fs)
f
#calculate derivatives
firstd<-D(f , "x")
secondd<-D(firstd,  "x")
secondd
sdf<-function(x){eval(secondd)}
#find where the second derivative is zero
x_lim<-uniroot(sdf,c(65, 80))$root
uniroot(sdf,c(80, 100))
y_lim<-ff(x_lim)
#values of ANI and AF of D. tunisiensis, the hard boundary
x_tun<-73.00
y_tun<-0.54
