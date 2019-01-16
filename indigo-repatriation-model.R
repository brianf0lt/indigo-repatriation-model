#######################################################################################
############## Modeling repatriation scenarios for Eastern Indigo Snakes ##############
############## 		B Folt, CP McGowan, DA Steen, C Guyer				 ##############
##############						12-11-2018							 ##############
#######################################################################################

### Clear the workspace
rm(list=ls())

install.packages("statmod", dependencies=TRUE)
library(statmod)

############## Part 1) Modeling demographic parameters for
############## Eastern Indigo Snakes (Drymarchon couperi)

# Model female populations of D. couperi as a population with five distinct life stages:
# 1) hatchlings/first-year snakes (ages 0-1)
# 2) juveniles (ca. 1-2 yr) 
# 3) subadults (ca. 2-3 yr)
# 4) first-year adults (ca. 3-4 yr)
# 5) adults (ca. 4<= yr)

# The survival within each life stage and transition between stages can be conceptualized with a population transition matrix:
# 
#		| 0			0				0				Sa1*Fa1*Ba1		Sa*Fa*Ba	|
#		| Sh*Thj	Sj(1-Tjsa)		0				0				0		|
# A =	| 0			Sj*Tjsa		Ssa(1-Tsaa1)			0	    			0		|, where:
#		| 0			0				Ssa*Tsaa1		Sa1(1-Ta1a)		0		|
#		| 0 		0				0				Sa1*Ta1a			Sa		|
#	
# S is the survival rate at stage i, where i can be h=hatchling, j=juvenile, sa=subadult, a1=primiparous adult, a=adult
# F is the fecundity at stages a1 and a, and 
# B is the likelihood of breeding at stages a1 and a.

### To simulate these parameters, we will use 1000 replications (r) that are projected over 30 years (t)
r = 1000
t = 30


######## Define and parameterize each of the above variables here using realistical values estimated from the literature: Hyslop et al. 2012 Population Ecology

### Sh -- survival of hatchlings
## This stage is comparable to the hatchling/juvenile stage from Hyslop et al. 2012, where they estimated survival in the first year by two components: hatchling survival during the first three months (0.49) and juvenile survival during the subsequent nine months (0.59). We modeled this as a conservative estimate of the two; 0.52
mSh = 0.52  		
varSh = 0.1							# Variance of mean
aSh = mSh*((mSh*(1-mSh)/(varSh^2))-1) 	 	
bSh = (1-mSh)*((mSh*(1-mSh)/(varSh^2))-1) 
Shi = matrix(rbeta(r,aSh,bSh),r,1)	# Parametric uncertainty 
SDmShi = matrix(rinvgauss(r,varSh^2,1),r,1)
AShi = matrix(0,r,1)				# beta distribution shape parameters
BShi = matrix(0,r,1)				# beta distribution shape parameters
Sht = matrix(0,r,t)					# annual variation survival


### Sj -- survival of juveniles 
## This stage is somewhat comparable to the subadult stage from Hyslop et al. 2012, where they estimated survival as 0.52 (0.20 SE). However, their sample size was low for individuals in this stage, and it seems unlikely that second year individuals would have lower survival than the hatchling stage. So we modeled survival in the second year (0.60) as a modest increase relative to survival in the first year.
mSj = 0.60  		
varSj = 0.1		# Variance of mean
aSj = mSj*((mSj*(1-mSj)/(varSj^2))-1) 	 	
bSj = (1-mSj)*((mSj*(1-mSj)/(varSj^2))-1) 
Sji = matrix(rbeta(r,aSj,bSj),r,1)			# Parametric uncertainty 
SDmSji = matrix(rinvgauss(r,varSj^2,1),r,1)
ASji = matrix(0,r,1)				# beta distribution shape parameters
BSji = matrix(0,r,1)				# beta distribution shape parameters
Sjt = matrix(0,r,t)					# annual variation survival

### Ssa -- survival of subadults  
## This stage roughly translates to the first-year adults from Hyslop et al. 2012's adult survival stage (0.74), so we estimated this as 0.70 as a slightly reduced estimate of that.
mSsa = 0.70  		
varSsa = 0.1		# Variance of mean 
aSsa = mSsa*((mSsa*(1-mSsa)/(varSsa^2))-1) 	 	
bSsa = (1-mSsa)*((mSsa*(1-mSsa)/(varSsa^2))-1) 
Ssai = matrix(rbeta(r,aSsa,bSsa),r,1)			# Parametric uncertainty 
SDmSsai = matrix(rinvgauss(r,varSsa^2,1),r,1)
ASsai = matrix(0,r,1)				# beta distribution shape parameters
BSsai = matrix(0,r,1)				# beta distribution shape parameters
Ssat = matrix(0,r,t)					# annual variation survival

### Sa1 -- survival of first-year adults 
## This stage roughly translates to the second-year adults from Hyslop et al. 2012's adult survival stage (0.74), so we estimated this as 0.70 as an increased estimate of that (0.80)
mSa1 = 0.80	
varSa1 = 0.1		# Variance of mean 
aSa1 = mSa1*((mSa1*(1-mSa1)/(varSa1^2))-1) 	 	
bSa1 = (1-mSa1)*((mSa1*(1-mSa1)/(varSa1^2))-1) 
Sa1i = matrix(rbeta(r,aSa1,bSa1),r,1)			# Parametric uncertainty 
SDmSa1i = matrix(rinvgauss(r,varSa1^2,1),r,1)
ASa1i = matrix(0,r,1)				# beta distribution shape parameters
BSa1i = matrix(0,r,1)				# beta distribution shape parameters
Sa1t = matrix(0,r,t)					# annual variation survival

### Sa -- survival of adults 
## Radiotelemetered individuals in CNF that survive one year proceed to have high 
## survival (0.90) (Stiles 2013). 
mSa = 0.85 		
varSa = 0.1		# Variance of mean
aSa = mSa*((mSa*(1-mSa)/(varSa^2))-1) 	 	
bSa = (1-mSa)*((mSa*(1-mSa)/(varSa^2))-1) 
Sai = matrix(rbeta(r,aSa,bSa),r,1)			# Parametric uncertainty 
SDmSai = matrix(rinvgauss(r,varSa^2,1),r,1)
ASai = matrix(0,r,1)				# beta distribution shape parameters
BSai = matrix(0,r,1)				# beta distribution shape parameters
Sat = matrix(0,r,t)					# annual variation survival

### The first four stages should theoretically transition into the next age stage after one year. However, there may be uncertainty if individuals vary in growth and individuals remain in stages, especially as they age. To model this uncertainty, we modeled transition probability declining from 0.99 with each successive stage up to the primiparous adult stage, until when primiparous females all transition to be adults.
mThj = 0.99  			# Mean transition from hatchling to juvenile
varThj = 0.01		
aThj = mSh*((mSh*(1-mSh)/varSh^2)-1)
bThj = (1-mSh)*((mSh*(1-mSh)/varSh^2)-1)
Thj = matrix(rbeta(r*t,aThj,bThj),r,t)

mTjsa = 0.90 			# Mean transition from juvenile to subadult
varTjsa = 0.03	
aTjsa = mSj*((mSj*(1-mSj)/varSj^2)-1)
bTjsa = (1-mSj)*((mSj*(1-mSj)/varSj^2)-1)
Tjsa = matrix(rbeta(r*t,aTjsa,bTjsa),r,t)

mTsaa1 = 0.80  			# Mean transition from subadult to primiparous adult
varTsaa1 = 0.05	
aTsaa1 = mSsa*((mSsa*(1-mSsa)/varSsa^2)-1)
bTsaa1 = (1-mSsa)*((mSsa*(1-mSsa)/varSsa^2)-1)
Tsaa1 = matrix(rbeta(r*t,aTsaa1,bTsaa1),r,t)

mTa1a = 0.99  			# Mean transition from primiparous adult to adult
varTa1a = 0.01	
aTa1a = mSa1*((mSa1*(1-mSa1)/varSa1^2)-1)
bTa1a = (1-mSa1)*((mSa1*(1-mSa1)/varSa1^2)-1)
Ta1a = matrix(rbeta(r*t,aTa1a,bTa1a),r,t)

### Fa -- fecundity of adults (clutch size)
## Mean Fa = 9 (4-12, range; Hyslop et al. 2012)
## Fa = 8.65 (6-12, range; husbandry data, C. Guyer pers. comm.)
muFa = 8.65
sdFa = 2
fecundShape2Fa = log((sdFa^2)/(muFa^2)+1)
fecundShape1Fa = log(muFa)-1/2*fecundShape2Fa
#hist(round(rlnorm(100,fecundShape1Fa,fecundShape2Fa)))	# for visualization purposes
round(rlnorm(1,fecundShape1Fa,fecundShape2Fa))  # Clutch size for matrix
eggsa = matrix(0,r,t)				# eggs from adults

### Pva -- proportion of viable eggs for adults
Pva = 0.85 
aPva = 100*Pva
bPva = 100*(1-Pva)
Pvat=matrix(rbeta(r*t,aPva,bPva),r,t)

### Fa1 -- fecundity (clutch size) of primiparous adults
## A primary reason why we modeled two stages of adults is because, in our experience, first-year breeding females lay clutches dominated by inviable eggs. So, we wanted to model these females as having different egg viability relative to older, more experienced individuals; clutch size was the same. 
muFa1 = 8.65							
sdFa1 = 2								
fecundShape2Fa1 = log((sdFa1^2)/(muFa1^2)+1)
fecundShape1Fa1 = log(muFa1)-1/2*fecundShape2Fa1
# hist(round(rlnorm(100,fecundShape1Fa1,fecundShape2Fa1)))	# for visualization purposes
round(rlnorm(1,fecundShape1Fa1,fecundShape2Fa1))  # Clutch size for matrix
eggsa1 = matrix(0,r,t)				# eggs from primiparous adults (a1)

### Pva1 -- proportion of viable eggs for primiparous adults
Pva1 = 0.35 
aPva1 = 100*Pva1
bPva1 = 100*(1-Pva1)
Pva1t=matrix(rbeta(r*t,aPva1,bPva1),r,t)

### Se -- survival of eggs (nests); following an estimate by Hyslop et al. (2012), but no data from the field have estimated this parameter
Se = 0.75 
varSe = 0.15
aSe=Se*((Se*(1-Se)/varSe^2)-1)
bSe=(1-Se)*((Se*(1-Se)/varSe^2)-1)
Set=matrix(0,r,t)

### SR -- Sex ratio of eggs
## A beta-distributed variable 0.5 (+/-0.04) to only model the proportion of females assuming a 1:1 sex ratio in clutches
mSR = 0.5
varSR = 0.04
aSR=mSR*((mSR*(1-mSR)/(varSR^2))-1)
bSR=(1-mSR)*((mSR*(1-mSR)/(varSR^2))-1)
SRi = matrix(rbeta(r,aSR,bSR),r,1)			 
SDmSRi = matrix(rinvgauss(r,aSR^2,1),r,1)
ASRi = matrix(0,r,1)				
BSRi = matrix(0,r,1)				
SRt = matrix(0,r,t)				

### Pbt -- proportion of individuals that breed
## An attempt to simulate some good-year/bad-year dynamics.
Pbt = matrix(0,r,t) 
PGY = 0.8
GY = matrix(rbinom(r*t,1,PGY),r,t) 
SDPb=0.2

### Proportion of individuals that breed during any given year.
PB = 0.8
varPB = 0.1
aPB = PB*((PB*(1-PB)/varPB^2)-1)
bPB = (1-PB)*((PB*(1-PB)/varPB^2)-1)
rbeta(1,aPB,bPB)

### Sr -- survival reduction for captives released into wild
## Individuals being released into the wild may not perform as well as wild-born individuals, so we modeled an acclimation effect on the survival of captive-bred snakes being released, where survival during the released life stage was multiplied by this variable. 
Sr = 0.5 						
Srt=matrix(0,r,t)
	
### pSamp -- the probability of sampling snakes
## Snakes are difficult to detect, especially fossorial species that have large homeranges. This parameter models how detection probability influences our ability to monitor population growth.
mPsamp=0.20 
SDmPsamp = (mPsamp*0.1)*(mPsamp*0.1)
aPs= mPsamp*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
bPs=(1-mPsamp)*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
Psamp = matrix(rbeta(r*t,aPs,bPs),r,t)
# We model pSamp = 0.20 here, but in the model for-loops we explore other values 


# Simulating parameters is now done
detach(package:statmod)



############## Part 2)
############## Project the population growth and extinction risk of D. couperi  
############## under different repatriation scenarios 

t = 30		# Number of years to project simulations
r = 1000	# Number of simulation replications 

### Population growth & quasi-extinction estimates
rlam=matrix(0,r,t)		# Real (true) population growth
olam=matrix(,r,t)		# Observed population growth, after sampling (detection) process
Pext=matrix(0,r,t)		# Probability of quasi-extinction

#### Preliminary modeling exercises created a repatriation strategy where 30 'head-started', 2 yr-old subadult snakes would be released each year for 10 years. This sought to create a sustainable population with relatively low risk of extinction. However, this model lacked parametric uncertainty, it is expensive to raise snakes through the second year, and also difficult to raise enough snakes to release that many for 10 years. So, managers were also interested in the feasability of releasing younger snakes, such as hatchlings or 1-yr old juveniles, and whether shorter release periods (e.g., 5 years) or irregular release periods (7 out of 10 years)  would be sufficient to generate viable populations with low risk of extinction. However, two problems potentially limiting success are that (1) snakes are difficult to sample and may have detection probabilities much lower than 1.0, and (2) released animals may have decreased survival in the first year (post-release acclimation effect) while they try to learn and figure out how to live in the wild.  

#### This provides us five factors to consider in our modeling exercise: 
## (1) annual release size (15 or 30 individuals)
## (2) snake release age (juvenile or subadult)
## (3) release duration (5 yr, 10 yr)
## (4) probability of sampling (observing) snakes (0.05, 0.15, 0.25, 0.35)
## (5) acclimation effects on post-release survival

#### 16 release scenarios

cols = c("nrelsa","nrelj","nyr","mPsamp","Sr")
rows = c("A","B","C","D","E","F","G","H","M","N","O","P","CNFc","ABRPc","CNFp","ABRPp")

scenarios = matrix(c(
	# Multiple snake release scenarios with different values for
	# subadult releases (col 1), juvenile releases (col 2), release duration in years (col 3), 
	# probability of releases (col 4), and detection probability (col 5)
	30,0,10,0.2,0.5,	#scenario A - evaluate different release number, age, and length
	0,30,10,0.2,0.5,	#scenario B 
	30,0,5,0.2,0.5,		#scenario C 
	0,30,5,0.2,0.5,		#scenario D
	15,0,10,0.2,0.5,	#scenario E 
	0,15,10,0.2,0.5,	#scenario F 
	15,0,5,0.2,0.5,		#scenario G 
	0,15,5,0.2,0.5,		#scenario H 
	
	30,0,10,0.35,0.8,		#scenario M -- evaluate different detection probabilities
	30,0,10,0.25,0.6,		#scenario N 		and acclimation effects on survival (Sr)
	30,0,10,0.15,0.4,		#scenario O
	30,0,10,0.05,0.2,		#scenario P
	
	0,0,0,0.15,0.5,		#scenario CNF	-- current population projections for CNF
	0,0,0,0.15,0.5,		#scenario ABRP	-- current population projections for ABRP	
	0,0,0,0.15,0.5,		#scenario CNF	-- projected population projections for CNF
	0,0,0,0.15,0.5),	#scenario ABRP	-- projected population projections for ABRP
	nrow=16, ncol=5, byrow=TRUE, dimnames=list(rows,cols))
					
n = length(scenarios[,1])

### Create an empty matrix to dump all the results into
results = matrix(0, nrow = 1, ncol = 32, dimnames=list(NA,c("Scenario","Stage",1:30)))


### Use a for-loop to iteratively calculate demography under different management scenarios
### by using each row of the scenarios object to provide unique combinations of parameters
 
for (h in 1:n){			# For-loop for each repatriation scenario
			
nrelsa=scenarios[h,1] 	# Mean no. of 2-yr old subadults releases (both sexes)
nrelj=scenarios[h,2]	# Mean number of 1 yr-old juveniles releases (both sexes)
nyr=scenarios[h,3] 		# Duration of the release program (years)

mPsamp=scenarios[h,4] 	# Probabilty of sampling snakes (i.e., detection probability)
SDmPsamp = (mPsamp*0.1)*(mPsamp*0.1)
aPs= mPsamp*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
bPs=(1-mPsamp)*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
Psamp = matrix(rbeta(r*t,aPs,bPs),r,t)

Sr = scenarios[h,5]		# Acclimation effect on survival (i.e., survival reduction) of released individuals	
Srt=matrix(0,r,t)

### Create some matrices and parameters for the model loop:
Nh =  matrix(0,r,t)		# Abundance of hatchlings 
Nj = matrix(0,r,t) 		# Juveniles
Nsa = matrix(0,r,t) 	# Subadults
Na1 = matrix(0,r,t) 	# Primiparous adults
Na = matrix(0,r,t) 		# Adults

Nh[,1]=0		# Initial abundance of hatchlings
Nj[,1]=0		# Juveniles
Nsa[,1]=0		# Subadults
Na1[,1]=0		# Primiparous adults
Na[,1]=0		# Adults

Njr = matrix(0,r,t) 	# Number of captive-reared juveniles (1-yr olds) released
Nsar = matrix(0,r,t) 	# Number of captive-reared subadults (2-yr olds) released

# Specify number of individuals currently released TO DATE in Conecuh and Apalachicola to date
rep.row<-function(x,n){matrix(rep(x,each=n),nrow=n)}	# function to build matrix

CNFc = rep.row(c(0,17*rbeta(1,aSR,bSR),30*rbeta(1,aSR,bSR),31*rbeta(1,aSR,bSR),20*rbeta(1,aSR,bSR),0,9*rbeta(1,aSR,bSR),0,26*rbeta(1,aSR,bSR),24*rbeta(1,aSR,bSR),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),1000)

ABRPc = rep.row(c(0,12*rbeta(1,aSR,bSR),20*rbeta(1,aSR,bSR),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),1000)

# Specify number of individuals PROJECTED for release in CNF and ABRP
CNFp = rep.row(c(0,17*rbeta(1,aSR,bSR),30*rbeta(1,aSR,bSR),31*rbeta(1,aSR,bSR),20*rbeta(1,aSR,bSR),0,9*rbeta(1,aSR,bSR),0,26*rbeta(1,aSR,bSR),24*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),23*rbeta(1,aSR,bSR),0,0,0,0,0,0,0,0,0,0,0,0,0,0),1000)

ABRPp = rep.row(c(0,12*rbeta(1,aSR,bSR),20*rbeta(1,aSR,bSR),27*rbeta(1,aSR,bSR),26*rbeta(1,aSR,bSR), 27*rbeta(1,aSR,bSR),26*rbeta(1,aSR,bSR),27*rbeta(1,aSR,bSR),26*rbeta(1,aSR,bSR),27*rbeta(1,aSR,bSR),26*rbeta(1,aSR,bSR),27*rbeta(1,aSR,bSR),26*rbeta(1,aSR,bSR),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),1000)

#this script simulates the same release schedule for all replicates and 
#is susceptible to stochastic variation

Nobreeders = matrix(0,r,t) 	# Number of adults observed/counted, given Psamp
Noimmatures = matrix(0,r,t)	# Number of immatures observed/counted, given Psamp

Nmax=100 			# Population ceiling for density dependence


#### Population-projection model that accounts for imperfect detection

for(i in 1:r){			# Replication loop; draws replicate-level means for survival at each stage

ASai[i] = 100*Sai[i]			# Adults 		
BSai[i] = 100*(1-Sai[i])		

ASa1i[i] = 100*Sa1i[i]		 	# Primiparous adults	
BSa1i[i] = 100*(1-Sa1i[i])

ASsai[i] = 100*Ssai[i]		 	# Subadults	
BSsai[i] = 100*(1-Ssai[i])

ASji[i] = 100*Sji[i]			# Juveniles  		
BSji[i] = 100*(1-Sji[i])

AShi[i] = 100*Shi[i]		 	# Hatchlings	
BShi[i] = 100*(1-Shi[i])

for(j in 1:t){			# Projection loop; drawing annual demographic rates
	
Sat[i,j]=rbeta(1,ASai[i],BSai[i])			# Adults
Sa1t[i,j]=rbeta(1,ASa1i[i],BSa1i[i])		# Primiparous adults
Ssat[i,j]=rbeta(1,ASsai[i],BSsai[i])		# Subadults
Sjt[i,j]=rbeta(1,ASji[i],BSji[i])			# Juveniles
Sht[i,j]=rbeta(1,aSh,bSh)					# Hatchlings
Set[i,j]=rbeta(1,aSe,bSe)					# Egg survival
Srt[i,j]=rbeta(1,100*Sr, 100*(1-Sr)) 		# Randomize the acclimation effect (Sr; survival reduction)

# Drawing the number of released subadults and juveniles, and account for sex ratio to remove males
if (j>0 && j<nyr+1 && scenarios[h,3]>0) Nsar[i,j]=runif(1,0.8*nrelsa,1.2*nrelsa)*rbeta(1,aSR,bSR)
if (j>0 && j<nyr+1 && scenarios[h,3]>0) Njr[i,j]=runif(1,0.8*nrelj,1.2*nrelj)*rbeta(1,aSR,bSR)

# Specify the number of individuals released in CNF and ABRP to date and projected
if (rownames(scenarios)[h] == "CNFc") Nsar=CNFc*rbeta(1,aSR,bSR)
if (rownames(scenarios)[h] == "ABRPc") Nsar=ABRPc*rbeta(1,aSR,bSR)
if (rownames(scenarios)[h] == "CNFp") Nsar=CNFp*rbeta(1,aSR,bSR)
if (rownames(scenarios)[h] == "ABRPp") Nsar=ABRPp*rbeta(1,aSR,bSR)

# Projection equation for adults
if (j>1) Na[i,j]=round((Na[i,j-1]*Sat[i,j-1]+Na1[i,j-1]*Sa1t[i,j-1]*Ta1a[i,j-1])) else Na[i,j]=0
round(Na[i,j],0)

# Projection equation for primiparous adults
if (j>1) Na1[i,j]=round((Na1[i,j-1]*(Sa1t[i,j-1]*(1-Ta1a[i,j-1])))+(Nsa[i,j-1]*Ssat[i,j-1]*Tsaa1[i,j-1])+(Nsar[i,j-1]*Ssat[i,j-1]*Tsaa1[i,j-1]*Srt[i,j-1])) else Na1[i,j]=0
round(Na1[i,j],0)

# Density-dependence; decrease subadult survival considerably when population sizes get too high;
# drop to 0.40 when abundance > nmax population size OR to 0.20 when abundance > 1.5*nmax
if (Na[i,j]>Nmax) Ssat[i,j]=runif(1,.3,.5)
if (Na[i,j]>Nmax*1.5) Ssat[i,j]=runif(1,.1,.3)

# Projection equation for subadults
if (j>1) Nsa[i,j]=round((Nsa[i,j-1]*(Ssat[i,j-1]*(1-Tsaa1[i,j-1])))+(Nj[i,j-1]*Sjt[i,j-1]*Tjsa[i,j-1])+(Njr[i,j-1]*Sjt[i,j-1]*Tjsa[i,j-1]*Srt[i,j-1])) else Nsa[i,j]=0
round(Nsa[i,j],0)

# Projection equation for juveniles
if (j>1) Nj[i,j]=round((Nj[i,j-1]*(Sjt[i,j-1]*(1-Tjsa[i,j-1])))+(Nh[i,j-1]*Sht[i,j-1]*Thj[i,j-1])) else Nj[i,j]=0
round(Nj[i,j],0)

# Good-year/bad-year function for probability of breeding
if (GY[i,j]==1) Pbt[i,j]=1 else Pbt[i,j]=runif(1,.7,.8)

# Calculate the number of eggs produced each year, while accounting for sex ratio (i.e., dividing by two)
eggsa[i,j] = sum(rpois(round(Na[i,j]*Pbt[i,j],0),muFa))/2
eggsa1[i,j] = sum(rpois(round(Na1[i,j]*Pbt[i,j],0),muFa1))/2
round(eggsa[i,j],0)
round(eggsa1[i,j],0)

# Calculate the number of hatchlings produced each year by accounting for viability (Pva & Pva1) and nest survival (i.e., egg survival; Se) 
if (j>1) Nh[i,j] = round((eggsa[i,j-1]*Pvat[i,j-1]*Set[i,j-1])+(eggsa1[i,j-1]*Pva1t[i,j-1]*Set[i,j-1])) else Nh[i,j]=0
round(Nh[i,j],0)

# Sum the number of reproductive and immature individuals to facilitate the construction of figures
Nbreeders = Na+Na1		# Sum all breeders at each time step
Nimmatures = Nh+Nj+Nsa	# Sum all immatures at each time step

# Implement sampling (detection) probability for adults and immatures
Nobreeders[i,j]=sum(rbinom(Nbreeders[i,j],1,Psamp[i,j]))
Noimmatures[i,j]=sum(rbinom(Nimmatures[i,j],1,Psamp[i,j]))

# Calculate population growth rate
if(j>5) rlam[i,j]=Nbreeders[i,j]/Nbreeders[i,j-1]
if(j>5 && Nobreeders[i,j-1]>0) olam[i,j]= Nobreeders[i,j]/Nobreeders[i,j-1]
if (Nbreeders[i,j]<5) Pext[i,j]=1 else Pext[i,j]=0 

}}	# Close replication and projection loops

# Create an categorical identifying variable of length r for each life stage
stages = c("H","J","SA","PA","A","Imm","Breeders","ObsImm","ObsBreeders")
for (l in 1:length(stages)){
	Stage = as.data.frame(rep(stages[l],r))
	assign(paste0("Stage",stages[l]),Stage)}

# Bind the results from a scenario projection into a data frame
scenario = rbind.data.frame(cbind.data.frame(StageH,Nh),cbind.data.frame(StageJ,Nj),
	cbind.data.frame(StageSA,Nsa),cbind.data.frame(StagePA,Na1),cbind.data.frame(StageA,Na),
	cbind.data.frame(StageImm,Nimmatures),cbind.data.frame(StageBreeders,Nbreeders), 
	cbind.data.frame(StageObsImm,Noimmatures),cbind.data.frame(StageObsBreeders,Nobreeders))
scenario = cbind.data.frame(rep(rownames(scenarios)[h],r*length(stages)),scenario)
colnames(scenario) = c("Scenario","Stage",1:30)


results = rbind(results,scenario)	# Bind results from the scenario to an object with all scenario results

medrlam = apply(rlam, 2, median, na.rm=TRUE)
medolam = apply(olam, 2, median, na.rm=TRUE)
PE = apply(Pext,2,sum)/r
PEt = PE[t]
muNsar = mean(apply(Nsar,1,sum))
muNjr = mean(apply(Njr,1,sum))

assign(paste0("medrlam", rownames(scenarios)[h]), medrlam) 
assign(paste0("medolam", rownames(scenarios)[h]), medolam) 
assign(paste0("PE", rownames(scenarios)[h]), PE) 
assign(paste0("PEt", rownames(scenarios)[h]), PEt) 
assign(paste0("muNsar", rownames(scenarios)[h]), muNsar) 
assign(paste0("muNjr", rownames(scenarios)[h]), muNjr)

} 	# Close scenario loop

results = results[-c(1),]	# Remove top NA row from results




############## Part 3)
############## Project the population growth and extinction risk of D. couperi  
############## under many combinations of repatriation scenarios to generate
############## a 3-D response surface of how extinction prob is influenced by
############## propagule size & release program duration

t = 30		# Number of years to project simulations
r = 100		# Number of simulation replications; 
			# r = 100 here for shorter computing times; ~12 mins

start_time <- Sys.time()

### Population growth & quasi-extinction estimates
rlam=matrix(0,r,t)		# Real (true) population growth
olam=matrix(,r,t)		# Observed population growth, after sampling (detection) process
Pext=matrix(0,r,t)		# Probability of quasi-extinction

# We want to create a 3-D response surface for how extinction probability is influenced by release age (juv vs subadult), propagule size (no. of snakes), and release program duration (no. of years). We will do simulations for both juveniles and subadults and then produce graphs for both stages. We will run simulations for propagule size varying from 1-20 females per year, release duration varying from 1-20 years.

## (1) snake release age (juvenile or subadult)
## (2) annual propagule size (4-300 individuals/year; will be split into females only later)
## (3) release duration (1-20 yr)

#### Generate a matrix with all combinations of the aforementioned variables. This 

require(utils)

releases = rep(seq(4,300,4),20)
zeros = rep(0,1500)
nyr = c(rep(1,75),rep(2,75),rep(3,75),rep(4,75),rep(5,75),	
			rep(6,75),rep(7,75),rep(8,75),rep(9,75),rep(10,75),
			rep(11,75),rep(12,75),rep(13,75),rep(14,75),rep(15,75),
			rep(16,75),rep(17,75),rep(18,75),rep(19,75),rep(20,75))
total = releases*nyr

# create a dataframe for subadult release scenarios
nrelsa = as.data.frame(cbind(releases,zeros,nyr,total))
colnames(nrelsa) = c("nrelsa","nrelj","nyr","total")

# now create a dataframe for juvenile release scenarios; slightly different order
nrelj = as.data.frame(cbind(zeros,releases,nyr,total))
colnames(nrelj) = c("nrelsa","nrelj","nyr","total")

# bind the subadult and juvenile scenarios
releases = rbind(nrelsa,nrelj)

# the multiplication of "total = releases*nyr" creates many scenarios that are impractically large; e.g., releasing 300 individuals/year for 20 years. This will cost a lot of computing time and is biologically unfeasible/impractical. So, we can truncate this to only evaluate scenarios that  release a total of 600 individuals (or, 300 females).
scenarios = releases[ which(releases[,4] < 600),]

n = length(scenarios[,1])
					
### Create an empty matrix to dump all the results into
RESULTS = matrix(0, nrow = 1, ncol = 32, dimnames=list(NA,c("Scenario","Stage",1:30)))
	# Population projection results
ProbExt = matrix(0, nrow = 1, ncol = 4, dimnames=list(NA, c("SubadRel","JuvRel","NoYears","ProbExt")))
	# Probability of extinction results -- for 3D graph -- most important here
	
	
### Use a for-loop to iteratively calculate demography under different management scenarios
### by using each row of the scenarios object to provide unique combinations of parameters
 
for (h in 1:n){			# For-loop for each repatriation scenario
			
nrelsa=scenarios[h,1] 	# Mean no. of 2-yr old subadults releases (both sexes)
nrelj=scenarios[h,2]	# Mean number of 1 yr-old juveniles releases (both sexes)
nyr=scenarios[h,3] 		# Duration of the release program (years)

mPsamp=0.2 	# Probabilty of sampling snakes (i.e., detection probability); 
			# this is hard-coded here for convenience b/c we aren't evaluating det. prob.
SDmPsamp = (mPsamp*0.1)*(mPsamp*0.1)
aPs= mPsamp*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
bPs=(1-mPsamp)*((mPsamp*(1-mPsamp)/SDmPsamp)-1)
Psamp = matrix(rbeta(r*t,aPs,bPs),r,t)

Sr = 0.5			# Acclimation effect on survival (survival reduction) of released individuals	
Srt=matrix(0,r,t)	# also hard coded here for convenience

### Create some matrices and parameters for the model loop:
Nh =  matrix(0,r,t)		# Abundance of hatchlings 
Nj = matrix(0,r,t) 		# Juveniles
Nsa = matrix(0,r,t) 	# Subadults
Na1 = matrix(0,r,t) 	# Primiparous adults
Na = matrix(0,r,t) 		# Adults

Nh[,1]=0		# Initial abundance of hatchlings
Nj[,1]=0		# Juveniles
Nsa[,1]=0		# Subadults
Na1[,1]=0		# Primiparous adults
Na[,1]=0		# Adults

Njr = matrix(0,r,t) 	# Number of captive-reared juveniles (1-yr olds) released
Nsar = matrix(0,r,t) 	# Number of captive-reared subadults (2-yr olds) released

Nobreeders = matrix(0,r,t) 	# Number of adults observed/counted, given Psamp
Noimmatures = matrix(0,r,t)	# Number of immatures observed/counted, given Psamp

Nmax=100 			# Population ceiling for density dependence


#### Population-projection model that accounts for imperfect detection

for(i in 1:r){			# Replication loop; draws replicate-level means for survival at each stage

ASai[i] = 100*Sai[i]			# Adults 		
BSai[i] = 100*(1-Sai[i])		

ASa1i[i] = 100*Sa1i[i]		 	# Primiparous adults	
BSa1i[i] = 100*(1-Sa1i[i])

ASsai[i] = 100*Ssai[i]		 	# Subadults	
BSsai[i] = 100*(1-Ssai[i])

ASji[i] = 100*Sji[i]			# Juveniles  		
BSji[i] = 100*(1-Sji[i])

AShi[i] = 100*Shi[i]		 	# Hatchlings	
BShi[i] = 100*(1-Shi[i])

for(j in 1:t){			# Projection loop; drawing annual demographic rates
	
Sat[i,j]=rbeta(1,ASai[i],BSai[i])			# Adults
Sa1t[i,j]=rbeta(1,ASa1i[i],BSa1i[i])		# Primiparous adults
Ssat[i,j]=rbeta(1,ASsai[i],BSsai[i])		# Subadults
Sjt[i,j]=rbeta(1,ASji[i],BSji[i])			# Juveniles
Sht[i,j]=rbeta(1,aSh,bSh)					# Hatchlings
Set[i,j]=rbeta(1,aSe,bSe)					# Egg survival
Srt[i,j]=rbeta(1,100*Sr, 100*(1-Sr)) 		# Randomize the acclimation effect (Sr; survival reduction)

# Drawing the number of released subadults and juveniles, and account for sex ratio to remove males
if (j>0 && j<nyr+1) Nsar[i,j]=runif(1,0.8*nrelsa,1.2*nrelsa)*rbeta(1,aSR,bSR)
if (j>0 && j<nyr+1) Njr[i,j]=runif(1,0.8*nrelj,1.2*nrelj)*rbeta(1,aSR,bSR)

# Projection equation for adults
if (j>1) Na[i,j]=round((Na[i,j-1]*Sat[i,j-1]+Na1[i,j-1]*Sa1t[i,j-1]*Ta1a[i,j-1])) else Na[i,j]=0
round(Na[i,j],0)

# Projection equation for primiparous adults
if (j>1) Na1[i,j]=round((Na1[i,j-1]*(Sa1t[i,j-1]*(1-Ta1a[i,j-1])))+(Nsa[i,j-1]*Ssat[i,j-1]*Tsaa1[i,j-1])+(Nsar[i,j-1]*Ssat[i,j-1]*Tsaa1[i,j-1]*Srt[i,j-1])) else Na1[i,j]=0
round(Na1[i,j],0)

# Density-dependence; decrease subadult survival considerably when population sizes get too high;
# drop to 0.40 when abundance > nmax population size OR to 0.20 when abundance > 1.5*nmax
if (Na[i,j]>Nmax) Ssat[i,j]=runif(1,.3,.5)
if (Na[i,j]>Nmax*1.5) Ssat[i,j]=runif(1,.1,.3)

# Projection equation for subadults
if (j>1) Nsa[i,j]=round((Nsa[i,j-1]*(Ssat[i,j-1]*(1-Tsaa1[i,j-1])))+(Nj[i,j-1]*Sjt[i,j-1]*Tjsa[i,j-1])+(Njr[i,j-1]*Sjt[i,j-1]*Tjsa[i,j-1]*Srt[i,j-1])) else Nsa[i,j]=0
round(Nsa[i,j],0)

# Projection equation for juveniles
if (j>1) Nj[i,j]=round((Nj[i,j-1]*(Sjt[i,j-1]*(1-Tjsa[i,j-1])))+(Nh[i,j-1]*Sht[i,j-1]*Thj[i,j-1])) else Nj[i,j]=0
round(Nj[i,j],0)

# Good-year/bad-year function for probability of breeding
if (GY[i,j]==1) Pbt[i,j]=1 else Pbt[i,j]=runif(1,.7,.8)

# Calculate the number of eggs produced each year, while accounting for sex ratio (i.e., dividing by two)
eggsa[i,j] = sum(rpois(round(Na[i,j]*Pbt[i,j],0),muFa))/2
eggsa1[i,j] = sum(rpois(round(Na1[i,j]*Pbt[i,j],0),muFa1))/2
round(eggsa[i,j],0)
round(eggsa1[i,j],0)

# Calculate the number of hatchlings produced each year by accounting for viability (Pva & Pva1) and nest survival (i.e., egg survival; Se) 
if (j>1) Nh[i,j] = round((eggsa[i,j-1]*Pvat[i,j-1]*Set[i,j-1])+(eggsa1[i,j-1]*Pva1t[i,j-1]*Set[i,j-1])) else Nh[i,j]=0
round(Nh[i,j],0)

# Sum the number of reproductive and immature individuals to facilitate the construction of figures
Nbreeders = Na+Na1		# Sum all breeders at each time step
Nimmatures = Nh+Nj+Nsa	# Sum all immatures at each time step

# Implement sampling (detection) probability for adults and immatures
Nobreeders[i,j]=sum(rbinom(Nbreeders[i,j],1,Psamp[i,j]))
Noimmatures[i,j]=sum(rbinom(Nimmatures[i,j],1,Psamp[i,j]))

# Calculate population growth rate
if(j>5) rlam[i,j]=Nbreeders[i,j]/Nbreeders[i,j-1]
if(j>5 && Nobreeders[i,j-1]>0) olam[i,j]= Nobreeders[i,j]/Nobreeders[i,j-1]
if (Nbreeders[i,j]<5) Pext[i,j]=1 else Pext[i,j]=0 

}}	# Close replication and projection loops

# Create an categorical identifying variable of length r for each life stage
stages = c("H","J","SA","PA","A","Imm","Breeders","ObsImm","ObsBreeders")
for (l in 1:length(stages)){
	Stage = as.data.frame(rep(stages[l],r))
	assign(paste0("Stage",stages[l]),Stage)}

# Bind the results from a scenario projection into a data frame
scenario = rbind.data.frame(cbind.data.frame(StageH,Nh),cbind.data.frame(StageJ,Nj),
	cbind.data.frame(StageSA,Nsa),cbind.data.frame(StagePA,Na1),cbind.data.frame(StageA,Na),
	cbind.data.frame(StageImm,Nimmatures),cbind.data.frame(StageBreeders,Nbreeders), 
	cbind.data.frame(StageObsImm,Noimmatures),cbind.data.frame(StageObsBreeders,Nobreeders))
scenario = cbind.data.frame(rep(rownames(scenarios)[h],r*length(stages)),scenario)
colnames(scenario) = c("Scenario","Stage",1:30)

RESULTS = rbind(RESULTS,scenario)	
# Bind results from the scenario to an object with all scenario results


medrlam = apply(rlam, 2, median, na.rm=TRUE)
medolam = apply(olam, 2, median, na.rm=TRUE)
PE = apply(Pext,2,sum)/r
PEt = PE[t]
muNsar = mean(apply(Nsar,1,sum))
muNjr = mean(apply(Njr,1,sum))

assign(paste0("medrlam", rownames(scenarios)[h]), medrlam) 
assign(paste0("medolam", rownames(scenarios)[h]), medolam) 
assign(paste0("PE", rownames(scenarios)[h]), PE) 
assign(paste0("PEt", rownames(scenarios)[h]), PEt) 
assign(paste0("muNsar", rownames(scenarios)[h]), muNsar) 
assign(paste0("muNjr", rownames(scenarios)[h]), muNjr)

# Create a probability of extinction (pe) matrix row that summarizes the total number of individuals released, the release program duration, and the resulting probability of extinction
pe = as.matrix(0, nrow=1, ncol=4, dimnames=list(NA, c("SubadRel","JuvRel","NoYears","ProbExt")))
pe[1:4] = c(muNsar, muNjr, nyr, PEt)

# Bind the probability of extinction (pe) matrix row to an object that has all of the probability of extinction results together
ProbExt = rbind(ProbExt, pe)

} 	# Close scenario loop

RESULTS = RESULTS[-c(1),]	# Remove top NA row from RESULTS
ProbExt = ProbExt[-c(1),]	# Remove top NA row from ProbExt

end_time <- Sys.time()
end_time - start_time		# See how long model takes to run




############## Part 4)
############## Visualizing the results with Tables and Figures

#########
######### TABLES
#########

### Table 1
# Number of snakes released per year in Conecuh National Forest and Apalachicola Bluffs & Ravines.
# Manually formatted in an Excel file, but contains the following information:
# Site	2010		2011		2012		2013		2014		2015		2016		2017		2018		Total
# CNF	17		30		31		20		0		9		0		26		24		157
# ABRP	-		-		-		-		-		-		-		12		20		32

### Table 2
# Tabulate total number released, true population growth (lambda), observed population growth (lambda obs), quasi-extinction probability for six reintroduction strategies 

cols = c("Scenario","ReleaseDuration","TotalReleased","Lambda","LambdaObs","Pext")
rows = c("A","C","E","G","B","D","F","H")
tableX1 = matrix(c("A", 10, round(muNsarA), mean(medrlamA[-c(1:5)]), mean(medolamA[-c(1:5)]),PEtA,
	"C", 5, round(muNsarC), mean(medrlamC[-c(1:5)]), mean(medolamC[-c(1:5)]),PEtC,
	"E", 10, round(muNsarE), mean(medrlamE[-c(1:5)]), mean(medolamE[-c(1:5)]),PEtE,
	"G", 5, round(muNsarG), mean(medrlamG[-c(1:5)]), mean(medolamG[-c(1:5)]),PEtG,
	"B", 10, round(muNjrB), mean(medrlamB[-c(1:5)]), mean(medolamB[-c(1:5)]),PEtB,
	"D", 5, round(muNjrD), mean(medrlamD[-c(1:5)]), mean(medolamD[-c(1:5)]),PEtD,
	"F", 10, round(muNjrF), mean(medrlamF[-c(1:5)]), mean(medolamF[-c(1:5)]),PEtF,
	"H", 5, round(muNjrH), mean(medrlamH[-c(1:5)]), mean(medolamH[-c(1:5)]),PEtH),
	nrow=8,ncol=6, byrow=TRUE, dimnames=list(rows,cols))

# Copy and paste top-model set into an Excel file to tidy for publication
clip = pipe("pbcopy","w")
write.table(tableX1, file=clip, sep="\t", row.names=FALSE)
close(clip)


### Table 3
# Summarize simulated population projections for CNF and ABRP given current and anticipated releases

cols = c("Site","Scenario","NoReleased","Lambda","LambdaObs","Pext","MedAbun30")
rows = c("CNFc","CNFp","ABRPc","ABRPp")
tableX2 = matrix(c(
"AL", "Current", round(2*sum(CNFc[1,])), mean(medrlamCNFc[6:30]), mean(medolamCNFc[6:30]), PECNFc[30], median(subset(results, results$Scenario == "CNFc" & Stage == "Breeders")[,32]),
"AL", "Current + futures", round(2*sum(CNFp[1,])), mean(medrlamCNFp[6:30]), mean(medolamCNFp[6:30]), PECNFp[30], median(subset(results, results$Scenario == "CNFp" & Stage == "Breeders")[,32]),
"FL", "Current", round(2*sum(ABRPc[1,])), mean(medrlamABRPc[6:30]), mean(medolamABRPc[6:30]), PEABRPc[30], median(subset(results, results$Scenario == "ABRPc" & Stage == "Breeders")[,32]),
"FL", "Current + futures", round(2*sum(ABRPp[1,])), mean(medrlamABRPp[6:30]), mean(medolamABRPp[6:30]), PEABRPp[30], median(subset(results, results$Scenario == "ABRPp" & Stage == "Breeders")[,32])),
	nrow=4, ncol=7, byrow=TRUE, dimnames=list(rows,cols))
tableX2

# Copy and paste top-model set into Excel file to tidy for publication
clip = pipe("pbcopy","w")
write.table(tableX2, file=clip, sep="\t", row.names=FALSE)
close(clip)


### Table 4
# Tabulate changes in metrics related to differences in acclimation effect
# Scenarios = A, 0.5; M, 0.8; N, 0.6; O, 0.4; P, 0.2

cols = c("Scenario","Lambda","LambdaChange","Pext","PextChange")
rows = c("A","M","N","O","P")
tableX2 = matrix(c("A", mean(medrlamA[-c(1:5)]), "NA",PEtA,"NA",
	"M", mean(medrlamM[-c(1:5)]), ((mean(medrlamM[-c(1:5)])/(mean(medrlamA[-c(1:5)])))-1)*100,PEtM,((PEtM/PEtA)-1)*100,
	"N", mean(medrlamN[-c(1:5)]), ((mean(medrlamN[-c(1:5)])/(mean(medrlamA[-c(1:5)])))-1)*100,PEtN,((PEtN/PEtA)-1)*100,
	"O", mean(medrlamO[-c(1:5)]), ((mean(medrlamO[-c(1:5)])/(mean(medrlamA[-c(1:5)])))-1)*100,PEtO,((PEtO/PEtA)-1)*100,
	"P", mean(medrlamP[-c(1:5)]), ((mean(medrlamP[-c(1:5)])/(mean(medrlamA[-c(1:5)])))-1)*100,PEtP,((PEtP/PEtA)-1)*100),
	nrow=5,ncol=5, byrow=TRUE, dimnames=list(rows,cols))

# Copy and paste top-model set into Excel file to tidy for publication
clip = pipe("pbcopy","w")
write.table(tableX2, file=clip, sep="\t", row.names=FALSE)
close(clip)



#########
######### FIGURES
#########


### Figure 1
# A picture of an Eastern Indigo Snake, a map of the current and historic distribution of the species, and a diagram of the population model built here. This figure was largely built in Microsoft Powerpoint, but the map was made using a separate script in R. 


### Figure 2
# Plot simulated adult abundance under different reintroduction scenarios using a multi-panel plot with different reintroduction scenarios
par(mfrow=c(4,2), oma=c(5,4,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7) 
scenarios = c("A","B","C","D","E","F","G","H")
max = 120	# set max for y axis & panel labels
			# will have to adjust manually under different simulation scenarios
for (i in 1:length(scenarios)){
	scenario = scenarios[i]
	plot(1:t, apply(subset(results, Scenario == scenario & Stage == "Breeders")[,-c(1:2)],2,median),
		type="l", lty=1, lwd=4, ylim=c(0,max), axes=FALSE, cex.lab=1.4, cex.axis=2, 
		xlab="Time (years)", ylab="Abundance")
	if (i > 6){	# restricts axis labels to bottom row
		axis(side=1, c(0,5,10,15,20,25,30), lwd=4, cex.axis=1.3)} 
		else {axis(side=1, labels=FALSE, lwd=4, cex.axis=1.3)}
	if (i %% 2 != 0){
		axis(side=2, lwd=4, cex.axis=1.3)} else {axis(side=2, labels=FALSE,lwd=4, cex.axis=1.3)}
	lines(1:t, apply(subset(results, Scenario == scenario & Stage == "Breeders")[,-c(1:2)],
		2,quantile,probs=c(0.025)), lty=3, lwd=4, col="darkgrey")
	lines(1:t, apply(subset(results, Scenario == scenario & Stage == "Breeders")[,-c(1:2)],
		2,quantile,probs=c(0.975)), lty=3, lwd=4, col="darkgrey")
	if(scenario == "B") {legend(4,0.97*max, c("Median abundance","95% CI"), inset=0.05, cex=1,
		lty=c(1,3), col=c("black","darkgrey"),
		lwd=c(3,3), box.lwd=2)}
	text(2,0.93*max, scenarios[i], cex=2)
	title(xlab = "Time (years)", ylab = "Abundance", outer=TRUE, line=3, cex.sub=2)
}

### Figure 3
# Single panel graph to evaluate how quasi-extinction probability varies through time under different management scenarios
plot(1:t, PEA, type="l", lty=1, lwd=4, col="black", ylim=c(0,1), axes=FALSE, cex.lab=1.2, cex.axi=2, xlab="Time (years)", ylab="Probability of quasi-extinction")	# 30 subadults; 10 yr
axis(1, c(0,5,10,15,20,25,30), lwd=4, cex.axis=1.3)
axis(2, lwd=4, cex.axis=1.3)
lines(1:t, PEC, lty=3, lwd=4, col="firebrick1")	# 30 subadults; 5 yr
lines(1:t, PEE, lty=5, lwd=4, col="cornflowerblue") # "15 subadults; 10 yr",
lines(1:t, PEG, lty=7, lwd=4, col="gold") # 15 subadults; 5 yr
lines(1:t, PEB, lty=2, lwd=4, col="gray48") # 30 juveniles; 10 yr
lines(1:t, PED, lty=4, lwd=4, col="chocolate1") # 30 juveniles; 5 yr
lines(1:t, PEF, lty=6, lwd=4, col="darkgreen") # 15 juveniles; 10 yr
lines(1:t, PEH, lty=8, lwd=4, col="blue1") # 15 juveniles; 5 yr
legend("topright", inset=0.01, c("30 SA; 10 yr","30 SA; 5 yr","15 SA; 10 yr","15 SA; 5 yr","30 J; 10 yr","30 J; 5 yr","15 J; 10 yr","15 J; 5 yr"),cex=1,lty=c(1,3,5,7,2,4,6,8),col=c("black","firebrick1","cornflowerblue","gold","gray48","chocolate1","darkgreen","blue1"), lwd=c(3,3), box.lwd=2)


#### Figure 4
# Use a 3-D surface plot to illustrate extinction probability under a variety of release scenarios with data from the object 'ProbExt'

### SCATTER3D
install.packages("plot3D", dependencies=TRUE)
library(plot3D)

### plot Probability of EXTINCTION all subadult and juvenile scenarios in two-panel graph
par(mfrow=c(1,2), oma=c(4,0,1,2)+0.5, cex.lab=1) 

res = subset(ProbExt, ProbExt[,1] != 0)	#subset to subadult releases
scatter3D(res[,1], res[,3], res[,4], clab=expression(italic(P[e])), pch=19, bty="g", 	
	theta=125, phi=5, xlab="No. of females", ylab="No. of years", zlab="Prob. of extinction",	
	ticktype="detailed",  lwd=1, cex.axis=0.6, cex.lab=0.8, type="h",
	colkey=list(length=0.6,width=1,cex.clab=1,cex.axis=0.7,dist=0.1))
text(locator(1),"A",cex=2)

res = subset(ProbExt, ProbExt[,2] != 0)	#subset to juvenile releases
scatter3D(res[,2], res[,3], res[,4], clab=expression(italic(P[e])), pch=19, bty="g", theta=125, phi=5, xlab="No. of females", ylab="No. of years", zlab="Prob. of extinction", ticktype="detailed", cex.axis=0.6, cex.lab=0.8, colkey=list(plot=FALSE), type="h")
text(locator(1),"B",cex=2)

### plot all scenarios that resulted in extinction < 0.1 threshold
par(mfrow=c(1,2), oma=c(4,0,1,2)+0.5, cex.lab=1) 

res = subset(ProbExt, ProbExt[,1] != 0)	#subset to subadult releases
res = subset(res, res[,4] < 0.2)
scatter3D(res[,1], res[,3], res[,4], clab=expression(italic(P[e])), pch=19, bty="g", 	
	theta=125, phi=5, xlab="No. of females", ylab="No. of years", zlab="Prob. of extinction",	
	ticktype="detailed",  lwd=1, cex.axis=0.6, cex.lab=0.8,
	colkey=list(length=0.6,width=1,cex.clab=1,cex.axis=0.7,dist=0.1), type="h")
text(locator(1),"C",cex=2)

res = subset(ProbExt, ProbExt[,2] != 0)	#subset to juvenile releases
res = subset(res, res[,4] < 0.2)
scatter3D(res[,2], res[,3], res[,4], clab=expression(italic(P[e])), pch=19, bty="g", theta=125, phi=5, xlab="No. of females", ylab="No. of years", zlab="Prob. of extinction", ticktype="detailed", cex.axis=0.6, cex.lab=0.8, colkey=list(plot=FALSE), type="h")
text(locator(1),"D",cex=2)

detach(package:plot3D)



### Figure 5
# Use a multipanel graph to do two things: (1) estimate ongoing repatration projections, & 
# (2) evaluate how detection limits our ability to monitor population growth in the field. 
# Estimate current population projection trends for CNF and ABRP and plot
# true abundance AND the number expected to be observed after detection process
# where detection = 0.05, 0.15, 0.25, 0.35. Assume subadult only releases & acclimation = 0.5. 
# Project the current scenarios for CNF and ABRP (A & C), 
# & also scenarios including tentatively planned releases at the respective sites (B & D).
par(mfrow=c(2,2), oma=c(5,4,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7)
scenarios = c("CNFc","CNFp","ABRPc","ABRPp")
labels = c("A","B","C","D")
max = 20	# set max for y axis & panel labels; may have to adjust manually to fit results
for (i in 1:length(scenarios)){
	scenario = scenarios[i]
	plot(1:t, apply(subset(results, Scenario == scenario & Stage == "Breeders")[,-c(1:2)],2,median),
		type="l", lty=1, lwd=4, ylim=c(0,max), axes=FALSE, cex.lab=1.4, cex.axis=2, 
		xlab="Time (years)", ylab="Abundance")
	lines(1:t, 0.5*apply(subset(results, Scenario == scenario & Stage == "ObsBreeders")[,-c(1:2)],
		2,median), lty=2, lwd=4, col="gray48")
	lines(1:t, apply(subset(results, Scenario == scenario & Stage == "ObsBreeders")[,-c(1:2)],
		2,median), lty=3, lwd=4, col="gray48")
	lines(1:t, 2*apply(subset(results, Scenario == scenario & Stage == "ObsBreeders")[,-c(1:2)],
		2,median), lty=4, lwd=4, col="gray48")
	lines(1:t, 3*apply(subset(results, Scenario == scenario & Stage == "ObsBreeders")[,-c(1:2)],
		2,median), lty=5, lwd=4, col="gray48")
	if (i > 2){			# restricts x-axis labels to bottom row
		axis(side=1, c(0,5,10,15,20,25,30), lwd=4, cex.axis=1.3)} 
		else {axis(side=1, labels=FALSE, lwd=4, cex.axis=1.3)}
	if (i %% 2 != 0){	# restricts y-axis labels to left column
		axis(side=2, lwd=4, cex.axis=1.3)} 
		else {axis(side=2, labels=FALSE,lwd=4, cex.axis=1.3)}
	if(i == 1){			# adds legend to first plot
		legend(2,0.85*max, 
		c(expression("True N"),expression(N[o]~' - '~italic( p)~'= 0.05'), 
		expression(N[o]~' - '~italic( p)~'= 0.15'),expression(N[o]~' - '~italic( p)~'= 0.25'),
		expression(N[o]~' - '~italic( p)~'= 0.35')), inset=0.01,
		cex=1,lty=c(1,2,3,4,5), col=c("black","gray48","gray48","gray48","gray48"),
		lwd=c(3,3,3,3,3), box.lwd=2)
		}
	# manually label extinction probabilities for each scenario; obtain probs from e.g., PECNFc[30]
	if(i == 1){text(24, 0.94*max, expression(italic(P[e])~' = 0.33'), cex=1.2)}
	if(i == 2){text(24, 0.94*max, expression(italic(P[e])~' = 0.18'), cex=1.2)}
	if(i == 3){text(24, 0.94*max, expression(italic(P[e])~' = 0.67'), cex=1.2)}
	if(i == 4){text(24, 0.94*max, expression(italic(P[e])~' = 0.20'), cex=1.2)}
	text(3,0.95*max, labels[i], cex=2)
	title(xlab = "Time (years)", ylab = "Abundance", outer=TRUE, line=3, cex.sub=2)
}


### How many simulations still observe 0 individuals through each year of the projects?
### for scenario A (subadults, 10 years, 30 ind/yr)

zeros = matrix(,1,30)
for (i in 1:30){
	column = subset(results, Scenario == "O" & Stage == "ObsBreeders")[,i+2]
	column[ column > 0] = 1
	zeros[1,i] = sum(column)/length(column)}
zeros = 1-zeros

zeros[5]*100		# percent of simulations with no observed individuals after 5 years
zeros[10]*100		# 10 years
zeros[15]*100 		# 15 years
zeros[20]*100 		# 20 years
zeros[25]*100 		# 25 years
zeros[30]*100 		# 30 years



























