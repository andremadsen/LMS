

### LMS digital growth curve modelling ###




# Prepare dataset (e.g. Excel) and Save as .csv format:
# R does not allow spaces [ ] or hyphens [-] in column names (i.e. Excel top row); instead use underlines: e.g. 'Albumin_G' (NOT 'Albumin G' or 'Albumin-G')
# R does not allow column names to start with numbers; do naming workarounds: e.g. 's_17OHP' (NOT '17OHP' or '17-OHP')

# Dataset must contain a column to specify child gender: Example below: <Gender> 1=male, 2=female
# Dataset must contain a column to specify child age: Example below: <Age_yrs> (numerical scale or decimal number)


# Import .csv dataset into R/RStudio
Data <- read.csv(file.choose(), header=T, sep=",", dec=".") # European .csv has sep=";" and dec=","
Data_Boys  <- Data[Data$Gender == 1,] # New data frame with male observations
Data_Girls <- Data[Data$Gender == 2,] # New data frame with female observations


# Required packages
install.packages("gamlss")
install.packages("ggplot")
library(gamlss) # or manually enable all 'gamlss', 'gamlss.data' and 'gamlss.dist' packages
library(ggplot2)


#################################### 
# LMS MODEL: FEMALE SHBG example
#################################### 

plot.default(Data_Girls$Age_yrs, Data_Girls$s_SHBG) #plot Age vs Hormone levels to visualize the setup and confirm correct variable names

# Establish new data frame ('b') with no missing values for girls' <Age_yrs> and <s_SHBG>
hormone <- Data_Girls$s_SHBG
age     <- Data_Girls$Age_yrs
keep    <- !is.na(hormone)&!is.na(age)
hormone <- hormone[keep]
age     <- age[keep]
b       <- data.frame(age, hormone)


# Apply LMS to the Hormone~Age data frame
Girls_SHBG <- gamlss(hormone~cs(age, df=2), sigma.fo=~cs(age, df=0), nu.fo=~cs(age, df=0), family=BCCG, data=b)

# Visualize the LMS model: the centiles below correspond to -2,-1,mean,+1 and +2 SD curves and can be customized
centiles(Girls_SHBG, age, cent=c(2.28, 15.87, 50, 84.13, 97.72), main = "", ylab ="SHBG, nmol/L", xlab= "Age, years", box(lwd=2),
         lwd.centiles = 2, col.centiles = c("black", "black","red","black", "black"))

# If the centiles are too loose/tight/volatile, re-run the LMS above and customize degrees of freedom (max df=4 for any parameter) 
# Generally, tune the degrees of freedom one at a time, in increments of 1, moving left to right in the gamlss(...) step above
# Continue when you have a visually good result (trial-and-error)
# If the hormone is changing exponentially with age, try gamlss(... , family=BCCGo , ...) instead


#Quality control
summary(Girls_SHBG)

rqres.plot(Girls_SHBG, howmany = 6, plot.type = "all", type = "QQ") # Detrend Q-Q plot

wp(Girls_SHBG) # Worm plot: a series of detrended Q-Q plots, split by covariate levels



#Z-score residuals and quality control
b$SHBG_sds  <- residuals(Girls_SHBG, what = "z-scores", type = "simple", terms=NULL)                      # retrieve z-scores for all observations in the LMS model
plot.default(b$age, b$SHBG_sds, xlab="Age, years", ylab="SHBG  z-scores")                                 # verify an even distribution of z-scores by age
hist(b$SHBG_sds, breaks=30, xlab="SHBG  z-scores", ylab="Frequency", main="Histogram of SHBG z-scores")   # verify an even distribution of z-scores in total 
b$LMSifelse <- as.factor(ifelse(b$SHBG_sds >= 1, 2, ifelse(b$SHBG_sds < -1, 1,NA)))                       # group z-scores by value for later visualization



#Z-score control verification
L <- predict(Girls_SHBG, what="nu", type="response", newdata=b)             # retrieve L parameter values in the model
M <- predict(Girls_SHBG, what="mu", type="response", newdata=b)             # retrieve M parameter values in the model
S <- predict(Girls_SHBG, what="sigma", type="response", newdata=b)          # retrieve S parameter values in the model
z <- residuals(Girls_SHBG, what = "z-scores", type = "simple", terms=NULL)  # retrieve z-scores for observations in the model [directly]
d <- data.frame(cbind(b, L,M,S,z))                                          # combine above columns with the previous dataframe ('b') 
d$calcZ <- (((d$hormone/d$M)^d$L)-1)/(d$L*d$S)                              # manually calculate z-scores from the L,M,S values [indirectly]
plot.default(d$z, d$calcZ)                                                  # should be a perfectly straight line: if not, reduce LMS model degrees of freedom





#GGPLOT INTERGRATION
newx <- seq(6,16,0.1) #Girls in my dataset are age 6.0-16.0 years, and we need a sequence in increments of 0.1 between those limits
mat2 <- centiles.pred(Girls_SHBG, xname="age", xvalues=newx, type="standard-centiles") #Retrieve centiles' XY coordinates from LMS model
tail(mat2)            #verify that you got the data and note the LAST row number (Example here has 101 rows of data)

#Create the ggplot figure: requires customization of axis limits and names depending on hormone
ggplot(b,aes(age,hormone,col=LMSifelse)) + geom_point(shape=15, fill="gray40", size=0.8) + theme_bw() + labs(title="Female SHBG") + 
  scale_x_discrete(name ="Age, y", limits=c(6:16)) + theme(text = element_text(size=14)) +  
  scale_y_continuous(name ="Serum SHBG, nmol/L", limits=c(0,200)) +
  guides(col = guide_legend(title = "Visual z-scores")) + scale_color_hue(labels = c("< -1 SDS", "> +1 SDS", "normal range")) +
  geom_smooth(aes(x=mat2$age, y=mat2$`-2`), data=mat2, inherit.aes = FALSE, stat="identity",linetype="dashed", color = "gray10") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`-1`),data=mat2, inherit.aes = FALSE,stat="identity",linetype="dashed",color = "gray10") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`0`),data=mat2, inherit.aes = FALSE,stat="identity" ,color = "black") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`1`),data=mat2, inherit.aes = FALSE,stat="identity",linetype="dashed",color = "gray10") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`2`),data=mat2, inherit.aes = FALSE,stat="identity",linetype="dashed",color = "gray10") +
  annotate("text", x=16.8, y=mat2[101,"-2"], label = "-2 SD", color = "gray10") +
  annotate("text", x=16.8, y=mat2[101,"-1"] , label = "-1 SD", color = "gray10") +
  annotate("text", x=16.8, y=mat2[101,"0"] , label = "Mean", color = "black") +
  annotate("text", x=16.8, y=mat2[101,"1"] , label = "+1 SD", color = "gray10") +
  annotate("text", x=16.8, y=mat2[101,"2"] , label = "+2 SD", color = "gray10") 


#Print the ggplot figure (also supports .pdf and .jpeg)
ggsave("Female_SHBG.tiff", dpi = 600)





#################################### 
# LMS MODEL FOR ORDINAL SCALE TANNER B STAGES B1-5
#################################### 

#Girls' Tanner B stage status was specified by dataset variable <Tanner_B> (as.numeric, NOT as.factor)
hormone <- Data_Girls$s_SHBG
TannerB <- Data_Girls$Tanner_B
keep <- !is.na(hormone)&!is.na(TannerB)
hormone <- hormone[keep]
TannerB <- TannerB[keep]
b <- data.frame(TannerB, hormone)
Girls_SHBG <- gamlss(hormone~cs(TannerB, df=2), sigma.fo=~cs(TannerB, df=0), nu.fo=~cs(TannerB, df=0), family=BCCG, data=b)

centiles(Girls_SHBG, TannerB, cent=c(2.28, 15.87, 50, 84.13, 97.72), main = "", ylab ="SHBG, nmol/L", xlab= "Age, years", 
         lwd.centiles = 2, col.centiles = c("black", "black","red","black", "black"))






#################################### 
# APPLY LMS ON AN ISOLATED DISTRIBUTION: EXAMPLE: SHBG LEVELS FOR GIRLS TANNER STAGE = B2
#################################### 

B2     <- Data_Girls[Data_Girls$Tanner_B == 2,]   # new data frame ('B2') to isolate girls in <Tanner_B> stage 2
mean   <- mean(B2$s_SHBG, na.rm=TRUE)             # SHBG Mean   = 84.07 
median <- median(B2$s_SHBG, na.rm=TRUE)           # SHBG Median = 78.00
sd     <- sd(B2$s_SHBG, na.rm=TRUE)               # SHBG SD     = 33.10



L <- (mean - median) / sd          # L parameter = Skew 
M <-  mean                         # M parameter = Mean
S <-  sd / mean                    # S parameter = Coefficient of variation

X <- 84                            # specify a hormone value to calculate its z-score
z_score <- (((X/M)^L)-1)/(L*S)     # formula to calculate z-score
z_score                            # print the z-score value


summary(B2$Age_yrs)                             # note the Min & Max values for Age in the distribution
stage <- data.frame(seq(7.7, 14.1, 0.1))        # draw a sequence between Age Min & Max values in increments of 0.1
colnames(stage) <- "Limits_Tanner_B2"           # rename the sequence 

stage$'-2' <- (M * (1+(L*S*-2))^(1/L))          # formula to calculate Y-coordinate of the -2 SD centile
stage$'-1' <- (M * (1+(L*S*-1))^(1/L))          # formula to calculate Y-coordinate of the -1 SD centile
stage$'0'  <- (M * (1+(L*S*0))^(1/L))           # formula to calculate Y-coordinate of the 0 SD <mean> centile
stage$'1'  <- (M * (1+(L*S*1))^(1/L))           # formula to calculate Y-coordinate of the +1 SD centile
stage$'2'  <- (M * (1+(L*S*2))^(1/L))           # formula to calculate Y-coordinate of the +2 SD centile

View(stage)   #contains the X,Y coordinates for centiles 


#Create the ggplot figure
ggplot(B2,aes(Age_yrs,s_SHBG)) + geom_point(shape=15, fill="#00AFBB", col="#00AFBB", size=1.4) + theme_bw() + 
  labs(title="Tanner B2 stage-specific reference: Female SHBG") + 
  scale_x_discrete(name ="Age, y", limits=c(6:16)) + theme(text = element_text(size=14)) +  
  scale_y_continuous(name ="Serum SHBG, nmol/L", limits=c(0,200)) + 
  geom_smooth(aes(x=stage$Limits_Tanner_B2, y=stage$`-2`), data=stage, inherit.aes = FALSE, stat="identity", color = "gray10") + 
  geom_smooth(aes(x=stage$Limits_Tanner_B2, y=stage$`-1`),data=stage, inherit.aes = FALSE,stat="identity", color = "gray10") + 
  geom_smooth(aes(x=stage$Limits_Tanner_B2, y=stage$`0`),data=stage, inherit.aes = FALSE,stat="identity", color = "red") + 
  geom_smooth(aes(x=stage$Limits_Tanner_B2, y=stage$`1`),data=stage, inherit.aes = FALSE,stat="identity", color = "gray10") + 
  geom_smooth(aes(x=stage$Limits_Tanner_B2, y=stage$`2`),data=stage, inherit.aes = FALSE,stat="identity", color = "gray10") +
  annotate("text", x=16.8, y=stage[1,"-2"], label = "-2 SD", color = "gray10") +
  annotate("text", x=16.8, y=stage[1,"-1"], label = "-1 SD", color = "gray10") +
  annotate("text", x=16.8, y=stage[1,"0"], label = "Mean", color = "red") +
  annotate("text", x=16.8, y=stage[1,"1"], label = "+1 SD", color = "gray10") +
  annotate("text", x=16.8, y=stage[1,"2"], label = "+2 SD", color = "gray10") +
  annotate("text", x=6, y=stage[1,"2"], label = ".", color = "white") +
  geom_point(Data_Girls, mapping=aes(Age_yrs, s_SHBG), mapping.aes = F, inherit.aes = F, 
             inherit.data = F, shape=20, fill="gray50", color="gray50", size=0.6)

#Print figure
ggsave("Female_B2_SHBG_LMS.jpg", dpi = 600)









