

### LMS digital reference curve modelling ###

Interactive reference curves are now available at the Anylite eHealth platform | https://anylite.io and https://anylite.no 

# Prepare dataset (e.g. Excel) and Save as .csv format:
# R does not allow spaces [ ] or hyphens [-] in column names (i.e. Excel top row); instead use underlines: e.g. 'C_reactive_protein'
# R does not allow column names to start with numbers; do naming workarounds: e.g. 's_17OHP' (NOT '17OHP' or '17-OHP')

# Dataset must contain a column to specify child gender: Example below: <Gender> 1=male, 2=female
# Dataset must contain a column to specify child age: Example below: <Age_yrs> (numerical scale or decimal number)
# Any number of columns representing numerical measures (e.g. blood sample results, hormone levels) are modelled sequentially as shown below


# Import .csv dataset into R/RStudio
Data <- read.csv(file.choose(), header=T, sep=",", dec=".")    # European .csv has sep=";" and dec=","
Data_Boys  <- Data[Data$Gender == 1,]                          # New data frame with male observations
Data_Girls <- Data[Data$Gender == 2,]                          # New data frame with female observations


# Required packages
install.packages("gamlss")
install.packages("ggplot2")
library(gamlss) 
library(ggplot2)


#################################### 
# LMS MODEL: FEMALE SHBG example
#################################### 

plot.default(Data_Girls$Age_yrs, Data_Girls$s_SHBG) #plot Age vs Hormone levels to visualize the setup and confirm correct variable names

# Establish new data frame ('b') with no missing values for female <Age_yrs> and <s_SHBG>
hormone <- Data_Girls$s_SHBG
age     <- Data_Girls$Age_yrs
keep    <- !is.na(hormone)&!is.na(age)&hormone>0
hormone <- hormone[keep]
age     <- age[keep]
b       <- data.frame(age, hormone)


# Apply LMS to the Hormone~Age data frame
LMS_obj <- gamlss(hormone~cs(age, df=2), sigma.fo=~cs(age, df=0), nu.fo=~cs(age, df=0), family=BCCG, data=b)

# Visualize the LMS model: the centiles below correspond to -2, -1, mean, +1 and +2 SD curves and can be customized
centiles(LMS_obj, age, cent=c(2.275, 15.865, 50, 84.134, 97.725), main = "", ylab ="SHBG, nmol/L", xlab= "Age, years", box(lwd=2),
         lwd.centiles = 2, col.centiles = c("black", "black","red","black", "black"))

# If the centiles are too underfitted/overfitted, re-run the LMS above and customize degrees of freedom (df)
# Optimize model with respect to residual distribution and Q-test (see below)


# Quality control
summary(LMS_obj)                                                        # LMS model summary
rqres.plot(LMS_obj, plot.type = "all", type = "QQ")                     # Detrend Q-Q plot
wp(LMS_obj, xvar = age, n.inter = 2, xcut.points = c(8,10,12,14))       # Worm plot
Q.stats(LMS_obj, xvar=age,plot=TRUE)                                    # Q-test [Z1=M; Z2=S; Z3=L; Z4=Kurtosis should oscillate around 0 and in the range -2 to +2]


# Z-score residuals quality control
b$SHBG_sds  <- residuals(LMS_obj, what = "z-scores", type = "simple", terms=NULL)                         # retrieve z-scores for all observations in the LMS model
plot.default(b$age, b$SHBG_sds, xlab="Age, years", ylab="SHBG  z-scores")                                 # verify normal distribution of z-scores by age
hist(b$SHBG_sds, breaks=30, xlab="SHBG  z-scores", ylab="Frequency", main="Histogram of SHBG z-scores")   # verify normal distribution of z-scores in total 


# Extract model L,M,S parameters for the relevant age range 6-16 years
age <- seq(6,16,0.1)
L <- predict(LMS_obj, what="nu", type="response", newdata=data.frame(age = age))            # retrieve L parameter values in the model
M <- predict(LMS_obj, what="mu", type="response", newdata=data.frame(age = age))            # retrieve M parameter values in the model
S <- predict(LMS_obj, what="sigma", type="response", newdata=data.frame(age = age))         # retrieve S parameter values in the model
LMS_Data <- data.frame(age,L,M,S)                                                           
write.csv(LMS_Data, file = "Female_SHBG_LMS_for_age.csv")


# Extract model centile coordinates for age
newx <- seq(6,16,0.1)
mat2 <- data.frame(centiles.pred(LMS_obj, xname="age", xvalues=newx, type="standard-centiles"))
write.csv(mat2, file = "Female_SHBG_Centile_coordinates_for_age.csv")


# Plot the LMS reference curve model
ggplot(b,aes(age,hormone2)) + geom_point(col="gray50", size=1) + theme_bw() + labs(title="") + 
  scale_x_discrete(name ="Age, y", limits=c(6:16)) + theme(text = element_text(size=14)) + 
  scale_y_continuous(name ="SHBG, nmol/L", limits=c(0,200)) +
  geom_smooth(aes(x=mat2$age, y=mat2$`-2`), data=mat2, inherit.aes = FALSE, stat="identity", linetype="dashed", color = "red") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`-1`),data=mat2, inherit.aes = FALSE,stat="identity", linetype="dashed", color = "black") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`0`),data=mat2, inherit.aes = FALSE,stat="identity", linetype="solid", color = "black") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`1`),data=mat2, inherit.aes = FALSE,stat="identity", linetype="dashed", color = "black") + 
  geom_smooth(aes(x=mat2$age, y=mat2$`2`),data=mat2, inherit.aes = FALSE,stat="identity", linetype="dashed", color = "red") +
  annotate("text", x=17, y=mat2[101,"-2"], label = "-2 SD", color = "red") +
  annotate("text", x=17, y=mat2[101,"-1"] , label = "-1 SD", color = "black") +
  annotate("text", x=17, y=mat2[101,"0"] , label = "Mean", color = "black") +
  annotate("text", x=17, y=mat2[101,"1"] , label = "+1 SD", color = "black") +
  annotate("text", x=17, y=mat2[101,"2"] , label = "+2 SD", color = "red") +
  annotate("text", x=17.4, y=mat2[101,"2"] , label = "", color = "white")


#Print the ggplot figure (also supports .pdf and .jpeg)
ggsave("Female_SHBG_LMS_model.tiff", dpi = 1000, compression = "lzw")



### Other tools ###
### Calculate percentiles reference curves instead of SD curves
newx <- seq(6,16,0.1)  
mat3 <- centiles.pred(LMS_obj, type = "centiles", xname = "age", xvalues = newx, cent = c(2.5, 10, 25, 50, 75, 90, 97.5))
colnames(mat3) <- c("age","p2.5","p10","p25","p50","p75","p90","p97.5")

### Extract Z-scores for individual observations comprising the current LMS model
z <- residuals(LMS_obj, what = "z-scores", type = "simple", terms=NULL)

### Calculate Z-score from L,M,S entries for age
(((X/M)^L)-1)/(L*S) # X = blood sample result

### Calculate percentile coordinates using L,M,S entries for a certain age
(M*(1+(L*S*X))^(1/L)) # X = number of SDs from the mean, e.g. c(-2, -1, 0, 1, 2)

### Convert Z-scores to percentile scale
pnorm(2) 	 # +2 SD is equivalent to p97.72499
qnorm(0.9772499)	 # p97.72499 is equivalent to +2 SD 
