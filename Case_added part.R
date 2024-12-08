rm(list=ls())
library(vars)
library(xts)

#load the dataset
FRED_url <- url("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv")
FRED <- read.csv(FRED_url)

#extract the transformation codes
TCode_Ind <- FRED$INDPRO[1]
TCode_Cpi <- FRED$CPIAUCSL[1]
TCode_Fed <- FRED$FEDFUNDS[1]

#reduce the dataset until 2020
FRED_clean <- FRED[-1, ]
FRED_clean$sasdate <- as.Date(FRED_clean$sasdate, format = "%m/%d/%Y")
FRED_clean <- FRED_clean[FRED_clean$sasdate <= as.Date("2019-12-31"), ]


# Exercise 1 --------------------------------------------------------------

#plot the series without transformation
INDPRO_xts <- xts(FRED_clean$INDPRO, order.by = FRED_clean$sasdate)
CPIAUCSL_xts <- xts(FRED_clean$CPIAUCSL, order.by = FRED_clean$sasdate)
FEDFUNDS_xts <- xts(FRED_clean$FEDFUNDS, order.by = FRED_clean$sasdate)
plot(INDPRO_xts) 
plot(CPIAUCSL_xts)
plot(FEDFUNDS_xts)

#Transformation of the series 
Ind_stat <- diff(log(INDPRO_xts)) #growth rate
Cpi_stat <- diff(diff(log(CPIAUCSL_xts))) #change in growth rate
Fed_stat <- diff(FEDFUNDS_xts) #difference 
plot(Ind_stat) 
plot(Cpi_stat)
plot(Fed_stat)


# Exercise 2 --------------------------------------------------------------

#creates the dataframe for the VAR
var_data <- cbind(Ind_stat, Cpi_stat, Fed_stat)
var_data <- var_data[-(1:2)] #Removing first two rows as they will contain NULL due to taken differences etc.

#Popular approaches for getting the VAR order: (still has to be worked out)
# 1) Sequential tests on the nullity of coefficient matrices
# 2) Information criteria
# 3) Autocorrelation tests on the residuals

#select the number of lags with Information Criteria
lag_selection <- VARselect(var_data, lag.max = 20, type = "const")
print(lag_selection)

#estimate the VAR model with taking BIC as information criterion
lag_order <- 3  
var_model <- VAR(var_data, p = lag_order, type = "const")
summary(var_model)
#Model: Where IND is the growth rate of the IP index, the CPI is the change in the growth rate of CPI (all items), and the FED is the change in the effective federal funds rate
#Ind_stat = 0.00106 + 0.2898Ind_stat.l1 + 0.2283Cpi_stat.l1 + 0.0005901Fed_stat.l1 + 0.0716926Ind_stat.l2 + 0.0913Cpi_stat.l2 + 0.00046Fed_stat.l2 + 0.0999Ind_stat.l3 + 0.0774Cpi_stat.l3 -0.000696Fed_stat.l3 + error
#Cpi_stat = -2.448e-05 -0.0125Ind_stat.l1 -0.4407Cpi_stat.l1 + 5.506e-04*Fed_stat.l1 + 0.01633Ind_stat.l2 -0.3323Cpi_stat.l2 + 3.683e-04*Fed_stat.l2 + 0.011Ind_stat.l3 -0.2628Cpi_stat.l3 -1.899e-04*Fed_stat.l3 + error 
#Fed_stat = -0.03107+ 8.905Ind_stat.l1 + 17.404Cpi_stat.l1 + 0.396Fed_stat.l1 +5.871Ind_stat.l2 + 2.798Cpi_stat.l2 + -0.183Fed_stat.l2 + -0.305Ind_stat.l3 + 7.895Cpi_stat.l3 + -0.051Fed_stat.l3 + error 

#Validate the VAR...
#test to check for serial correlation in the residuals in a VAR model
serial_test <- serial.test(var_model, type = "BG")
print(serial_test)
plot(serial_test)


#Test for heteroscedasticity of residuals
arch_test <- arch.test(var_model, multivariate.only = FALSE)
print(arch_test)
plot(arch_test)

#Test Residuals for normality
normality_test <- normality.test(var_model,  multivariate.only = FALSE)
print(normality_test)
plot(normality_test)



# Exercise 3 --------------------------------------------------------------

#Explanation what Granger causality means...

#Granger tests
#Null: variable Y does NOT granger cause variable X for grangertest(X~Y)
#We reject this Null when the p-value of the grangertest is below 0.05
#We choose order=3 as we have a VAR(3)

Cpi_causes_Ind <- grangertest(Ind_stat ~ Cpi_stat, data = var_data, order = 3)
print(Cpi_causes_Ind) #NO REJECTION, SO NO GRANGER CAUSALITY
Ind_causes_Cpi <- grangertest(Cpi_stat ~ Ind_stat, data = var_data, order = 3)
print(Ind_causes_Cpi) #NO REJECTION, SO NO GRANGER CAUSALITY

Ind_causes_Fed <- grangertest(Fed_stat ~ Ind_stat, data = var_data, order = 3)
print(Ind_causes_Fed) #REJECTION, SO GRANGER CAUSALITY
Fed_causes_Ind <- grangertest(Ind_stat ~ Fed_stat, data = var_data, order = 3)
print(Fed_causes_Ind) #NO REJECTION, SO NO GRANGER CAUSALITY

Fed_causes_Cpi <- grangertest(Cpi_stat ~ Fed_stat, data = var_data, order = 3)
print(Fed_causes_Cpi) #REJECTION, SO GRANGER CAUSALITY
Cpi_causes_Fed <- grangertest(Fed_stat ~ Cpi_stat, data = var_data, order = 3)
print(Cpi_causes_Fed) #NO REJECTION, SO NO GRANGER CAUSALITY

#Hence, the growth rate of the IP index granger causes the change in the effective federal funds rate 
#And, the change in the effective federal funds rate granger causes the change in the growth rate of CPI (all items)


# Exercise 4 --------------------------------------------------------------

#What do Impulse Response functions measure?

#irf
irf_Cpi_Ind <- irf(var_model, impulse = "Ind_stat", response = "Cpi_stat", n.ahead = 10, , boot = TRUE)
irf_Ind_Cpi <- irf(var_model, impulse = "Cpi_stat", response = "Ind_stat", n.ahead = 10, , boot = TRUE)
irf_Fed_Ind <- irf(var_model, impulse = "Fed_stat", response = "Ind_stat", n.ahead = 10, , boot = TRUE)
irf_Ind_Fed <- irf(var_model, impulse = "Ind_stat", response = "Fed_stat", n.ahead = 10, , boot = TRUE)
irf_Cpi_Fed <- irf(var_model, impulse = "Cpi_stat", response = "Fed_stat", n.ahead = 10, , boot = TRUE)
irf_Fed_Cpi <- irf(var_model, impulse = "Fed_stat", response = "Cpi_stat", n.ahead = 10, , boot = TRUE)

plot(irf_Cpi_Ind)
plot(irf_Ind_Cpi)
plot(irf_Fed_Ind)
plot(irf_Ind_Fed)
plot(irf_Cpi_Fed)
plot(irf_Fed_Cpi)




