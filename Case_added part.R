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
## Impulse Response Functions

# With the var package: irf() function
irf_cpi_ind <- irf(var_model, impulse = "Cpi_stat", response = "Ind_stat", n.ahead = 10, boot = TRUE)
irf_ind_ind <- irf(var_model, impulse = "Ind_stat", response = "Ind_stat", n.ahead = 10, boot = TRUE)
irf_fed_ind <- irf(var_model, impulse = "Fed_stat", response = "Ind_stat", n.ahead = 10, boot = TRUE)

irf_cpi_fed <- irf(var_model, impulse = "Cpi_stat", response = "Fed_stat", n.ahead = 10, boot = TRUE)
irf_fed_fed <- irf(var_model, impulse = "Fed_stat", response = "Fed_stat", n.ahead = 10, boot = TRUE)
irf_ind_fed <- irf(var_model, impulse = "Ind_stat", response = "Fed_stat", n.ahead = 10, boot = TRUE)

irf_cpi_cpi <- irf(var_model, impulse = "Cpi_stat", response = "Cpi_stat", n.ahead = 10, boot = TRUE)
irf_ind_cpi <- irf(var_model, impulse = "Ind_stat", response = "Cpi_stat", n.ahead = 10, boot = TRUE)
irf_fed_cpi <- irf(var_model, impulse = "Fed_stat", response = "Cpi_stat", n.ahead = 10, boot = TRUE)


# Plotting
library(ggplot2)

response_vars <- c("Ind_stat", "Fed_stat", "Cpi_stat") # Response variables
impulse_vars <- c("Cpi_stat", "Fed_stat", "Ind_stat")  # Impulse variables
all_fits <- list()

# Here we collect IRF fits for each response variable
for (response in response_vars) {
  fits <- lapply(impulse_vars, function(impulse) {
    irf(var_model, response = response, impulse = impulse,
        n.ahead = 10, boot = TRUE)
  })
  names(fits) <- impulse_vars
  all_fits[[response]] <- fits
}

plotdf_list <- list()

for (response in response_vars) {
  plotdf <- lapply(names(all_fits[[response]]), function(impulse) {
    data.frame(
      index = 1:nrow(all_fits[[response]][[impulse]]$irf[[1]]),
      value = all_fits[[response]][[impulse]]$irf[[1]][, 1],
      Lower = all_fits[[response]][[impulse]]$Lower[[1]][, 1],
      Upper = all_fits[[response]][[impulse]]$Upper[[1]][, 1],
      Impulse = impulse,
      Response = response
    )
  })
  plotdf_list[[response]] <- do.call(rbind, plotdf)
}

# We combine all data frames into one
final_plotdf <- do.call(rbind, plotdf_list)

# We create and save separate plots for each response variable
for (response in response_vars) {
  # Filter data for the current response variable
  response_plotdf <- final_plotdf[final_plotdf$Response == response, ]
  
  # Create the IRF plot for the current response variable
  p <- ggplot(response_plotdf, aes(x = index, y = value)) +
    geom_line() +
    facet_wrap(~Impulse, scales = "free_y") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = NA, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "red") +
    theme_bw() +
    labs(title = paste("Impulse Response Functions for", response),
         x = "Horizon",
         y = paste("Response", response))
  
  # Print the plot to view it
  print(p)
  
  # Save the plot as a file (optional)
  ggsave(filename = paste0("IRF_", response, ".png"), plot = p, width = 10, height = 6)
}


## Price Puzzle
plot(irf_fed_cpi)



## IRFs for the original time series in levels
# Starting levels
# Because when we compute the log differences, we effectively lose the 
# information about the original levels of the time series. To reconstruct
# the levels from the cumulative log differences, we need a reference point
# which is the last observed level from the original series before differencing.

starting_level_ind <- tail(as.numeric(INDPRO_xts), 1)
starting_level_cpi <- tail(as.numeric(CPIAUCSL_xts), 1)
starting_level_fed <- tail(as.numeric(FEDFUNDS_xts), 1)


#Response Variable: Industrial Production
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

# IND_stat → IND_stat
irf_values_ind_ind <- irf_ind_ind$irf$Ind_stat #impulse response coefficients
print(irf_values_ind_ind)  

cumulative_irf_ind_ind <- cumsum(irf_values_ind_ind) 
print(cumulative_irf_ind_ind)  

level_irf_ind_ind <- exp(cumulative_irf_ind_ind) * starting_level_ind
print(level_irf_ind_ind) 

plot(
  0:(length(level_irf_ind_ind) - 1), 
  level_irf_ind_ind, 
  type = "l", 
  col = "blue",
  xlab = "Lag", 
  ylab = "Level Response - IND_stat",
  main = "IND_stat -> IND_stat (levels)",
)
grid()

# CPI_stat → IND_stat
# Extract IRF values
irf_values_cpi_ind <- irf_cpi_ind$irf$Cpi_stat
print(irf_values_cpi_ind)

# Accumulate the IRFs
cumulative_irf_cpi_ind <- cumsum(irf_values_cpi_ind)
print(cumulative_irf_cpi_ind)

# Convert to levels
level_irf_cpi_ind <- exp(cumulative_irf_cpi_ind) * starting_level_ind
print(level_irf_cpi_ind)

# Plot
plot(
  0:(length(level_irf_cpi_ind) - 1), 
  level_irf_cpi_ind, 
  type = "l", 
  col = "blue",
  main = "CPI_stat -> IND_stat (levels)",
  xlab = "Lag", 
  ylab = "Level Response - IND_stat"
)
grid()  


# FED_stat → IND_stat
irf_values_fed_ind <- irf_fed_ind$irf$Fed_stat  
print(irf_values_fed_ind)

cumulative_irf_fed_ind <- cumsum(irf_values_fed_ind)  
print(cumulative_irf_fed_ind)

# Convert to levels
level_irf_fed_ind <- exp(cumulative_irf_fed_ind) * starting_level_ind
print(level_irf_fed_ind)

plot(
  0:(length(level_irf_fed_ind) - 1), 
  level_irf_fed_ind, 
  type = "l", 
  col = "blue",
  main = "FED_stat -> IND_stat (levels)",
  xlab = "Lag", 
  ylab = "Level Response - IND_stat"
)
grid()

mtext(
  "IRFs (levels) for Industrial Production (IND_stat)", 
  outer = TRUE, 
  cex = 1.5
)



# RESPONSE: COnsumer Price Index
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

# FED_stat → CPI_stat
irf_values_fed_cpi <- irf_fed_cpi$irf$Fed_stat  
print(irf_values_fed_cpi)

cumulative_irf_fed_cpi <- cumsum(irf_values_fed_cpi)  
print(cumulative_irf_fed_cpi)

level_irf_fed_cpi <- exp(cumulative_irf_fed_cpi) * starting_level_cpi
print(level_irf_fed_cpi)


plot(
  0:(length(level_irf_fed_cpi) - 1), 
  level_irf_fed_cpi, 
  type = "l", 
  col = "blue",
  main = "FED_stat -> CPI_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - CPI_stat"
)
grid()

# IND_stat → CPI_stat
irf_values_ind_cpi <- irf_ind_cpi$irf$Ind_stat  
print(irf_values_ind_cpi)

cumulative_irf_ind_cpi <- cumsum(irf_values_ind_cpi)  
print(cumulative_irf_ind_cpi)

level_irf_ind_cpi <- exp(cumulative_irf_ind_cpi) * starting_level_cpi
print(level_irf_ind_cpi)


plot(
  0:(length(level_irf_ind_cpi) - 1), 
  level_irf_ind_cpi, 
  type = "l", 
  col = "blue",
  main = "IND_stat -> CPI_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - CPI_stat"
)
grid()

# CPI_stat → CPI_stat
irf_values_cpi_cpi <- irf_cpi_cpi$irf$Cpi_stat  
print(irf_values_cpi_cpi)

cumulative_irf_cpi_cpi <- cumsum(irf_values_cpi_cpi)  
print(cumulative_irf_cpi_cpi)

level_irf_cpi_cpi <- exp(cumulative_irf_cpi_cpi) * starting_level_cpi
print(level_irf_cpi_cpi)

plot(
  0:(length(level_irf_cpi_cpi) - 1), 
  level_irf_cpi_cpi, 
  type = "l", 
  col = "blue",
  main = "CPI_stat -> CPI_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - CPI_stat"
)
grid()

mtext(
  "IRFs (levels) for Consumer Price Index (CPI_stat)", 
  outer = TRUE, 
  cex = 1.5
)


# Response Variable: Federal Funds Rates
# CPI_stat → FED_stat
irf_values_cpi_fed <- irf_cpi_fed$irf$Cpi_stat  
print(irf_values_cpi_fed)

cumulative_irf_cpi_fed <- cumsum(irf_values_cpi_fed)  
print(cumulative_irf_cpi_fed)

level_irf_cpi_fed <- exp(cumulative_irf_cpi_fed) * starting_level_fed
print(level_irf_cpi_fed)

plot(
  0:(length(level_irf_cpi_fed) - 1), 
  level_irf_cpi_fed, 
  type = "l", 
  col = "blue",
  main = "CPI_stat -> FED_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - FED_stat"
)
grid()

# IND_stat → FED_stat
irf_values_ind_fed <- irf_ind_fed$irf$Ind_stat  
print(irf_values_ind_fed)

cumulative_irf_ind_fed <- cumsum(irf_values_ind_fed)  
print(cumulative_irf_ind_fed)

level_irf_ind_fed <- exp(cumulative_irf_ind_fed) * starting_level_fed
print(level_irf_ind_fed)

plot(
  0:(length(level_irf_ind_fed) - 1), 
  level_irf_ind_fed, 
  type = "l", 
  col = "blue",
  main = "IND_stat -> FED_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - FED_stat"
)
grid()


# FED_stat → FED_stat
irf_values_fed_fed <- irf_fed_fed$irf$Fed_stat  
print(irf_values_fed_fed)

cumulative_irf_fed_fed <- cumsum(irf_values_fed_fed)  
print(cumulative_irf_fed_fed)

level_irf_fed_fed <- exp(cumulative_irf_fed_fed) * starting_level_fed
print(level_irf_fed_fed)

plot(
  0:(length(level_irf_fed_fed) - 1), 
  level_irf_fed_fed, 
  type = "l", 
  col = "blue",
  main = "FED_stat -> FED_stat (Levels)",
  xlab = "Lag", 
  ylab = "Level Response - FED_stat"
)
grid()

mtext(
  "IRFs (levels) for Federal Funds Rate (FED_stat)", 
  outer = TRUE, 
  cex = 1.5
)

