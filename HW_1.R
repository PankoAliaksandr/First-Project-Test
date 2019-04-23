# Assignment 1: Minimum Variance Portfolios

# Glossary:
# PC - The Principal Component Estimator
# SAC - Shrinkage Towards the Average-Correlation Estimator 

library(zoo)
library(matlib)
library(matrixcalc)
library(quadprog)
library(ggplot2)
library(reshape2)

# Constants
number_of_stocks <- 500
number_of_months <- 60
one_month <- 1/12
last_available_date <- "Dec 2017"
five_years <- 5
pattern <- "%d-%m-%Y"
end_date <- as.yearmon("31-12-2017", pattern)

# Load provided data
load("D:/data_ex1_covariance_20181122.Rdata")

# Functions
prepare_data <- function(){

  # Create unique keys:
  
  # Delete leading zeros
  cnam$iid <<- sub("^[0]+", "", cnam$iid)
  DT.members$iid <<- sub("^[0]+", "", DT.members$iid)
  returns$iid <<- sub("^[0]+", "", returns$iid)
  
  # Unite gvkey and iid to create unique keys:
  cnam$id <<- paste(cnam$gvkey, cnam$iid, sep="")
  DT.members$id <<- paste(DT.members$gvkey, DT.members$iid, sep="")
  returns$id <<- paste(returns$gvkey, returns$iid, sep="")
  
  # Drop unrequired columns:
  cnam <<- subset(cnam, select = -c(iid, gvkey, sich))
  DT.members <<- subset(DT.members, select = -c(iid, gvkey))
  returns <<- subset(returns, select = -c(iid, gvkey))
}
get_relevant_companies <- function(end_date){
  
  # Find companies that were in index in June 2018 at least 1 day
  sp500_stocks <- DT.members[DT.members$ym == end_date,]
  
  # Reset index in the dataframe
  rownames(sp500_stocks) <- NULL
  
  # Add complany name
  sp500_stocks <- merge(x = sp500_stocks, y = cnam, by = "id", all.x = TRUE)
  
  # Drop unrequired columns: Date and id industry classification
  sp500_stocks <- subset(sp500_stocks, select = -c(ym))
  
  # Update number of stocks in the index
  number_of_stocks <<- nrow(sp500_stocks)
  
  return(sp500_stocks)
}
get_monthly_returns <- function(sp500_stocks, end_date){
  
  # Create a time period column
  time_period <- seq((end_date - five_years + one_month), end_date, by = one_month )
  
  monthly_returns <- data.frame(matrix(ncol = 1, nrow = number_of_months))
  colnames(monthly_returns) <- "Time Period" 
  monthly_returns["Time Period"] <- time_period
  
  # Create returns data frame size = 60x501
  for(stock_id in sp500_stocks$id){
    stock_returns_df <- returns[(returns$id == stock_id), c("ym", "trt1m")]
    colnames(stock_returns_df) <- c("Time Period", stock_id)
    monthly_returns <- merge(x = monthly_returns, y = stock_returns_df,
                             by = "Time Period", all.x = TRUE)
  }
  
  # Set index and drop date column
  rownames(monthly_returns) <- monthly_returns$`Time Period`
  monthly_returns$`Time Period`<- NULL
  
  # Fill NA with mean
  for(col in colnames(monthly_returns)){
    ratio_of_na <- sum(is.na(monthly_returns[[col]])) / length(monthly_returns[[col]])
    if(ratio_of_na > 0.2){
      monthly_returns[[col]] <- NULL
    }
    else{
      col_mean <- mean(monthly_returns[[col]], na.rm = TRUE)
      monthly_returns[[col]][is.na(monthly_returns[[col]])] <- col_mean
    }
  }
  number_of_stocks <<- min(number_of_stocks, ncol(monthly_returns))
  return(monthly_returns)
  
}
calculate_PC_covariance_matrix <- function(monthly_returns, sample_cov_monthly_returns){
  
  n <- 3
  
  # Eigenvalue Decomposition:
  # Vectors matrix here is orthogonal
  # Eigenvalues already sorted in descending order
  eigenS <- eigen(sample_cov_monthly_returns, symmetric = TRUE)
  
  Lambda <- diag(eigenS$values)
  diag(Lambda)[(n + 1):number_of_stocks] <- 0
  SigmaPC1 <- eigenS$vectors %*% Lambda %*% t(eigenS$vectors)
  
  # Find weights in a portfolio using formula in hint
  w1 <- eigenS$vectors[,1] / sum(eigenS$vectors[,1]) 
  w2 <- eigenS$vectors[,2] / sum(eigenS$vectors[,2])
  w3 <- eigenS$vectors[,3] / sum(eigenS$vectors[,3])
  
  # Determine corresponding returns using formula in hint
  r1 <- as.matrix(monthly_returns) %*% w1 
  r2 <- as.matrix(monthly_returns) %*% w2
  r3 <- as.matrix(monthly_returns) %*% w3
  
  resVar <- rep(0, times = number_of_stocks)
  
  for (j in 1:number_of_stocks){
    model <- lm(monthly_returns[,j] ~ (1 + r1 + r2 + r3))
    resVar[j] <- var(model$residuals)
  }
  
  SigmaPC2 <- diag(resVar)
  SigmaPC <- SigmaPC1 + SigmaPC2
  
  return(SigmaPC)
}
calculate_SAC_covariance_matrix <- function(monthly_returns, sample_cov_monthly_returns){
  
  sample_cor_monthly_returns <- cor(monthly_returns)
  
  mean_corr_value <- mean(sample_cor_monthly_returns[upper.tri(sample_cor_monthly_returns)])
  
  # 500x500
  C <- matrix(data = mean_corr_value,
              nrow = number_of_stocks,
              ncol = number_of_stocks)
  
  diag(C) <- 1
  
  # 500x500
  Delta <- diag(sqrt(diag(sample_cov_monthly_returns)))
  
  # 500x500
  SigmaAC <- Delta %*% C %*% Delta
  
  # 500x500
  SigmaSAC <- 0.5 * sample_cov_monthly_returns + 0.5 * SigmaAC
  
  return(SigmaSAC)
}
matrix_check <- function(matr){
  matr <- round(matr, 8)
  cat("Symmetric:", isSymmetric(matr))
  cat('\n')
  cat("Singular:", is.singular.matrix(matr, tol = 1e-8))
  cat('\n')
  cat("Positive semi-definite:", is.positive.semi.definite(matr, tol=1e-8))
  cat('\n')
}
check_covariance_matrices <- function(sample_cov_monthly_returns, SigmaSAC, SigmaPC){
  cat("\nFor Original covariance matrix:\n")
  matrix_check(sample_cov_monthly_returns)
  cat("\nFor SigmaPC covariance matrix:\n")
  matrix_check(SigmaPC)
  cat("\nFor SigmaSAC covariance matrix:\n")
  matrix_check(SigmaSAC)
}
analize_portfolio_weights <- function(portfolio_weights){
  portfolio_weights_rounded <- as.vector(round(portfolio_weights, digit = 4))
  
  cat("\nNumber of weights > 0:\n", sum(portfolio_weights_rounded > 0))
  cat("\nNumber of weights < 0:\n", sum(portfolio_weights_rounded < 0))
  cat("\nNumber of weights == 0:\n", sum(portfolio_weights_rounded == 0))
  
  return(portfolio_weights_rounded)
  
}
calculate_portfolio_weights <- function(SigmaPC, SigmaSAC){
  #  Minimum Variance Portfolio Construction 
  
  # Case 1: PC: No restriction (short sales allowed)
  
  # Find inverse matrix
  SigmaPC_inv <- solve(SigmaPC)
  
  # Creating unit-vector of length 500
  Unit <- c(rep(1,number_of_stocks)) 
  
  # Weight of stocks in minimum variance portfolio using the formula
  w_PC <- SigmaPC_inv %*% (Unit) / as.numeric(t(Unit) %*% SigmaPC_inv %*% Unit)
  
  # Case 2: PC: Restriction on short sales (short sales are not allowed)
  w_PC_c <- solve.QP( Dmat = SigmaPC,
                      dvec = rep(0, number_of_stocks),
                      Amat = cbind(rep(1, number_of_stocks), diag(number_of_stocks)),
                      bvec = c(1, rep(0, number_of_stocks)),
                      meq = 1)$solution
  
  # Case 3: SAC: No restriction (short sales allowed)
  
  # Find inverse matrix
  SigmaSAC_inv <- solve(SigmaSAC)
  
  # creating unit-vector of length 500
  Unit <- c(rep(1,number_of_stocks)) 
  
  # Weight of stocks in minimum variance portfolio using the formula
  w_SAC <- SigmaSAC_inv %*% (Unit) / as.numeric(t(Unit) %*% SigmaSAC_inv %*% Unit)
  
  # Case 4: SAC: Restriction on short sales (short sales are not allowed)
  w_SAC_c <- solve.QP( Dmat = SigmaSAC,
                       dvec = rep(0, number_of_stocks),
                       Amat = cbind(rep(1, number_of_stocks), diag(number_of_stocks)),
                       bvec = c(1, rep(0, number_of_stocks)),
                       meq = 1)$solution
  
  # Case 5: Equally weighted
  weight <- as.double(1/number_of_stocks)
  w_EQ <- rep(weight, times = number_of_stocks)
  
  weights <- data.frame("w_PC" = w_PC,
                        "w_PC_c" = w_PC_c,
                        "w_SAC" = w_SAC,
                        "w_SAC_c" = w_SAC_c,
                        "w_EQ" = w_EQ)
  
  return(weights)
}
analyze_portfolios_weights <- function(weights){
  # PC: No restrictions: How many of the stocks have positive, zero and negative weights?
  cat("\nFor PC non-restricted:\n")
  rounded_w_PC <- analize_portfolio_weights(weights$w_PC)
  # PC: Restricted: How many of the stocks have positive, zero and negative weights?
  cat("\nFor PC restricted:\n")
  rounded_w_PC_c <- analize_portfolio_weights(weights$w_PC_c)
  # SAC: No restrictions: How many of the stocks have positive, zero and negative weights?
  cat("\nFor SAC non-restricted:\n")
  rounded_w_SAC <- analize_portfolio_weights(weights$w_SAC)
  # SAC: Restricted: How many of the stocks have positive, zero and negative weights?
  cat("\nFor SAC restricted:\n")
  rounded_w_SAC_c <- analize_portfolio_weights(weights$w_SAC_c)
  
  rounded_weights <- data.frame("w_PC" = rounded_w_PC,
                        "w_PC_c" = rounded_w_PC_c,
                        "w_SAC" = rounded_w_SAC,
                        "w_SAC_c" = rounded_w_SAC_c)
  return(rounded_weights)
}
find_most_important_stocks <- function(sp500_stocks, rounded_weights) { 
 
  # PC: No restrictions: Output weights data frame 
  PC_weight_no_rest_df <- cbind(sp500_stocks, weights = rounded_weights$w_PC)
  # PC: Restricted: Output weights data frame 
  PC_weight_rest_df <- cbind(sp500_stocks, weights = rounded_weights$w_PC_c)
  # SAC: No restrictions: Output weights data frame 
  SAC_weight_no_rest_df <- cbind(sp500_stocks, weights = rounded_weights$w_SAC)
  # SAC: No restrictions: Output weights data frame 
  SAC_weight_rest_df <- cbind(sp500_stocks, weights = rounded_weights$w_SAC_c)
  
  ### SORT to find the most important ####
  
  # PC: No restrictions: Output weights data frame 
  PC_weight_no_rest_df_s <- PC_weight_no_rest_df[order(PC_weight_no_rest_df$weights,decreasing = TRUE),]
  # PC: Restricted: Output weights data frame 
  PC_weight_rest_df_s <- PC_weight_rest_df[order(PC_weight_rest_df$weights, decreasing = TRUE),]
  # SAC: No restrictions: Output weights data frame 
  SAC_weight_no_rest_df_s <- SAC_weight_no_rest_df[order(SAC_weight_no_rest_df$weights, decreasing = TRUE),]
  # SAC: No restrictions: Output weights data frame 
  SAC_weight_rest_df_s <- SAC_weight_rest_df[order(SAC_weight_rest_df$weights, decreasing = TRUE),]

  top_10 <- data.frame( "w_PC" = PC_weight_no_rest_df_s[1:10,],
                        "w_PC_c" = PC_weight_rest_df_s[1:10,],
                        "w_SAC" = SAC_weight_no_rest_df_s[1:10,],
                        "w_SAC_c" = SAC_weight_rest_df_s[1:10,])
  return(top_10)
} 
calculate_portfolios_returns <- function(monthly_returns, weights, end_date){
  # Calculate monthly return
  returns_PC_no_rest <- as.matrix(monthly_returns) %*% weights$w_PC
  returns_PC_rest <- as.matrix(monthly_returns) %*% weights$w_PC_c
  returns_SAC_no_rest <- as.matrix(monthly_returns) %*% weights$w_SAC
  returns_SAC_rest <- as.matrix(monthly_returns) %*% weights$w_SAC_c
  returns_EQ <- as.matrix(monthly_returns) %*% weights$w_EQ
  
  sp500_returns <- spx[spx$ym >= (end_date - five_years) &  spx$ym <= end_date]
  returns_SP500 <- sp500_returns$ret * 100 # in percents
  
  portfolios_returns <- data.frame( "PC" = returns_PC_no_rest,
                                    "PC_c" = returns_PC_rest,
                                    "SAC" = returns_SAC_no_rest,
                                    "SAC_c" = returns_SAC_rest,
                                    "EQ" = returns_EQ,
                                    "SP500" = returns_SP500)
  return(portfolios_returns)
}
calculate_avg_returns <- function(portfolios_returns){
  
  # Calculate average monthly return
  avg_returns_PC_no_rest <- mean(portfolios_returns$PC)
  avg_returns_PC_rest <- mean(portfolios_returns$PC_c)
  avg_returns_SAC_no_rest <- mean(portfolios_returns$SAC)
  avg_returns_SAC_rest <- mean(portfolios_returns$SAC_c)
  avg_returns_EQ <- mean(portfolios_returns$EQ)
  avg_returns_SP500 <- mean(portfolios_returns$SP500)
  
  avg_returns <- data.frame("PC" = avg_returns_PC_no_rest,
                            "PC_c" = avg_returns_PC_rest,
                            "SAC" = avg_returns_SAC_no_rest,
                            "SAC_c" = avg_returns_SAC_rest,
                            "EQ" = avg_returns_EQ,
                            "SP500" = avg_returns_SP500)
  return(avg_returns)
}
calculate_std_returns <- function(portfolios_returns){
  
  # Standard Deviations
  std_returns_PC_no_rest <- sd(portfolios_returns$PC)
  std_returns_PC_rest <- sd(portfolios_returns$PC_c)
  std_returns_SAC_no_rest <- sd(portfolios_returns$SAC)
  std_returns_SAC_rest <- sd(portfolios_returns$SAC_c)
  std_returns_EQ <- sd(portfolios_returns$EQ)
  std_returns_SP500 <- sd(portfolios_returns$SP500)
  
  std_returns <- data.frame("PC" = std_returns_PC_no_rest,
                            "PC_c" = std_returns_PC_rest,
                            "SAC" = std_returns_SAC_no_rest,
                            "SAC_c" = std_returns_SAC_rest,
                            "EQ" = std_returns_EQ,
                            "SP500" = std_returns_SP500)
  return(std_returns)
}
plot_average_returns <- function(avg_returns_df){
  
  column_names <- colnames(avg_returns)
  column_names <- c(column_names, "index")
  avg_returns_df$index = seq(1:516)
  colnames(avg_returns_df) <- column_names
  D = melt(avg_returns_df, id = 'index')
  ggplot(D, aes(index, value, group = variable, color = variable)) +
    geom_line() + ggtitle("Average Monthly Returns(%)") + xlab("Rolling window number") +
    ylab("Return")
  
}
plot_std_returns <- function(std_returns_df){
  
  column_names <- colnames(std_returns)
  column_names <- c(column_names, "index")
  std_returns_df$index = seq(1:516)
  colnames(std_returns_df) <- column_names
  R = melt(std_returns_df, id = 'index')
  ggplot(R, aes(index, value, group = variable, color = variable)) +
    geom_line() + ggtitle("Monthly Standard Deviations(%)") + xlab("Rolling window number") +
    ylab("Standard Deviation")
}
main_part_1 <- function(){
  prepare_data()
  sp500_stocks <- get_relevant_companies(end_date)
  monthly_returns <- get_monthly_returns(sp500_stocks, end_date)
  sample_cov_monthly_returns <- cov(monthly_returns)
  SigmaSAC <- calculate_SAC_covariance_matrix(monthly_returns, sample_cov_monthly_returns)
  SigmaPC <- calculate_PC_covariance_matrix(monthly_returns, sample_cov_monthly_returns)
  check_covariance_matrices(sample_cov_monthly_returns, SigmaSAC, SigmaPC)
  weights <- calculate_portfolio_weights(SigmaPC, SigmaSAC)
  rounded_weights <- analyze_portfolios_weights(weights)
  top_10 <- find_most_important_stocks(sp500_stocks, rounded_weights)
  portfolio_returns <- calculate_portfolios_returns(monthly_returns, weights, end_date)
  avg_returns <- calculate_avg_returns(portfolio_returns)
  std_returns <- calculate_std_returns(portfolio_returns)
}
main_part_2 <- function(){
  end_date <- DT.members[1]$ym + five_years
  i <- 1
  avg_returns_df <- data.frame(matrix(NA, nrow = 517, ncol = 6))
  std_returns_df <- data.frame(matrix(NA, nrow = 517, ncol = 6))
  
  while(end_date <= last_available_date){
    cat(i)
    start_date <- end_date - five_years
    
    sp500_stocks <- get_relevant_companies(end_date)
    monthly_returns <- get_monthly_returns(sp500_stocks, end_date)
    sample_cov_monthly_returns <- cov(monthly_returns)
    SigmaSAC <- calculate_SAC_covariance_matrix(monthly_returns, sample_cov_monthly_returns)
    SigmaPC <- calculate_PC_covariance_matrix(monthly_returns, sample_cov_monthly_returns)
    #check_covariance_matrices(sample_cov_monthly_returns, SigmaSAC, SigmaPC)
    weights <- calculate_portfolio_weights(SigmaPC, SigmaSAC)
    #rounded_weights <- analyze_portfolios_weights(weights)
    portfolio_returns <- calculate_portfolios_returns(monthly_returns, weights)
    avg_returns <- calculate_avg_returns(portfolio_returns)
    std_returns <- calculate_std_returns(portfolio_returns)
    
    avg_returns_df[i,]<- avg_returns
    std_returns_df[i,]<- std_returns
    
    i <- i + 1
    
    end_date <- end_date + one_month
  }

  plot_average_returns(avg_returns_df)
  plot_std_returns(std_returns_df)
  
}

main_part_1()
main_part_2()
