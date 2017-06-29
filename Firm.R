#Solve the firm's problem for given wages and endogenous measures

# Libraries ---------------------------------------------------------------

library(stargazer, quietly = TRUE)
library(lmtest, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(dplyr)
library(reshape2)
library(ggfortify)
library(data.table)
library(RMySQL)
library(plotly)
library(latex2exp)
library(truncnorm)
Sys.setenv("plotly_username"="tlarrouc")
Sys.setenv("plotly_api_key"="bdqUsPxMwTYk9lI7KmY4")
dir = '~/Dropbox/Technology/Codes/Robots/'
setwd(dir)


# Parameters --------------------------------------------------------------
#This is just to see the shape of the policy functions
eta = 0.5
zeta = 0.5
gamma = 5 #This is not a parameter, but just to try
rho = 0.8
chi0 = 0.5
Y = 1
sigma = 0.5
psi = 0.1
w0 = 0.1
w1 = 0.05


# Functions ---------------------------------------------------------------

#Firm i>I* (non-automated task)

y0 = function(eta, zeta, gamma, chi0, rho, inputs){
  q = inputs[1]
  l = inputs[2]
  B = (1-eta)^(zeta/(1-zeta))
  aux = B*(eta*q^((zeta-1)/zeta) + (1-eta)*(gamma*l*chi0 + rho*gamma*l*(1-chi0))^((zeta-1)/zeta))^(zeta/(zeta-1))
  return(aux)  
}

Profits_0 = function(Y, sigma, psi, w1, w0, eta, zeta, gamma, chi0, inputs){
  q = inputs[1]
  l = inputs[2]
  y = y0(eta, zeta, gamma, chi0, rho, inputs)
  aux = Y^(1/sigma)*y^(1-1/sigma) - psi*q - w0*l
  return(aux)
}

y1 = function(eta, zeta, gamma, chi1, rho, inputs){
  q = inputs[1]
  l = inputs[2]
  B = (1-eta)^(zeta/(1-zeta))
  aux = B*(eta*q^((zeta-1)/zeta) + (1-eta)*(gamma*l*chi1 + rho*gamma*l*(1-chi1))^((zeta-1)/zeta))^(zeta/(zeta-1))
  return(aux)  
}

# Simulation --------------------------------------------------------------

opt_object = optim(par=c(0.3,0.2), fn = Profits_0, Y = Y, sigma = sigma, psi = psi,
                   w1 = w1, w0 = w0, eta = eta, zeta = zeta, gamma = gamma, chi0 = chi0,
                   control = list(fnscale = -1), method = "L-BFGS-B")
inputs_opt = opt_object$par
q_opt = inputs_opt[1]
l = inputs_opt[2]
Max_Profits0 = opt_object$value














