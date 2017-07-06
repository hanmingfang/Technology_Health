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
eta = 0.5
psi = 1
w_0 = 10
w_1 = 9.46
theta = 0.5
mu = 50
P = 0
sd = 30
a = 0
b = 100
# Functions ---------------------------------------------------------------
#With truncated Normal
U_0 = function(theta,w_0,psi,eta){
  l_0 = (w_0/psi)^(1/eta)
 #the function of the integrand must be written in vector form 
  u_0 = function(x){
    aux = vector(length = length(x))
    for(i in 1:length(x)){
      pdf_x = dtruncnorm(x[i], a = a, b = b, mean = mu, sd = sd)
      aux[i] = (1/(1-theta))*((w_0*l_0 - x[i] - 
                               psi*((l_0^(1+eta))/(1+eta)))^(1-theta)-1)*pdf_x
    }
    return(aux)
  }
  
  integrand = u_0
  integral = integrate(integrand, lower = 0, upper = 100)
  return(integral$value)
}

U_1 = function(theta,w_1,psi,eta, P){
  l_1 = (w_1/psi)^(1/eta)
  aux = (1/(1-theta))*((w_1*l_1 - P - psi*((l_1^(1+eta))/(1+eta)))^(1-theta)-1)
  return(aux)

}
#Finding the willingness to pay for Health insurance by equating expected 
#utilities (see the document)
P_will = function(theta,w_0,w_1,psi,eta){
  l_0 = (w_0/psi)^(1/eta)
  l_1 = (w_1/psi)^(1/eta)
  #These values must be admitable by the function so I set such that argument = 0
  # and a arbitrary big negative number in the other case
  initial = -1e6
  final = w_1*l_1 - psi*((l_1^(1+eta))/(1+eta))
  #Function
  fun = function (P){
    return(U_1(theta,w_1,psi,eta, P)-U_0(theta,w_0,psi,eta))
  } 
  aux = uniroot(fun, c(initial,final), tol = 1e-13)$root
  return(aux)
} 

# Plots -------------------------------------------------------------------

# Willingness to pay ------------------------------------------------------

#Varying theta between 0 and 1
#Update the directory every time
dir = '~/Dropbox/Technology/Documents/Notes Model 06-05-17/figures/'
setwd(dir)
theta_grid = seq(from = 0, to = 0.99, by = 0.01)
will_grid = vector(length = length(theta_grid))

for(i in 1:length(theta_grid)){
  theta = theta_grid[i]
  will_grid[i] = P_will(theta,w_0,w_1,psi,eta)
}
#Transform to dataframe
d = as.data.frame(cbind(theta_grid, will_grid))
#Plot
ggplot(d, aes(x = theta_grid)) +
  geom_line(aes(y = will_grid, colour = "Willigness to pay")) + 
  xlab("Theta") +
  ylab("Willigness to pay") +
ggsave(file="will_pay_neg.pdf", width=8, height=5)





