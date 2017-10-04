# Primitive Functions ---------------------------------------------------------------
#Distribution objects
#Distribution for Positive part of Medical expenditure
D_mg = Exp(rate_g)
D_mb = Exp(rate_b)
#Distribution for Risk aversion parameter
D_theta_gH = Gammad(scale = scale_gH,shape = shape_gH)
D_theta_gL = Gammad(scale = scale_gL,shape = shape_gL)
D_theta_bH = Gammad(scale = scale_bH,shape = shape_bH)
D_theta_bL = Gammad(scale = scale_bL,shape = shape_bL)
#Conditional cdf of Households' risk aversion parameter
F_gH  = function(theta){ 
  aux = p(D_theta_gH)(theta) #Expectation is shape*rate, Good health
  return(aux)
}
F_gL  = function(theta){ 
  aux = p(D_theta_gL)(theta) #Expectation is shape*rate, Good health
  return(aux)
} 
F_bH  = function(theta){
  aux = p(D_theta_bH)(theta) #Bad health
  return(aux)
}
F_bL  = function(theta){
  aux = p(D_theta_bL)(theta) #Bad health
  return(aux)
}
#Conditional cdf of Medical Expenditure (Positive part!)
H_g  = function(m){
  aux =  p(D_mg)(m) #Good health
  return(aux)
}     
H_b  = function(m){
  aux =  p(D_mb)(m) #Bad health
  return(aux)
}
#Parametrized functions
#Labor productivity (per skill unit)
gamma_prod = function(s,i){   
  if(s == sL){aux = A_0*exp(A*i)}
  else{aux = A_0*exp(A*i)*exp(delta_H*i)}
  return(aux)
}
#Sorting of workers
#Primitive function
sorting_f = function(i){
  #aux = 1
  #aux = i
  aux = exp(lambda_d*i)
  return(aux)
}
#Normalizing constant
norm_const_delta_sort = integrate(Vectorize(sorting_f),lower = N-1, upper = N)$value 
#Sorting function
delta_sort = function(i){
  aux = sorting_f(i)/norm_const_delta_sort
  return(aux)
}
#Capital productivity
z_prod = function(i){
  aux = 1
  return(aux)
}
#Function for elasticity of subst
zeta_elas = function(i){
  aux = zeta
  return(aux)
}
#Normalizing constant in production function
B = function(i){
  zetai = zeta_elas(i)
  aux = (1-eta)^(zetai/(1-zetai))
  return(aux)
}
#Automation Fixed Cost
C_A = function(i){
  aux = 0.1*exp(D*i)
  return(aux)
}
