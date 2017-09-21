
# Firm equations ----------------------------------------------------------
#Endogenous proportion of healthy workers for no health insurance
#TODO: Change beliefs out of path when appropiate. Right now we are not interested in corner equilibria
#So Im taking these beliefs to be the same
Chi_0g = function(s,w0,w1,i){
  #Assign the same beliefs for the moment if there is no Labor supply
  if((L0_s('g',s,w0,w1) + L0_s('b',s,w0,w1))==0){
    aux = Chi_1g(s,w0,w1,i)
  }
  else{
    aux = (delta_sort(i)*L0_s('g',s,w0,w1)/(delta_sort(i)*L0_s('g',s,w0,w1) + L0_s('b',s,w0,w1))) 
  }
  return(aux)
}
#Endogenous proportion of healthy workers for health insurance
Chi_1g = function(s,w0,w1,i){
  #Assign the same beliefs for the moment if there is no Labor supply
  if((L1_s('g',s,w0,w1) + L1_s('b',s,w0,w1))==0){
    aux = Chi_0g(s,w0,w1,i)
  }
  else{
    aux = (delta_sort(i)*L1_s('g',s,w0,w1)/(delta_sort(i)*L1_s('g',s,w0,w1) + L1_s('b',s,w0,w1))) 
  }
  return(aux)
}
#Expected expenditure shock
#TODO: it has the analytical formula for exponential dist
E_m = function(h){
  if(h == 'g'){
    mu_h = mu_g
    sigma_h = sigma_g
    P_0h = P_0g
  }
  else{
    mu_h = mu_b
    sigma_h = sigma_b
    P_0h = P_0b
  }
  aux = (1-P_0h)*(mu_h + (sigma_h*dnorm(-mu_h/sigma_h)/(1-pnorm(-mu_h/sigma_h))))
  return(aux)
}
#Expected firm's medical expenditure
M = function(s,w0,w1,i){
  aux = (E_m('g')*Chi_1g(s,w0,w1,i)+E_m('b')*(1-Chi_1g(s,w0,w1,i)))/l1_s(w1)
  return(aux)
}
#Average labor productivity for contract without health insurance
gamma_prod_bar_0 = function(s,w0,w1,i){
  aux = s*gamma_prod(s,i)*((1-rho)*Chi_0g(s,w0,w1,i)+rho)
  return(aux)
}
#Average labor productivity for contract with health insurance
gamma_prod_bar_1 = function(s,w0,w1,i){
  aux = s*gamma_prod(s,i)*((1-rho)*Chi_1g(s,w0,w1,i)+rho)
  return(aux)
}
#Effective wage without health insurance
w_hat0 = function(s,w0,w1,i){
  aux = w0/gamma_prod_bar_0(s,w0,w1,i)
  return(aux)
}
#Effective wage with health insurance
w_hat1 = function(s,w0,w1,i){
  aux = (w1+M(s,w0,w1,i))/gamma_prod_bar_1(s,w0,w1,i)
  return(aux)
}
#Effective price of capital
R_hat = function(R,i){
  aux = R/z_prod(i)
  return(aux)
}
#Sufficient statistic for effective prices without insurance
p_0 = function(s,w0,w1,i){
  aux = (eta*w_hat0(s,w0,w1,i))/((1-eta)*psi)
  return(aux)
}
#Sufficient statistic for effective prices with insurance
p_1 = function(s,w0,w1,i){
  aux = (eta*w_hat1(s,w0,w1,i))/((1-eta)*psi)
  return(aux)
}
#Sufficient statistic for effective prices with capital
p_k = function(R,i){
  aux = (eta*R_hat(R,i))/((1-eta)*psi)
  return(aux)
}
#Conditional output function without health insurance
y_0 = function(s,w0,w1,i,Y){
  p0i = p_0(s,w0,w1,i)
  zetai = zeta_elas(i)
  w_hat0i = w_hat0(s,w0,w1,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(w_hat0i + psi*p0i^zetai))^sigma
  return(aux)
}
#Conditional output function with health insurance
y_1 = function(s,w0,w1,i,Y){
  p1i = p_1(s,w0,w1,i)
  zetai = zeta_elas(i)
  w_hat1i = w_hat1(s,w0,w1,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(w_hat1i + psi*p1i^zetai))^sigma
  return(aux)
}
#Conditional output function with capital
y_k = function(R,i,Y){
  pki = p_k(R,i)
  zetai = zeta_elas(i)
  R_hati = R_hat(R,i)
  aux = (Y*((sigma-1)/sigma)^sigma)*((B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))/(R_hati + psi*pki^zetai))^sigma
  return(aux)
}
#Conditional effective labor wihout health insurance 
l_hat0 = function(s,w0,w1,i,Y){
  zetai = zeta_elas(i)
  p0i = p_0(s,w0,w1,i)
  aux = y_0(s,w0,w1,i,Y)/(B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
}
#Conditional effective labor with health insurance 
l_hat1 = function(s,w0,w1,i,Y){
  zetai = zeta_elas(i)
  p1i = p_1(s,w0,w1,i)
  aux = y_1(s,w0,w1,i,Y)/(B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
}
#Conditional capital 
k = function(R,i,Y){
  zetai = zeta_elas(i)
  pki = p_k(R,i)
  aux = y_k(R,i,Y)/(B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates without health insurance
q_0 = function(s,w0,w1,i,Y){
  zetai = zeta_elas(i)
  p0i = p_0(s,w0,w1,i)
  aux = (y_0(s,w0,w1,i,Y)*p0i^zetai)/(B(i)*(eta*p0i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates with health insurance
q_1 = function(s,w0,w1,i,Y){
  zetai = zeta_elas(i)
  p1i = p_1(s,w0,w1,i)
  aux = (y_1(s,w0,w1,i,Y)*p1i^zetai)/(B(i)*(eta*p1i^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional intermediates with capital
q_k = function(R,i,Y){
  zetai = zeta_elas(i)
  pki = p_k(R,i)
  aux = (y_k(R,i,Y)*pki^zetai)/(B(i)*(eta*pki^(zetai-1)+(1-eta))^(zetai/(zetai-1)))
  return(aux)
} 
#Conditional Profit for Capital
Pi_k = function(R,Y,i){
  R_hati = R/z_prod(i)
  p_ki = (eta*R_hati)/((1-eta)*psi)
  aux =  (((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)*
    (((eta*p_ki^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/
       (R_hati+psi*p_ki^zeta_elas(i)))^(sigma-1) - C_A(i)
  return(aux)
}
#Conditional Profit for No Insurance
Pi_0 = function(s,w0,w1,R,Y,i){
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s_memo('g',s,w0,w1)
  L0_s_b_var = L0_s_memo('b',s,w0,w1)
  L1_s_g_var = L1_s_memo('g',s,w0,w1)
  L1_s_b_var = L1_s_memo('b',s,w0,w1)
  #Assigning the same beliefs if no labor supply
  if((L0_s_g_var + L0_s_b_var)==0){
    Chi_0gi = (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if((L1_s_g_var + L1_s_b_var)==0){
    Chi_0gi = (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = s*gamma_prod(s,i)*((1-rho)*Chi_0gi+rho)
  w_hat0i = w0/gamma_prod_bar_0i
  p_0i = (eta*w_hat0i)/((1-eta)*psi)
  ###
  aux = (((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)*
    (((eta*p_0i^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/
       (w_hat0i+psi*p_0i^zeta_elas(i)))^(sigma-1)
  return(aux)
}
#Conditional Profit for Insurance
Pi_1 = function(s,w0,w1,R,Y,i){
  #Do not call the same function more than one time if is not neccesary
  L0_s_g_var = L0_s_memo('g',s,w0,w1)
  L0_s_b_var = L0_s_memo('b',s,w0,w1)
  L1_s_g_var = L1_s_memo('g',s,w0,w1)
  L1_s_b_var = L1_s_memo('b',s,w0,w1)
  #Assigning the same beliefs if no labor supply
  if((L0_s_g_var + L0_s_b_var)==0){
    Chi_0gi = (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if((L1_s_g_var + L1_s_b_var)==0){
    Chi_0gi = (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_1i = s*gamma_prod(s,i)*((1-rho)*Chi_1gi+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = (E_mg*Chi_1gi+E_mb*(1-Chi_1gi))/l1_s(w1)
  w_hat1i = (w1+Mi)/gamma_prod_bar_1i
  p_1i = (eta*w_hat1i)/((1-eta)*psi)
  ###
  aux = (((B(i)*(sigma-1)/sigma)^(sigma-1))*Y/sigma)*
    (((eta*p_1i^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1)))/
       (w_hat1i+psi*p_1i^zeta_elas(i)))^(sigma-1)-C_IN
  return(aux)
}
###Functions with probabilities
ccp = function(w0H,w0L,w1H,w1L,R,Y,i){
  #Remember that if we want to vectorize this fn we need to change max by pmax
  Pi_k_val = Pi_k(R,Y,i)
  Pi_0H_val = Pi_0(sH,w0H,w1H,R,Y,i)
  Pi_0L_val = Pi_0(sL,w0L,w1L,R,Y,i)
  Pi_1H_val = Pi_1(sH,w0H,w1H,R,Y,i)
  Pi_1L_val = Pi_1(sL,w0L,w1L,R,Y,i)
  Pi_max = max(Pi_k_val,Pi_0H_val,Pi_0L_val,Pi_1H_val,Pi_1L_val)
  Pi_k_val_norm = (Pi_k_val-Pi_max)/s_ccp
  Pi_0H_val_norm = (Pi_0H_val-Pi_max)/s_ccp
  Pi_0L_val_norm = (Pi_0L_val-Pi_max)/s_ccp
  Pi_1H_val_norm = (Pi_1H_val-Pi_max)/s_ccp
  Pi_1L_val_norm = (Pi_1L_val-Pi_max)/s_ccp
  ccp_k = exp(Pi_k_val_norm)/(exp(Pi_k_val_norm)+exp(Pi_0H_val_norm)+exp(Pi_0L_val_norm)+exp(Pi_1H_val_norm)+exp(Pi_1L_val_norm))
  ccp_0H = exp(Pi_0H_val_norm)/(exp(Pi_k_val_norm)+exp(Pi_0H_val_norm)+exp(Pi_0L_val_norm)+exp(Pi_1H_val_norm)+exp(Pi_1L_val_norm))
  ccp_0L = exp(Pi_0L_val_norm)/(exp(Pi_k_val_norm)+exp(Pi_0H_val_norm)+exp(Pi_0L_val_norm)+exp(Pi_1H_val_norm)+exp(Pi_1L_val_norm))
  ccp_1H = exp(Pi_1H_val_norm)/(exp(Pi_k_val_norm)+exp(Pi_0H_val_norm)+exp(Pi_0L_val_norm)+exp(Pi_1H_val_norm)+exp(Pi_1L_val_norm))
  ccp_1L = exp(Pi_1L_val_norm)/(exp(Pi_k_val_norm)+exp(Pi_0H_val_norm)+exp(Pi_0L_val_norm)+exp(Pi_1H_val_norm)+exp(Pi_1L_val_norm))
  aux = c(ccp_k, ccp_0H, ccp_0L, ccp_1H, ccp_1L)
  return(aux)  
}
#Memosie version
ccp_memo = memoise(ccp)
