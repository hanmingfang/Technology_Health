# Excess demands with ccp -------------------------------------------------
#Using exponentials to ensure positivity of prices
#Excess demand for capital
k_excess_d_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  ###
  R_hati = function(i) R/z_prod(i)
  p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
  ###
  ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
  ###
  yki = function(i) {
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i) + psi*p_ki(i)^zeta_elas(i)))^sigma
    return(aux)
  }
  integrand_k = function(i) ccp_i(i)[1]*(yki(i)/(B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))
  integral = integrate(Vectorize(integrand_k), lower = N-1, upper = N) #Bounds specified in the document
  #We need to vectorize the function to use integrate
  aux = integral$value-K
  return(aux)
}
#Excess demand for labor without health insurance (High skill)
l0H_excess_d_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  L0_s_g_var = L0_s_memo('g',sH,w0H,w1H)
  L0_s_b_var = L0_s_memo('b',sH,w0H,w1H)
  L1_s_g_var = L1_s_memo('g',sH,w0H,w1H)
  L1_s_b_var = L1_s_memo('b',sH,w0H,w1H)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = function(i) sH*gamma_prod(sH,i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0H/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
  ###
  y0i = function(i){
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i) + psi*p_0i(i)^zeta_elas(i)))^sigma
    return(aux)
  }  
  #The integrand of l_0i should be l_hat_0i divided by gamma gamma_prod_bar_0i
  integrand_l0 = function(i) ccp_i(i)[2]*((y0i(i)/(B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_0i(i))
  integral = integrate(Vectorize(integrand_l0), lower = N-1, upper = N) #Bounds specified in the document, Vectorizing
  aux = integral$value - (L0_s_g_var + L0_s_b_var) #subtracting Labor supplied for no insurance
  return(aux)  
}
#Excess demand for labor without health insurance (Low skill)
l0L_excess_d_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  L0_s_g_var = L0_s_memo('g',sL,w0L,w1L)
  L0_s_b_var = L0_s_memo('b',sL,w0L,w1L)
  L1_s_g_var = L1_s_memo('g',sL,w0L,w1L)
  L1_s_b_var = L1_s_memo('b',sL,w0L,w1L)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_0i = function(i) sL*gamma_prod(sL,i)*((1-rho)*Chi_0gi(i)+rho)
  w_hat0i = function(i) w0L/gamma_prod_bar_0i(i)
  p_0i = function(i) (eta*w_hat0i(i))/((1-eta)*psi)
  ###
  ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
  ###
  y0i = function(i){
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0i(i) + psi*p_0i(i)^zeta_elas(i)))^sigma
    return(aux)
  }  
  #The integrand of l_0i should be l_hat_0i divided by gamma gamma_prod_bar_0i
  integrand_l0 = function(i) ccp_i(i)[3]*((y0i(i)/(B(i)*(eta*p_0i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_0i(i))
  integral = integrate(Vectorize(integrand_l0), lower = N-1, upper = N) #Bounds specified in the document, Vectorizing
  aux = integral$value - (L0_s_g_var + L0_s_b_var) #subtracting Labor supplied for no insurance
  return(aux)  
}
#Excess demand for labor with health insurance (High skill)
l1H_excess_d_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  L0_s_g_var = L0_s_memo('g',sH,w0H,w1H)
  L0_s_b_var = L0_s_memo('b',sH,w0H,w1H)
  L1_s_g_var = L1_s_memo('g',sH,w0H,w1H)
  L1_s_b_var = L1_s_memo('b',sH,w0H,w1H)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_1i = function(i) sH*gamma_prod(sH,i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1H)
  w_hat1i = function(i) (w1H+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  ###
  ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
  ###
  y1i = function(i){
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i) + psi*p_1i(i)^zeta_elas(i)))^sigma
    return(aux)
  }
  integrand_l1 = function(i) ccp_i(i)[4]*((y1i(i)/(B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_1i(i))
  integral = integrate(Vectorize(integrand_l1), lower = N-1, upper = N) #Bounds specified in the document, vectorize
  aux = integral$value - (L1_s_g_var + L1_s_b_var) #subtracting Labor supplied for insurance
  return(aux)  
}
#Excess demand for labor with health insurance (Low skill)
l1L_excess_d_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  L0_s_g_var = L0_s_memo('g',sL,w0L,w1L)
  L0_s_b_var = L0_s_memo('b',sL,w0L,w1L)
  L1_s_g_var = L1_s_memo('g',sL,w0L,w1L)
  L1_s_b_var = L1_s_memo('b',sL,w0L,w1L)
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_var)/(L0_s_g_var + L0_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  else if(is.na((L1_s_g_var)/(L1_s_g_var + L1_s_b_var))){
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var))
  }
  else{
    Chi_0gi = function(i) (delta_sort(i)*L0_s_g_var)/(delta_sort(i)*L0_s_g_var + L0_s_b_var)
    Chi_1gi = function(i) ((delta_sort(i)*L1_s_g_var)/(delta_sort(i)*L1_s_g_var + L1_s_b_var))
  }
  gamma_prod_bar_1i = function(i) sL*gamma_prod(sL,i)*((1-rho)*Chi_1gi(i)+rho)
  E_mg = E_m('g')
  E_mb = E_m('b')
  Mi = function(i) (E_mg*Chi_1gi(i)+E_mb*(1-Chi_1gi(i)))/l1_s(w1L)
  w_hat1i = function(i) (w1L+Mi(i))/gamma_prod_bar_1i(i)
  p_1i = function(i) (eta*w_hat1i(i))/((1-eta)*psi)
  ###
  ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
  ###
  y1i = function(i){
    aux = (Y*((sigma-1)/sigma)^sigma)*
      ((B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^
          (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1i(i) + psi*p_1i(i)^zeta_elas(i)))^sigma
    return(aux)
  }
  integrand_l1 = function(i) ccp_i(i)[5]*((y1i(i)/(B(i)*(eta*p_1i(i)^(zeta_elas(i)-1)+(1-eta))^(zeta_elas(i)/(zeta_elas(i)-1))))/gamma_prod_bar_1i(i))
  integral = integrate(Vectorize(integrand_l1), lower = N-1, upper = N) #Bounds specified in the document, vectorize
  aux = integral$value - (L1_s_g_var + L1_s_b_var) #subtracting Labor supplied for insurance
  return(aux)  
}
#Consistency of aggregate output
#TODO: check for this where to write the ccps
#TODO: Check again this function (is too long)
Y_excess_s_prob = function(p){
  w0H = exp(p[1])
  w0L = exp(p[2])
  w1H = exp(p[3])
  w1L = exp(p[4])
  R = exp(p[5])
  Y = exp(p[6])
  
  ###High skill
  L0_s_g_varH = L0_s_memo('g',sH,w0H,w1H)
  L0_s_b_varH = L0_s_memo('b',sH,w0H,w1H)
  L1_s_g_varH = L1_s_memo('g',sH,w0H,w1H)
  L1_s_b_varH = L1_s_memo('b',sH,w0H,w1H)
  
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_varH)/(L0_s_g_varH + L0_s_b_varH))){
    Chi_0giH = function(i) (delta_sort(i)*L1_s_g_varH)/(delta_sort(i)*L1_s_g_varH + L1_s_b_varH)
    Chi_1giH = function(i) ((delta_sort(i)*L1_s_g_varH)/(delta_sort(i)*L1_s_g_varH + L1_s_b_varH))
  }
  else if(is.na((L1_s_g_varH)/(L1_s_g_varH + L1_s_b_varH))){
    Chi_0giH = function(i) (delta_sort(i)*L0_s_g_varH)/(delta_sort(i)*L0_s_g_varH + L0_s_b_varH)
    Chi_1giH = function(i) ((delta_sort(i)*L0_s_g_varH)/(delta_sort(i)*L0_s_g_varH + L0_s_b_varH))
  }
  else{
    Chi_0giH = function(i) (delta_sort(i)*L0_s_g_varH)/(delta_sort(i)*L0_s_g_varH + L0_s_b_varH)
    Chi_1giH = function(i) ((delta_sort(i)*L1_s_g_varH)/(delta_sort(i)*L1_s_g_varH + L1_s_b_varH))
  }
  ###Low skill
  L0_s_g_varL = L0_s_memo('g',sL,w0L,w1L)
  L0_s_b_varL = L0_s_memo('b',sL,w0L,w1L)
  L1_s_g_varL = L1_s_memo('g',sL,w0L,w1L)
  L1_s_b_varL = L1_s_memo('b',sL,w0L,w1L)
  
  #Assigning the same beliefs if no labor supply
  if(is.na((L0_s_g_varL)/(L0_s_g_varL + L0_s_b_varL))){
    Chi_0giL = function(i) (delta_sort(i)*L1_s_g_varL)/(delta_sort(i)*L1_s_g_varL + L1_s_b_varL)
    Chi_1giL = function(i) ((delta_sort(i)*L1_s_g_varL)/(delta_sort(i)*L1_s_g_varL + L1_s_b_varL))
  }
  else if(is.na((L1_s_g_varL)/(L1_s_g_varL + L1_s_b_varL))){
    Chi_0giL = function(i) (delta_sort(i)*L0_s_g_varL)/(delta_sort(i)*L0_s_g_varL + L0_s_b_varL)
    Chi_1giL = function(i) ((delta_sort(i)*L0_s_g_varL)/(delta_sort(i)*L0_s_g_varL + L0_s_b_varL))
  }
  else{
    Chi_0giL = function(i) (delta_sort(i)*L0_s_g_varL)/(delta_sort(i)*L0_s_g_varL + L0_s_b_varL)
    Chi_1giL = function(i) ((delta_sort(i)*L1_s_g_varL)/(delta_sort(i)*L1_s_g_varL + L1_s_b_varL))
  }  
  
  Y_s = function(Y){
    R_hati = function(i) R/z_prod(i)
    p_ki = function(i) (eta*R_hati(i))/((1-eta)*psi)
    ###
    ccp_i = function(i) ccp_memo(w0H,w0L,w1H,w1L,R,Y,i)
    ###
    yki = function(i) {
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_ki(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(R_hati(i) + psi*p_ki(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    ###No insurance
    ###High skill
    gamma_prod_bar_0iH = function(i) sH*gamma_prod(sH,i)*((1-rho)*Chi_0giH(i)+rho)
    w_hat0iH = function(i) w0H/gamma_prod_bar_0iH(i)
    p_0iH = function(i) (eta*w_hat0iH(i))/((1-eta)*psi)
    y0iH = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_0iH(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0iH(i) + psi*p_0iH(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    ###Low skill
    gamma_prod_bar_0iL = function(i) sL*gamma_prod(sL,i)*((1-rho)*Chi_0giL(i)+rho)
    w_hat0iL = function(i) w0L/gamma_prod_bar_0iL(i)
    p_0iL = function(i) (eta*w_hat0iL(i))/((1-eta)*psi)
    y0iL = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_0iL(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat0iL(i) + psi*p_0iL(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    ###With insurance
    ###High skill
    gamma_prod_bar_1iH = function(i) sH*gamma_prod(sH,i)*((1-rho)*Chi_1giH(i)+rho)
    E_mg = E_m('g')
    E_mb = E_m('b')
    MiH = function(i) (E_mg*Chi_1giH(i)+E_mb*(1-Chi_1giH(i)))/l1_s(w1H)
    w_hat1iH = function(i) (w1H+MiH(i))/gamma_prod_bar_1iH(i)
    p_1iH = function(i) (eta*w_hat1iH(i))/((1-eta)*psi)
    y1iH = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_1iH(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1iH(i) + psi*p_1iH(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    ###Low skill
    gamma_prod_bar_1iL = function(i) sL*gamma_prod(sL,i)*((1-rho)*Chi_1giL(i)+rho)
    E_mg = E_m('g')
    E_mb = E_m('b')
    MiL = function(i) (E_mg*Chi_1giL(i)+E_mb*(1-Chi_1giL(i)))/l1_s(w1L)
    w_hat1iL = function(i) (w1L+MiL(i))/gamma_prod_bar_1iL(i)
    p_1iL = function(i) (eta*w_hat1iL(i))/((1-eta)*psi)
    y1iL = function(i){
      aux = (Y*((sigma-1)/sigma)^sigma)*
        ((B(i)*(eta*p_1iL(i)^(zeta_elas(i)-1)+(1-eta))^
            (zeta_elas(i)/(zeta_elas(i)-1)))/(w_hat1iL(i) + psi*p_1iL(i)^zeta_elas(i)))^sigma
      return(aux)
    }
    ###Computing the integral with expected output
    integrand = function(i) (ccp_i(i)[1]*yki(i)+ccp_i(i)[2]*y0iH(i)+
                               ccp_i(i)[3]*y0iL(i)+ccp_i(i)[4]*y1iH(i)+
                               ccp_i(i)[5]*y1iL(i))^((sigma-1)/sigma) #integrand of y_1 in the consumption good production
    integral = integrate(Vectorize(integrand), lower = N-1, upper = N) #Vectorizing
    aux = integral$value
    return(aux)
  }
  aux = Y - (Y_s(Y))^(sigma/(sigma-1))
  return(aux) 
}
#Excess demand funcns stacked
F_zeros = function(p){
  aux = c(k_excess_d_prob(p), l0H_excess_d_prob(p), l0L_excess_d_prob(p), 
          l1H_excess_d_prob(p), l1L_excess_d_prob(p), Y_excess_s_prob(p))
  return(aux)
}