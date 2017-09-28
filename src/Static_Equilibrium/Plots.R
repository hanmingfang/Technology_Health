# Plots -------------------------------------------------------------------
#Some plots to check some results
#Most of the plots are just for the deterministic version of the model
dir = '~/Technology_Health/plots/'
setwd(dir)

#Plot to show that FOSD in Assumption 1 holds for this case
ggplot(data.frame(x=c(0, 4)), aes(x=x)) + 
  stat_function(fun=H_g, geom="line", aes(colour = "H_g")) + xlab("x") + 
  ylab("y") + stat_function(fun=H_b, geom="line",aes(colour = "H_b"))
ggsave(file="H_FOSD.pdf", width=8, height=5)

#Plot to show that FOSD in Assumption 2 holds for this case
ggplot(data.frame(x=c(0, 4)), aes(x=x)) + xlab("x") +  ylab("y") +
  stat_function(fun=F_gH, geom="line", aes(colour = "F_gH")) +  
  stat_function(fun=F_gL, geom="line", aes(colour = "F_gL")) + 
  stat_function(fun=F_bH, geom="line",aes(colour = "F_bH")) +
  stat_function(fun=F_bL, geom="line",aes(colour = "F_bL"))
ggsave(file="F_FOSD.pdf", width=8, height=5)

#Plot gamma_prod
gamma_prodH_plot = function(i) gamma_prod(sH,i)
gamma_prodL_plot = function(i) gamma_prod(sL,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x))+ xlab("i") + ylab("") +
  stat_function(fun=gamma_prodH_plot, geom="line", aes(colour = "gamma_prodH")) + 
  stat_function(fun=gamma_prodL_plot, geom="line", aes(colour = "gamma_prodL")) 
ggsave(file="gamma_prod.pdf", width=8, height=5)

#Plot delta_sort
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") + ylab("") +
  stat_function(fun=delta_sort, geom="line") 
ggsave(file="delta_sort.pdf", width=8, height=5)

#Plot C_A
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + 
  stat_function(fun=C_A , geom="line") + xlab("i") + ylab("")
ggsave(file="automation_cost.pdf", width=8, height=5)

#Plot Chi_0g and Chi_1g  
Chi_0gH_plot = function(i) Chi_0g(sH, w0=w0H, w1=w1H,i)
Chi_0gL_plot = function(i) Chi_0g(sL, w0=w0L, w1=w1L,i)
Chi_1gH_plot = function(i) Chi_1g(sH, w0=w0H, w1=w1H,i)
Chi_1gL_plot = function(i) Chi_1g(sL, w0=w0L, w1=w1L,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") + ylab("") +
  stat_function(fun=Chi_0gH_plot, geom="line", aes(colour = "Chi_0gH")) + 
  stat_function(fun=Chi_0gL_plot, geom="line", aes(colour = "Chi_0gL")) + 
  stat_function(fun=Chi_1gH_plot, geom="line",aes(colour = "Chi_1gH")) +
  stat_function(fun=Chi_1gL_plot, geom="line",aes(colour = "Chi_1gL"))
ggsave(file="endogenous_proportion_healthy.pdf", width=8, height=5)

#Plot difference for adverse selection
Delta_Chi_gH = function(i) Chi_0gH_plot(i) - Chi_1gH_plot(i)
Delta_Chi_gL = function(i) Chi_0gL_plot(i) - Chi_1gL_plot(i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") + ylab("") +
  stat_function(fun=Delta_Chi_gH, geom="line", aes(colour = "Delta_Chi_gH")) + 
  stat_function(fun=Delta_Chi_gL, geom="line", aes(colour = "Delta_Chi_gL")) +
  ggtitle("Degree of Adverse selection conditional on Skill type")
ggsave(file="delta_endogenous_proportion_healthy.pdf", width=8, height=5)

#PlotEvolution of Expected medical expenditure across i 
M_H_plot = function(i) M(sH, w0H, w1H,i)
M_L_plot = function(i) M(sL, w0L, w1L,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") +  ylab("") +
  stat_function(fun=M_H_plot, geom="line", aes(colour = "M_H")) +
  stat_function(fun=M_L_plot, geom="line", aes(colour = "M_L")) 
ggsave(file="expected_medical_expenditure_skill_type.pdf", width=8, height=5)

#Plot Labor average productivity
gamma_prod_bar_0H_plot = function(i) gamma_prod_bar_0(sH, w0H, w1H,i)
gamma_prod_bar_0L_plot = function(i) gamma_prod_bar_0(sL, w0L, w1L,i)
gamma_prod_bar_1H_plot = function(i) gamma_prod_bar_1(sH, w0H, w1H,i)
gamma_prod_bar_1L_plot = function(i) gamma_prod_bar_1(sL, w0L, w1L,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") +  ylab("") + 
  stat_function(fun=gamma_prod_bar_0H_plot, geom="line", aes(colour = "gamma_bar0H")) + 
  stat_function(fun=gamma_prod_bar_1H_plot, geom="line", aes(colour = "gamma_bar1H")) +
  stat_function(fun=gamma_prod_bar_0L_plot, geom="line", aes(colour = "gamma_bar0L")) + 
  stat_function(fun=gamma_prod_bar_1L_plot, geom="line", aes(colour = "gamma_bar1L"))
ggsave(file="average_labor_productivity.pdf", width=8, height=5)

#Plot effective wages and prices
w_hat0H_plot = function(i) w_hat0(sH, w0H, w1H,i)
w_hat0L_plot = function(i) w_hat0(sL, w0L, w1L,i)
w_hat1H_plot = function(i) w_hat1(sH, w0H, w1H,i)
w_hat1L_plot = function(i) w_hat1(sL, w0L, w1L,i)
R_hat_plot = function(i) R_hat(R,i)
ggplot(data.frame(x=c(N-1,2)), aes(x=x)) + xlab("i") +  ylab("") + 
  stat_function(fun=w_hat0H_plot, geom="line", aes(colour = "w_hat0H")) + 
  stat_function(fun=w_hat1H_plot, geom="line", aes(colour = "w_hat1H")) +
  stat_function(fun=w_hat0L_plot, geom="line", aes(colour = "w_hat0L")) + 
  stat_function(fun=w_hat1L_plot, geom="line", aes(colour = "w_hat1L")) +
  stat_function(fun=R_hat_plot, geom="line", aes(colour = "R_hat"))
ggsave(file="effective_wages_equilibrium.pdf", width=8, height=5)

#Plot ccp
ccp_k_plot = function(i){
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
  aux = ccp_k
  return(aux)  
}
ccp_0H_plot = function(i){
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
  aux = ccp_0H
  return(aux)  
}
ccp_0L_plot = function(i){
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
  aux = ccp_0L
  return(aux)  
}
ccp_1H_plot = function(i){
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
  aux = ccp_1H
  return(aux)  
}
ccp_1L_plot = function(i){
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
  aux = ccp_1L
  return(aux)  
}
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + ylab("") + xlab("i") +
  stat_function(fun=ccp_k_plot, geom="line", aes(colour = "ccp_k")) + 
  stat_function(fun=ccp_0H_plot, geom="line", aes(colour = "ccp_0H")) +
  stat_function(fun=ccp_0L_plot, geom="line", aes(colour = "ccp_0L")) +
  stat_function(fun=ccp_1H_plot, geom="line", aes(colour = "ccp_1H")) +
  stat_function(fun=ccp_1L_plot, geom="line", aes(colour = "ccp_1L")) +
  ggtitle(paste("(w0H,w0L,w1H,w1L,R,Y) = (",round(w0H,2),",",round(w0L,2),",",
                round(w1H,2),",",round(w1L,2),",",round(R,2),",",
                round(Y,2),") \n","(w0,w1) = ","(",round(w0,2),",",round(w1,1),")"))
ggsave(file="ccp.pdf", width=8, height=5)

#Plot equilibrium conditional profits (deterministic part)
Pi_k_plot = function(i) Pi_k(R,Y,i)
Pi_0L_plot = function(i) Pi_0(sL,w0L,w1L,R,Y,i)
Pi_0H_plot = function(i) Pi_0(sH,w0H,w1H,R,Y,i)
Pi_1L_plot = function(i) Pi_1(sL,w0L,w1L,R,Y,i)
Pi_1H_plot = function(i) Pi_1(sH,w0H,w1H,R,Y,i)
ggplot(data.frame(x=c(N-1,N)), aes(x=x)) + xlab("i") + ylab("") +
  stat_function(fun = Vectorize(Pi_k_plot), geom="line", aes(colour = "Pi_k")) + 
  stat_function(fun = Vectorize(Pi_0H_plot), geom="line", aes(colour = "Pi_0H")) +
  stat_function(fun = Vectorize(Pi_1H_plot), geom="line", aes(colour = "Pi_1H")) +
  stat_function(fun = Vectorize(Pi_0L_plot), geom="line", aes(colour = "Pi_0L")) +
  stat_function(fun = Vectorize(Pi_1L_plot), geom="line", aes(colour = "Pi_1L")) +
  ggtitle(paste("(w0H,w0L,w1H,w1L,R,Y) = (",round(w0H,2),",",round(w0L,2),",",
                round(w1H,2),",",round(w1L,2),",",round(R,2),",",
                round(Y,2),") \n","(w0,w1) = ","(",round(w0,2),",",round(w1,1),")"))
ggsave(file="conditional_profits_prob_equilibrium.pdf", width=8, height=5)

#Plot theta_ins with MGF
w0H = 104
w1H = 4
Delta_w = w0H-w1H
#Delta_w = w0L-w1L

theta_ins_MGF_g_w_approx = function(theta){
  rate_h = rate_g
  P_0h = P_0g
  if(theta!= rate_h){
    if(exp((theta-rate_h)*M_trunc)==Inf){ #Doing a Taylor approximation in this case
      aux = ((theta-rate_h)*M_trunc +log(rate_h*(1-P_0h)) -log((1-exp(-rate_h*M_trunc))*(theta-rate_h)))/theta - Delta_w
    }
    else{
      aux = log(rate_h*(exp((theta-rate_h)*M_trunc)-1)*(1-P_0h)/((1-exp(-rate_h*M_trunc))*(theta-rate_h)) + P_0h)/theta - Delta_w
    }
  }
  else{
    aux = log(rate_h*M_trunc*(1-P_0h)/(1-exp(-rate_h*M_trunc)) + P_0h)/rate_h - Delta_w
  }
  return(aux)
}

ggplot(data.frame(x=c(0.001,10)), aes(x=x)) + xlab("theta") + ylab("") +
  stat_function(fun = Vectorize(theta_ins_MGF_g_w_approx), geom="line", aes(colour = "P(theta)")) 
#  stat_function(fun = theta_ins_MGF_b, geom="line", aes(colour = "Theta_MGF_b")) 
ggsave(file="P(theta).pdf", width=8, height=5)

















