#Plots for comparative statics. Run first "Comparative_statics.R" to update the data
dir = '~/Technology_Health/'
setwd(dir)
load("data/data.Rda")
dir = '~/Technology_Health/plots/comparative_statics_plots/'
setwd(dir)

param_changes =   c('D', 'lambda_d','rho', 'delta_H','lambda_H','lambda_g')
for(j in param_changes){
d = data_eq[[j]]
d$x = d$grid_par..j..
ggplot(d) + xlab(j) + ylab("Values") +
  geom_line(aes(x=d$x,y=d$w0H, colour = "w0H")) + 
  geom_line(aes(x=d$x,y=d$w1H, colour = "w1H"))
ggsave(file=paste(j,"wH",".pdf"), width=8, height=5)
ggplot(d) + xlab(j) + ylab("Values") +
  geom_line(aes(x=d$x,y=d$w0L, colour = "w0L")) + 
  geom_line(aes(x=d$x,y=d$w1L, colour = "w1L")) 
ggsave(file=paste(j,"wL",".pdf"), width=8, height=5)
ggplot(d) + xlab(j) + ylab("Values") +
  geom_line(aes(x=d$x,y=d$R, colour = "R"))
ggsave(file=paste(j,"R",".pdf"), width=8, height=5)
ggplot(d) + xlab(j) + ylab("Values") +
  geom_line(aes(x=d$x,y=d$m_k, colour = "m_k")) + 
  geom_line(aes(x=d$x,y=d$m_l0H, colour = "m_l0H")) + 
  geom_line(aes(x=d$x,y=d$m_l0L, colour = "m_l0L")) +  
  geom_line(aes(x=d$x,y=d$m_l1H, colour = "m_l1H")) + 
  geom_line(aes(x=d$x,y=d$m_l1L, colour = "m_l1L"))   
ggsave(file=paste(j,"measure_firms",".pdf"), width=8, height=5)
}
