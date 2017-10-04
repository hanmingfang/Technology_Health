#Plots for comparative statics
dir = '~/Technology_Health/'
setwd(dir)
load("data/data.Rda")


data_aux = data.frame(nodes_seq, laguerre_vec, legendre_vec) 
ggplot(data_aux) + geom_line(aes(x=nodes_seq,y=laguerre_vec, colour = "Laguerre")) + 
  geom_line(aes(x=nodes_seq,y=legendre_vec, colour = "Legendre")) + 
  xlab("Number of nodes") + ylab("Integral value")
ggsave(file="numerical_integration_example.pdf", width=8, height=5)  