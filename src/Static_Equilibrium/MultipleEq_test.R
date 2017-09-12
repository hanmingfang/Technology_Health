# Multiple Equilibria -----------------------------------------------------

N_sims = 10
set.seed(3)
F_zero_mat = function(x){
  aux = matrix(nrow = N_sims, ncol = 6)
  for(i in 1:N_sims){
    aux[i,] = F_zeros(x[i,])
  }
  return(aux)
}
x_initial = matrix(runif(6*N_sims,min=0.1,max=5), N_sims, 6) # N initial guesses, each of length 6
#This code runs nleqslv iteratively for each initial guess
#So far I've found just one equilibrium
#TODO: Sometimes I get a problem of Non-finite value in he integral, check this maybe with tryCatch
ans = searchZeros(log(x_initial),F_zeros, method="Broyden",global="dbldog")
ans$x
w0H = exp(ans$x[1])
w0L = exp(ans$x[2])
w1H = exp(ans$x[3])
w1L = exp(ans$x[4])
R = exp(ans$x[5])
Y = exp(ans$x[6])
val = c(w0H,w0L,w1H,w1L,R,Y)
val
