dpostb <- function(y,x,z,betas,gammas,bpri,Bpri,model,m) {
  
  L <- veros(y,x,z,betas,gammas,model,m)
  P <- dmvnorm(t(betas), bpri, Bpri) # priori
  value <- L*P
  value
}