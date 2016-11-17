dpostg <- function(y,x,z,betas,gammas,gpri,Gpri,model,m) {
  L <- veros(y,x,z,betas,gammas,model,m)
  P <- dmvnorm(t(gammas), gpri, Gpri) # priori
  value <- L*P
  value
}