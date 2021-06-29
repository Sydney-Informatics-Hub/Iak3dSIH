
testf <- function(){
  x  <- 1
  y <- 2
  z <- 3
  return(list(x=x,y=y,z=z))
}

printnew <- function(parameters){
  print(parameters$data)
}

ModeltestOptions <- list(FitCubits=TRUE, LoadModel = list(a=1,b=2),data = mtcars, 
                    new=FALSE)
#ModeltestOptions$FitCubits
#ModeltestOptions$LoadModel
#ModeltestOptions$data

print(ModeltestOptions$LoadModel$a)
print(ModeltestOptions$data)