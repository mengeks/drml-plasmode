getPrimeNumbers <- function(n) {  
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
  {
    primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
    last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
  }
  which(primes)
}


library(doParallel)  
no_cores <- detectCores() - 1  


{
  cl <- makeCluster(no_cores, type="FORK")  
  registerDoParallel(cl)  
  ptm0 = proc.time()
  result <- foreach(i=10:10000) %dopar% getPrimeNumbers(i)
  stopCluster(cl) 
  ptm1 = proc.time()
  ptm1 - ptm0
}

{ 
  registerDoParallel(cores = no_cores)
  ptm0 = proc.time()
  result <- foreach(i=10:10000) %dopar% getPrimeNumbers(i)
  ptm1 = proc.time()
  ptm1 - ptm0
}


merge.by.time <- function(a, b) {
  merge(a, b, by='timestamp', suffixes=c('', ncol(a)))
}

results = foreach(i=1:5, .combine=merge.by.time) %dopar% {
  data.frame(timestamp=sample(1:10), feature=rnorm(10))
}

print(results)

result <- foreach(i=10:1000, .combine=function(a,b) rbind(a,b)) %dopar% {
   k <- i+1
   t <- i-1
   c(k,t)
}
result
matrix(unlist(result),ncol = 2)

do.call(rbind,lapply(file,function(x) x))

lapply(file,function(x) x)
rbind(c(1,2),c(3,4))

my_list <- list(l1 = c(1, 3, 5, 7),                 
                l2 = c(1, 2, 3),                     
                l3 = c(1, 1, 10, 5, 8, 65, 90))    

# Apply unlist R function 
print(unlist(my_list))       
