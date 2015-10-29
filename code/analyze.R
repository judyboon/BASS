n = 100
p = 500
output = "../output/MCMC0PX20ik15n100/"

# read estimated latent variable
X = scan(paste(output,'ETA',sep="")); 

# determine the number of factors left
k = length(X)/n;
X = matrix(X,k,n,byrow=T);

# read estimated loading
LAMBDA = scan(paste(filedir,'LAMBDA',sep=""));
Lambda = matrix(LAMBDA,p,k,byrow=T);

# read indicator variable Z
Z = scan(paste(filedir,'ZS',sep=""));
Z = matrix(Z,nrow=v,byrow=T);

# read error variance
Sigma2inv = scan(paste(filedir,'SIGMA2INV',sep=""));
Sigma2 = diag(1/Sigma2inv)
