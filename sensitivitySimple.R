#index - numerical index denoting matched set membership
#q - vector q formed based upon observed responses to form test statistic t = Z'q
#Z - treatment indicator

sensitivitySimple = function(index, q, Z, alpha = 0.05, alternative = "two.sided", Gamma=1, calculatePvalue = T, continuousRelax = F, maximizeGamma = T)
{
	gur = suppressWarnings(require(gurobi))
	require(Matrix)
	Z = 1*Z
	ns = table(index)
	ms = table(index[Z==1])
	if(any(ms!=1 & ns-ms!=1))
	{
		stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
	}
	if(any(Z!=0 & Z!=1))
	{
		stop("Treatment Vector (Z) Must be Binary")
	}
	if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
	{
		stop("Alternative options are two.sided, greater, or less")
		
	}
	Gammachange = 1
	Gamma.vec = Gamma
	res = sensitivity(index, q, Z, alpha = alpha, alternative = alternative, Gamma.vec=1)
	pval = res$pval
	null.expec = res$null.expec
	Tobs = res$Tobs
	
	
	pvalsens = rep(0, length(Gamma[Gamma>1]))
	rejectvec = rep(0, length(Gamma[Gamma>1]))
	if(maximizeGamma == T)
	{
	  Gammachange = 1
	  if(pval < alpha)
	  {
	    Gammachange = uniroot(sensitivity2, c(1,1.05), index=index, q=q, Z = Z, alpha = alpha, alternative = alternative, continuousRelax = T, extendInt = "downX")$root
	    if(continuousRelax == F)
	    {
	      Gammachange = uniroot(sensitivity2, c(Gammachange,Gammachange+.05), index=index, q=q, Z=Z, alpha = alpha, alternative = alternative, continuousRelax = F, extendInt = "downX")$root
	    }
	  }
    
    
	}
	
	Gamma = Gamma.vec	
	Gamma = Gamma[Gamma>1]
  u = matrix(0,sum(ns), length(Gamma))
	
	if(length(Gamma) > 0)
	{
		for(i in 1:length(Gamma))
		{
		  res = sensitivity(index, q, Z, alpha, alternative, Gamma.vec = Gamma[i], calculate.pval = calculatePvalue, continuous.relax = continuousRelax)
		  u[,i] = res$u
      if(calculatePvalue == T)
		{
		
      pvalsens[i] = res$pval
		}
		else
		{
			rejectvec[i]=res$Reject
		}
		}
	}
	
	Pvalues = cbind(c(1, Gamma.vec[Gamma.vec>1]), c(pval, pvalsens))		
	colnames(Pvalues) = c("Gamma", "P-value")
	rejectvec = c((pval < alpha), rejectvec)
	rejectvec = cbind(c(1, Gamma.vec[Gamma.vec>1]), rejectvec)	
	colnames(rejectvec) = c("Gamma", "Reject")
	Gamma = c(1, Gamma)
if(maximizeGamma == T)
{
	if(calculatePvalue == T)
	{
	return(list(observed = Tobs, nullExpectation = null.expec, alternative = alternative, Pvalues = Pvalues, maxGamma = Gammachange, u = u))
	}
	else
	{
		return(list(observed = Tobs, nullExpectation = null.expec, alternative = alternative, Reject = rejectvec, maxGamma = Gammachange, u = u))
	}
}
else
{
	if(calculatePvalue == T)
	{
	return(list(observed = Tobs, nullExpectation = null.expec, alternative = alternative, Pvalues = Pvalues, maxGamma = NULL, u = u))
	}
	else
	{
		return(list(observed = Tobs, nullExpectation = null.expec, alternative = alternative, Reject = rejectvec, maxGamma = Gammachange, u = u))
	}
}

}



sensitivity = function(index, q, Z, alpha = .05, alternative = "two.sided", Gamma.vec=1, calculate.pval = T, continuous.relax = F)
{
PVAL = calculate.pval
sdq = sd(q)
q = q/sd(q)
Z = 1*Z
ns = table(index)
ms = table(index[Z==1])
ns.types = (ns-1)
N.total = length(q)
pval = 0
nostratum = length(unique(index))

for(i in 1:nostratum)
{
	ind = which(index==i)
	if(ms[i] > 1)
	{
		qsum = sum(q[ind])
		Z[ind] = 1-Z[ind]
		q[ind] = qsum - q[ind]
	}
}



treatment = Z

N.vars = sum((ns-1))
index = index
null.expec = sum(q/ns[index])
Tobs = sum(treatment*q)
max.e = (sum(treatment*q) > null.expec)
index.symm = rep(1:nostratum,ns.types)

PM = rep(0, N.vars)
PV = rep(0, N.vars)
row.ind = rep(0, 2*N.vars + 1)
col.ind = row.ind
values = row.ind
b = rep(0, nostratum+1)
for(kk in 1:nostratum)
{
  row.ind[which(index.symm==kk)] = rep(kk, (ns.types[kk]))
  col.ind[which(index.symm==kk)] = which(index.symm==kk)  
  values[which(index.symm==kk)] = rep(1, ns.types[kk])
  b[kk] = 1
}
row.ind[(N.vars+1):(2*N.vars+1)] = rep(nostratum + 1, N.vars+1)
col.ind[(N.vars+1):(2*N.vars+1)] = 1:(N.vars+1)

opt.expec = rep(0, length(Gamma.vec))
opt.var = opt.expec
zscore = opt.expec
pvalvec = zscore
Rejectvec = zscore
kappavec = opt.expec
U = matrix(0, sum(ns), length(Gamma.vec))
for(ee in 1:length(Gamma.vec))
{
  Gamma.sens = Gamma.vec[ee]
  
  for(kk in 1:nostratum)
  {
    ind = which(index==kk)
    i=kk
	Q = q[ind]
      
      
      
      qi = Q*max.e - Q*(!max.e)
      ord = order(qi)
      qi.sort = sort(qi)
      
      
      mu = rep(0, length(ind)-1)
      sigma2 = rep(0, length(ind)-1)
      
      
      for(j in 1:(length(ind)-1))
      {
        mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
        sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
      }
      mu[abs(mu) < 1e-8] = 0
      sigma2[sigma2 < 1e-8] = 0
      PM[index.symm == kk] = mu*(max.e) - mu*(!max.e)
      PV[index.symm == kk] = (sigma2)

      
    }
  
  
values[(N.vars+1):(2*N.vars+1)] = c(-PM, sign(Tobs))
  b[nostratum+1] = 0
alpha.opt = alpha
  if(alternative != "two.sided")
  {
  	alpha.opt = 2*alpha
  }

	
  const.dir = c(rep("=", nostratum+1))
  model = list()
if(Gamma.sens==1)
{
  V.test = sum(tapply(PV, index.symm, mean))  
  
   tstat = ((Tobs- null.expec)/sqrt(V.test))
   zed = tstat
  kappavec[ee] = ((Tobs- null.expec)^2 - qchisq(1-alpha.opt,1)*V.test)
    pval = 0
    tstat = zed
  if(alternative == "two.sided")
  	{
  		pval = 2*pnorm(-abs(tstat))
  	}
  	if(alternative == "greater")
  	{
  		pval = 1 - pnorm((tstat))
  	}
  	if(alternative == "less")
  	{
  		pval = pnorm((tstat))
  	}
  	Reject = (pval < alpha)

  }
  if(Gamma.sens != 1)
  {
  	diff = 10
  	kappa = qchisq(1-alpha.opt, 1)
  	count=0
  	while(diff > 1e-8)
  	{
  	Plin = -2*Tobs*PM - kappa*PV 
  	rowind.q =  1:(N.vars+1)
  colind.q = 1:(N.vars+1)
  values.q = c(rep(0, N.vars),1)
  Q = sparseMatrix(rowind.q, colind.q, x=values.q)
  	model$A = sparseMatrix(row.ind, col.ind, x=values)
  	model$obj = c(Plin,0)
  	model$Q = Q
  	model$sense = const.dir
  	model$rhs = b
  	 model$vtype = c(rep("I", N.vars), "C")
  if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
  model$lb = c(rep(0, N.vars), -Inf)


  model$modelsense = "min"
  
  
  solm = gurobi(model, params = list(OutputFlag = 0))
  x = solm$x[1:N.vars]
  kappa.new = (Tobs - sum(PM*x))^2/sum(PV*x)
  kappavec[ee] = (Tobs - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
  diff = abs(kappa.new - kappa)
  pval = 0
  if(PVAL == F)
  {
  	diff = 0
  }
  kappa = kappa.new

  }
  if(continuous.relax == F)
  {
  for(kk in 1:nostratum)
  {
    ind = which(index==kk)
    i=kk
    Q = q[ind]
    xind = x[index.symm == kk]
    uvec = c(rep(0, which(xind ==1)), rep(1, length(Q) - which(xind==1)))
    qi = Q*max.e - Q*(!max.e)
    ord = order(qi)
    U[ind[ord],ee] = uvec 
  }  
    
  }
  zed = sqrt((Tobs - sum(PM*x))^2/sum(PV*x))
  if(alternative == "less")
  {
  	zed = -zed
  }
  zscore[ee] = zed
  tstat = zed
  
   
    if(alternative == "two.sided")
  	{
  		pval = 2*pnorm(-abs(tstat))
  	}
  	if(alternative == "greater")
  	{
  		pval = 1 - pnorm((tstat))
  	}
  	if(alternative == "less")
  	{
  		pval = pnorm((tstat))
  	}
  	Reject = (pval < alpha)
  
  
  
  if(sign(Tobs- sum(PM*x))!=sign(Tobs - null.expec))
  {
  	Reject = F
 	pval = 0.5
    
 	if(alternative == "two.sided")
 	{
 		pval = 1
 	}
  }
  
  	if(alternative == "greater" & sum(PM*x) < null.expec)
  	{
  		pval = .5
  	}
  	if(alternative == "less" & sum(PM*x) > null.expec)
  	{
  		pval = .5
  	}
  	
}
pvalvec[ee] = pval
Rejectvec[ee] = Reject   
}    

if(PVAL == F)
{
return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, Tobs = Tobs*sdq, null.expec = null.expec*sdq, pval = NULL, kappa = kappavec, u = U))
}
if(PVAL == T)
{
	return(list(Gamma.vec = Gamma.vec, pval = pvalvec, Tobs = Tobs*sdq, null.expec = null.expec*sdq, u = U))
}

}

sort.new = function(x)
{ 
  temp = sort(unique(x))
  new = 1:length(temp)
  ret = rep(0,length(x))
  for(i in new)
  {
    ret[x == temp[i]] = i
  }
  ret
}

sensitivity2 = function(Gamma.vec, index, q, Z, alpha = 0.05, alternative = "two.sided", continuousRelax = T)
{
  sensitivity(index, q, Z, alpha = alpha, alternative = alternative, Gamma.vec=Gamma.vec, calculate.pval = F, continuous.relax = continuousRelax)$kappa
}   

