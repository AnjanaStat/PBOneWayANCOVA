#' Finds test statistic value and critical value of asymptotic Min-T test for testing equality  of treatment effects against ordered alternatives in one-way ANCOVA model with single covariate
#'
#' More detailed description
#'
#' @param Y a real dataset
#' @param X a real dataset
#' @param k a positive integer
#' @param alpha real number between 0 and 1 called significance level
#'
#' @return numeric values
#'
#' @examples
#' k=4;N=c(20,30,20,50);S=c(1,1,1,2);alpha=0.05
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=rnorm(N[i],2,sqrt(S[i]))
#' }
#' Y=g
#' G=NULL
#' for(i in 1:k)
#' {
#' G[[i]]=rnorm(N[i],3,sqrt(S[i]))
#' }
#' X=G
#' AMinT1(Y,X,4,0.05)
#' @export

AMinT1<-function(Y,X,k,alpha)
{
  fun1<-function(data1,data2,k)
  {
    Y=data1
    X=data2
    N=unlist(rbind(lapply(Y,length)))
    yM=unlist(rbind(lapply(Y,mean)))
    xM=unlist(rbind(lapply(X,mean)))
    t2=NULL
    for(i in 1:k)
    {
      tm1=sum((Y[[i]]-yM[i])*(X[[i]]-xM[i]));tm2=sum((X[[i]]-xM[i])^2)
      i2=(sum((Y[[i]]-yM[i])^2)-tm1^2/tm2)/(N[i]-2)
      t2[i]=i2
    }
    Sxx=NULL
    for(i in 1:k)
    {
      tm3=sum((X[[i]]-xM[i])^2)
      Sxx[i]=tm3
    }
    s=sum(Sxx)
    t3=NULL
    for(i in 1:k)
    {
      tm4=t2[i]*Sxx[i]
      t3[i]=tm4
    }
    psi=sum(t3)/s^2
    t4=NULL
    for(i in 1:k)
    {
      tm5=sum((Y[[i]]-yM[i])*(X[[i]]-xM[i]))
      t4[i]=tm5
    }
    b=sum(t4)/s
    a=NULL
    for(i in 1:k)
    {
      tm6=yM[i]-b*xM[i]
      a[i]=tm6
    }
    T=NULL
    for(i in 1:k-1)
    {
      tm6=(t2[i]/N[i])+(t2[i+1]/N[i+1])+psi*(xM[i+1]-xM[i])^2
      tm7=sqrt(tm6)
      T[i]=(a[i+1]-a[i])/tm7
    }
    value=min(T)
    return(value)
  }
  fun2<-function(S,N,nu,d,k,alpha)
  {
    x=NULL
    for ( i in 1:k-2)
    {
      temp=-sqrt(nu[i+2]*nu[i])*S[i+1]
      x[i]=temp
    }
    m <- diag(d^2)
    m[row(m)-col(m)== 1]<-x
    m[row(m)-col(m)==-1]<-x
    D=diag(d)
    R=solve(D)%*%m%*%solve(D)
    q=mvtnorm::qmvnorm(alpha,tail="upper.tail",mean=0,sigma=R)$quantile
    return(q)
  }
  data1<-lapply(Y, function(col)col[!is.na(col)])
  data2<-lapply(X, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  yM=unlist(rbind(lapply(data1,mean)))
  xM=unlist(rbind(lapply(data2,mean)))
  S=NULL
  for(i in 1:k)
  {
    tm9=sum((data1[[i]]-yM[i])*(data2[[i]]-xM[i]));tm10=sum((data2[[i]]-xM[i])^2)
    i5=(sum((data1[[i]]-yM[i])^2)-tm9^2/tm10)/(N[i]-2)
    S[i]=i5
  }
  Ntot=sum(N)
  nu=N/Ntot
  d=NULL
  for(i in 1:k-1)
  {
    temp=sqrt((nu[i]*S[i+1])+(nu[i+1]*S[i]))
    d[i]=temp
  }
  set.seed(466)
  test_statistic<-fun1(data1,data2,k)
  crit_value<-fun2(S,N,nu,d,k,alpha)
  result<-c(test_statistic,crit_value)
  print("test statistic value and critical value")
  print(result)
  r1=result[1];r2=result[2]
  if(r1>r2)
  {
    print("Null hypothesis is rejected")
  }
  else
  {
    print("Null hypothesis is not rejected")
  }
}
