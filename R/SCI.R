#' It gives lower bound of simultaneous confidence intervals for successive differences of treatment effects.
#'
#' More detailed description
#'
#' @param Y a real data set
#' @param X a real data set
#' @param k a positive integer
#' @param q a positive integer
#' @param alpha real number between 0 and 1 called significance level
#'
#' @return numeric vector
#'
#' @examples
#' k=4;q=3
#' N=c(20,10,10,20,20,10,10,20,20,10,10,20);S=c(1,1,2,1,1,2,3,1,2,4,6,3,2,4)
#' g=NULL
#' for(i in 1:(k*q))
#' {
#'   g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' X=g
#' G2=NULL
#' N=c(20,10,10,20);a=c(1,2,3,4)
#' for(i in 1:k)
#' {
#'  G2[[i]]=rnorm(N[i],a[i],sqrt(S[i]))
#' }
#' Y=G2
#' SCI(Y,X,k,q,0.05)
#' @export


SCI<-function(Y,X,k,q,alpha)
{
  fun1<-function(data1,data2,k,q)
  {
    Y=data1
    X=data2
    N=unlist(rbind(lapply(Y,length)))
    yM=unlist(rbind(lapply(Y,mean)))
    xM=unlist(rbind(lapply(X,mean)))
    Xbar=matrix(xM,nrow=k,ncol=q,byrow=FALSE)
    t2=NULL
    tm1=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(i in 1:k)
        {
          t1=(X[[(h-1)*k+i]]-xM[(h-1)*k+i])*(X[[(l-1)*k+i]]-xM[(l-1)*k+i])
          tm1[i]=sum(t1)
        }
        i2=sum(tm1)
        t2[l+(h-1)*q]=i2
      }
    }
    Sxx=matrix(t2,nrow=q,ncol=q,byrow=TRUE)
    t3=NULL
    tm2=NULL
    t1=NULL
    for(h in 1:q)
    {
      for(i in 1:k)
      {
        t1=(X[[(h-1)*k+i]]-xM[(h-1)*k+i])*(Y[[i]]-yM[i])
        tm2[i]=sum(t1)
      }
      i2=sum(tm2)
      t3[h]=i2
    }
    Sxy=matrix(t3,nrow=q,ncol=1,byrow=TRUE)
    b0=solve(Sxx)%*%Sxy
    t4=NULL
    t5=NULL
    for(i in 1:k)
    {
      for(h in 1:q)
      {
        t1=b0[h]*xM[(h-1)*k+i]
        t4[h]=t1
      }
      t5[i]=sum(t4)
    }
    A=yM-t5
    t6=NULL
    t7=NULL
    for(i in 1:k)
    {
      t6=NULL
      for(j in 1:N[i])
      {
        t6[j]=0
      }
      for(h in 1:q)
      {
        t1=b0[h]*(X[[(h-1)*k+i]]-xM[(h-1)*k+i])
        t6=t6+t1
      }
      i2=(Y[[i]]-yM[i])-t6
      t7[i]=sum(i2^2)/(N[i]-q-1)
    }
    S1=t7
    t8=NULL
    tm4=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(i in 1:k)
        {
          t1=(X[[(h-1)*k+i]]-xM[(h-1)*k+i])*(X[[(l-1)*k+i]]-xM[(l-1)*k+i])
          tm4[i]=S1[i]*sum(t1)
        }
        i2=sum(tm4)
        t8[l+(h-1)*q]=i2
      }
    }
    S=matrix(t8,nrow=q,ncol=q,byrow=TRUE)
    t9=S1/N
    S2=diag(t9)+Xbar%*%solve(Sxx)%*%S%*%solve(Sxx)%*%t(Xbar)
    V=NULL
    T=NULL
    for(i in 1:k-1)
    {
      V[i]=sqrt(S2[i,i]+S2[i+1,i+1]-2*S2[i,i+1])
      T[i]=(A[i+1]-A[i])/V[i]
    }
    value=max(T)
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
    q=mvtnorm::qmvnorm(1-alpha,tail="lower.tail",mean=0,sigma=R)$quantile
    return(q)
  }
  data1<-lapply(Y, function(col)col[!is.na(col)])
  data2<-lapply(X, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  yM=unlist(rbind(lapply(data1,mean)))
  xM=unlist(rbind(lapply(data2,mean)))
  t15=NULL;tm5=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=N[i]*xM[(h-1)*k+i]
      tm5[i]=t1
    }
    t15[h]=sum(tm5)/sum(N)
  }
  t16=N*yM
  totY=sum(t16)/sum(N)
  t17=NULL
  tm6=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(i in 1:k)
      {
        t1=(data2[[(h-1)*k+i]]-t15[h])*(data2[[(l-1)*k+i]]-t15[l])
        tm6[i]=sum(t1)
      }
      i2=sum(tm6)
      t17[l+(h-1)*q]=i2
    }
  }
  sxx=matrix(t17,nrow=q,ncol=q,byrow=TRUE)
  t18=NULL;tm7=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=(data2[[(h-1)*k+i]]-t15[h])*(data1[[i]]-totY)
      tm7[i]=sum(t1)
    }
    i2=sum(tm7)
    t18[h]=i2
  }
  sxy=matrix(t18,nrow=q,ncol=1,byrow=TRUE)
  b=solve(sxx)%*%sxy
  t19=NULL
  t20=NULL
  for(i in 1:k)
  {
    t19=NULL
    for(j in 1:N[i])
    {
      t19[j]=0
    }
    for(h in 1:q)
    {
      t1=b[h]*(data2[[(h-1)*k+i]]-xM[(h-1)*k+i])
      t19=t19+t1
    }
    i2=(data1[[i]]-yM[i])-t19
    t20[i]=sum(i2^2)/(N[i]-q-1)
  }
  S=t20
  Ntot=sum(N)
  nu=N/Ntot
  d=NULL
  for(i in 1:k-1)
  {
    temp=sqrt((nu[i]*S[i+1])+(nu[i+1]*S[i]))
    d[i]=temp
  }
  set.seed(465)
  C<-fun2(S,N,nu,d,k,alpha)
  Xbar=matrix(xM,nrow=k,ncol=q,byrow=FALSE)
  t2=NULL
  tm1=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(i in 1:k)
      {
        t1=(data2[[(h-1)*k+i]]-xM[(h-1)*k+i])*(data2[[(l-1)*k+i]]-xM[(l-1)*k+i])
        tm1[i]=sum(t1)
      }
      i2=sum(tm1)
      t2[l+(h-1)*q]=i2
    }
  }
  Sxx=matrix(t2,nrow=q,ncol=q,byrow=TRUE)
  t3=NULL
  tm2=NULL
  t1=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=(data2[[(h-1)*k+i]]-xM[(h-1)*k+i])*(data1[[i]]-yM[i])
      tm2[i]=sum(t1)
    }
    i2=sum(tm2)
    t3[h]=i2
  }
  Sxy=matrix(t3,nrow=q,ncol=1,byrow=TRUE)
  b0=solve(Sxx)%*%Sxy
  t4=NULL
  t5=NULL
  for(i in 1:k)
  {
    for(h in 1:q)
    {
      t1=b0[h]*xM[(h-1)*k+i]
      t4[h]=t1
    }
    t5[i]=sum(t4)
  }
  A=yM-t5
  t6=NULL
  t7=NULL
  for(i in 1:k)
  {
    t6=NULL
    for(j in 1:N[i])
    {
      t6[j]=0
    }
    for(h in 1:q)
    {
      t1=b0[h]*(data2[[(h-1)*k+i]]-xM[(h-1)*k+i])
      t6=t6+t1
    }
    i2=(data1[[i]]-yM[i])-t6
    t7[i]=sum(i2^2)/(N[i]-q-1)
  }
  S1=t7
  t8=NULL
  tm4=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(i in 1:k)
      {
        t1=(data2[[(h-1)*k+i]]-xM[(h-1)*k+i])*(data2[[(l-1)*k+i]]-xM[(l-1)*k+i])
        tm4[i]=S1[i]*sum(t1)
      }
      i2=sum(tm4)
      t8[l+(h-1)*q]=i2
    }
  }
  S=matrix(t8,nrow=q,ncol=q,byrow=TRUE)
  t9=S1/N
  S2=diag(t9)+Xbar%*%solve(Sxx)%*%S%*%solve(Sxx)%*%t(Xbar)
  for(i in 1:k-1)
  {
    print("lower confidence limit for difference between the treatment effects")
    tm5=sqrt(S2[i,i]+S2[i+1,i+1]-2*S2[i,i+1])
    lb=A[i+1]-A[i]-C*tm5
    l=i+1;u=i
    print(l);print(u)
    print(lb)
  }
}

