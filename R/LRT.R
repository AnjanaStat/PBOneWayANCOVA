#' It gives the test statistic value and critical value of likelihood ratio test for testing homogeneity of treatment effects against ordered alternatives in one-way ANCOVA model
#'
#' More detailed description
#'
#' @param Y a real data set
#' @param X a real data set
#' @param k a positive integer
#' @param q a positive integer
#' @param alpha real number between 0 and 1 called significance level
#' @param B a positive integer
#'
#' @return numeric vector
#'
#' @examples
#' k=4;q=3;B=5000;alpha=0.05
#' N=c(20,10,10,20,20,10,10,20,20,10,10,20);S=c(1,1,2,1,1,2,3,1,2,4,6,3,2,4)
#' g=NULL
#' for(i in 1:(k*q))
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' X=g
#' G2=NULL
#' N=c(20,10,10,20);a=c(1,2,3,4)
#' for(i in 1:k)
#' {
#'  G2[[i]]=rnorm(N[i],a[i],sqrt(S[i]))
#' }
#' Y=G2
#' LRT(Y,X,k,q,alpha,B)
#' @export
LRT<-function(Y,X,k,q,alpha,B)
{
  fun1<-function(data1,data2,k,q)
  {
    Y=data1
    X=data2
    N=unlist(rbind(lapply(Y,length)))
    yM=unlist(rbind(lapply(Y,mean)))
    xM=unlist(rbind(lapply(X,mean)))
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
    a0=yM-t5
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
    B0=b0
    S1=t7
    ###Algorithm (MLEs under full parameters space)
    repeat
    {
      t9=NULL;W=NULL;t8=NULL
      for(i in 1:k)
      {
        for(h in 1:q)
        {
          t1=B0[h]*xM[(h-1)*k+i]
          t8[h]=t1
        }
        i3=yM[i]-sum(t8)
        t9[i]=i3
        W[i]=N[i]/S1[i]
      }
      An=Iso::pava(t9,W)
      t10=NULL
      t11=NULL;tm3=NULL
      for(h in 1:q)
      {
        for(l in 1:q)
        {
          for(i in 1:k)
          {
            t1=(X[[(h-1)*k+i]]*X[[(l-1)*k+i]])
            tm3[i]=sum(t1)/S1[i]
          }
          i4=sum(tm3)
          t11[l+(h-1)*q]=i4
        }
      }
      Pxx=matrix(t11,nrow=q,ncol=q,byrow=TRUE)
      t12=NULL;tm4=NULL
      for(h in 1:q)
      {
        for(i in 1:k)
        {
          t1=X[[(h-1)*k+i]]*(Y[[i]]-An[i])
          tm4[i]=sum(t1)/S1[i]
        }
        i2=sum(tm4)
        t12[h]=i2
      }
      Qxy=matrix(t12,nrow=q,ncol=1,byrow=TRUE)
      bn=solve(Pxx)%*%Qxy
      t13=NULL
      t14=NULL
      for(i in 1:k)
      {
        t13=NULL
        for(j in 1:N[i])
        {
          t13[j]=0
        }
        for(h in 1:q)
        {
          t1=bn[h]*X[[(h-1)*k+i]]
          t13=t13+t1
        }
        i2=(Y[[i]]-An[i])-t13
        t14[i]=sum(i2^2)/N[i]
      }
      dif1=max(abs(a0-An));dif2=max(abs(B0-bn))
      if(dif1<=0.00001&dif2<=0.00001)
      {
        break
      }
      a0=An;B0=bn;S1=t14
    }

    ## Finding MLEs under null parameter space
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
          t1=(X[[(h-1)*k+i]]-t15[h])*(X[[(l-1)*k+i]]-t15[l])
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
        t1=(X[[(h-1)*k+i]]-t15[h])*(Y[[i]]-totY)
        tm7[i]=sum(t1)
      }
      i2=sum(tm7)
      t18[h]=i2
    }
    sxy=matrix(t18,nrow=q,ncol=1,byrow=TRUE)
    b01=solve(sxx)%*%sxy
    a01=totY-sum(b01*t15)
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
        t1=b01[h]*(X[[(h-1)*k+i]]-xM[(h-1)*k+i])
        t19=t19+t1
      }
      i2=(Y[[i]]-yM[i])-t19
      t20[i]=sum(i2^2)/(N[i]-q-1)
    }
    B01=b01
    S2=t20
    ###Algorithm (MLEs under null parameter space)
    repeat
    {
      U=NULL
      for(i in 1:k)
      {
        U[i]=N[i]/S2[i]
      }
      tm8=NULL;tm9=NULL
      for(i in 1:k)
      {
        for(h in 1:q)
        {
          t1=B01[h]*xM[(h-1)*k+i]
          tm8[h]=t1
        }
        i3=yM[i]-sum(tm8)
        tm9[i]=i3
      }
      an1=sum(U*tm9)/sum(U)
      t21=NULL
      t22=NULL
      for(h in 1:q)
      {
        for(l in 1:q)
        {
          for(i in 1:k)
          {
            t1=X[[(h-1)*k+i]]*X[[(l-1)*k+i]]
            t21[i]=sum(t1)/S2[i]
          }
          i2=sum(t21)
          t22[l+(h-1)*q]=i2
        }
      }
      pxx=matrix(t22,nrow=q,ncol=q,byrow=TRUE)
      t23=NULL;t24=NULL
      for(h in 1:q)
      {
        for(i in 1:k)
        {
          t1=(Y[[i]]-an1)*X[[(h-1)*k+i]]
          t23[i]=sum(t1)/S2[i]
        }
        i2=sum(t23)
        t24[h]=i2
      }
      qxy=matrix(t24,nrow=q,ncol=1,byrow=TRUE)
      bn1=solve(pxx)%*%qxy
      t25=NULL
      t26=NULL
      for(i in 1:k)
      {
        t25=NULL
        for(j in 1:N[i])
        {
          t25[j]=0
        }
        for(h in 1:q)
        {
          t1=bn1[h]*X[[(h-1)*k+i]]
          t25=t25+t1
        }
        i2=(Y[[i]]-an1)-t25
        t26[i]=sum(i2^2)/N[i]
      }
      dif3=abs(a01-an1);dif4=max(abs(B01-bn1))
      if(dif3<=0.00001&dif4<=0.000001)
      {
        break
      }
      a01=an1;B01=bn1;S2=t26
    }
    rto=S1/S2;ratio=NULL
    for(i in 1:k)
    {
      rslt=(rto[i])^(N[i]/2)
      ratio[i]=rslt
    }
##likelihood ratio value
    value=prod(ratio)
    return(value)
  }

  fun2<-function(data1,N,S,b,k,q)
  {
    X=data1
    t25=NULL
    t26=NULL
    for(i in 1:k)
    {
      t25=NULL
      for(j in 1:N[i])
      {
        t25[j]=0
      }
      for(h in 1:q)
      {
        t1=b[h]*X[[(h-1)*k+i]]
        t25=t25+t1
      }
      t26[[i]]=t25
    }
    g=NULL
    for(i in 1:k)
    {
      g[[i]]=rnorm(N[i],t26[[i]],sqrt(S[i]))
    }
    value=fun1(g,data1,k,q)
    return(value)
  }
  fun3<-function(data1,N,S,b,k,q,alpha,B)
  {
    x<-replicate(B,fun2(data1,N,S,b,k,q))
    y<-sort(x,decreasing=FALSE)
    m=(alpha)*B
    c<-y[m]
    return(c)
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
  set.seed(419)
  test_statistic<-fun1(data1,data2,k,q)
  crit_value<-fun3(data2,N,S,b,k,q,alpha,B)
  result<-c(test_statistic,crit_value)
  print("test statistic value and critical value")
  print(result)
  r1=result[1];r2=result[2]
  if(r1<r2)
  {
    print("Null hypothesis is rejected")
  }
  else
  {
    print("Null hypothesis is not rejected")
  }
}
