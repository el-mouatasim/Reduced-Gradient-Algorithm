%The program of reduced gradient method for the problem:
% min f(X)
% s.t. AX'=b
%       X >= 0
% A. El Mouatasim 2016
function [FX,X,i]=SPRGB_v2a(F,gradF,A,b,X0,maxiter,ksto)
[m,n]=size(A);
XX(:,1)=X0';
%disp(F(X0),'F=')
i=1;
 [B,N,IB,IN]=basic(X0,A,2); %artificiel variables
while i<maxiter   
 lm=lmbda(i);  
  %disp('>>>>>>>>>>--Iteration--<<<<<<<<<<<');
  %disp(i);
   X=XX(:,i)';
    FX=F(X);
   [phi,hors]=cont(A,b,X);
   if hors==1
       disp('hours contraints')
       break
   end
   sb=find(X(IB)==0);
   ns=size(sb,1);
   if ns>0
        [B,N,IB,IN]=pivot(A,B,N,IB,IN,X,sb);
    end
    XB=X(IB); 
    XN=X(IN);  
   g=gradF(X)'; 
   gB=g(IB); 
   gN=g(IN); 
   [d,dB,dN]=direction(B,N,XN,gB,gN,n,m,IB,IN);
   %ndN=norm(dN)
   % if ndN< 10e-10 
    %    disp('X is a Kuhn-Tucker point')
    %end
   [PX]=perturbation(A,b,X,F,ksto,d,lm,IN,IB,B,N,gradF);
        i=i+1;
        XX(:,i)=PX;
   % if  XX(:,i)==XX(:,i-1)
   %         break
   % end
       
end
function [B,N,IB,IN]=basic(X,A,p)
[m,n]=size(A);
base=[];
if p==1 
    [Y,in]=sort(X);
else
    in=1:n;
end
in=1:n;
k1=0;
k2=0;
cb=[];
cn=[];
for j=in
    famille=[base,A(:,j)];
    if rank(famille)>rank(base) & X(j)>0 % Y
        base=famille;
        k1=k1+1;
        cb(k1)=j;
    else
        k2=k2+1;
       cn(k2)=j; 
    end
end
   IB=cb(1:m);
   B=A(:,IB);
   IN=[cb(m+1:k1),cn]; 
   N=A(:,IN);
function [B,N,IB,IN]=pivot(A,B,N,IB,IN,X,sb)
    [m,m]=size(B);
    base=B;
    base(:,sb)=[];
    xn=X(IN);
    [y,in]=sort(xn);
k1=0;
cb=[];
for j=in
    famille=[base,A(:,IN(j))];
    if rank(famille)>rank(base) && xn(j)>0 && rank(famille)<=m
        base=famille;
        k1=k1+1;
        cb(k1)=j;
    end
end
           r = IN(cb);
           s=IB(sb);
           IB(sb)=r;
           IN(cb)=s;
            Fix=B;
            B(:,sb)=N(:,cb);
            N(:,cb)=Fix(:,sb);       
function [d,dB,dN]=direction(B,N,XN,gB,gN,n,m,IB,IN)
    IBN=-inv(B)*N;
    Nr=gB*IBN+gN;
%   for j=1:n-m % B4 ,hs48, hs62
%        if 0 <= Nr(j) && XN(j)==0 
%            dN(j)=0;
%        else 
%            dN(j)=-Nr(j);
%        end
%    end
  for j=1:n-m   %ref book nonlinear...
        if Nr(j) <=  0 
            dN(j)=Nr(j);
        else 
            dN(j)=-XN(j)*Nr(j);
        end
  end
    dB=IBN*dN';
    d(IB)=dB;
    d(IN)=dN;
function [tt]=step(F,d,X,gradF)
    ip = find(d<0);
    sp = length(ip);
    if sp == 0
        lmax = Inf;
    else
    S=-X(ip)./d(ip);
    lmax=min(S);
    end
   if (gradF(X)*d)*(gradF(X+lmax*d)*d)<0 
        lopt=bisection(lmax,F,X,d);
   else 
        lopt=lmax;
   end
    lk = lopt;
    tt=X+lk*d;
function l=bisection(b,F,X,d)
    prec=0.0001; 
    if b>10e+3
        prec = 1;
    end
    xa=0; xb=b;xm=0.5*(b); l=xm;
    bd=xb*d; md=xm*d; xxb=X+bd; xxm=X+md;
    fia=F(X); fib=F(xxb); fim=F(xxm);
    while abs(xb-xa)>prec 
        xl=0.5*(xa+xm); xr=0.5*(xb+xm);
        fil=F(X+xl*d); fir=F(X+xr*d);
        f(1)=fia; f(2)=fib; f(3)=fim; f(4)=fil; f(5)=fir;
        fimin=min(f);
        if fimin==fia || fimin==fil 
            xb=xm; xm=xl;
            fib=fim; fim=fil;
            l=xl;
        elseif fimin==fim 
            xa=xl; xb=xr;
            fia=fil; fib=fir;
            l=xm;
        elseif fimin==fir || fimin==fib
            xa=xm; xm=xr;
            fia=fim; fim=fir;
            l=xr;
        end
    end
    % stochastic perturbation
function lm=lmbda(k)
    lm = 1/sqrt(log(k+2));
function z = gauss(n_pert,n)
randn('state', sum(100*clock));
z = randn(n, n_pert);
%rand('state', sum(100*clock));
%z = random('Normal',0,1,n,n_pert);

           
function [phi,hors]=cont(A,b,X)
    m=length(b);
    hors=0;
    phi=A*X'-b;
    xp=find(X<-0.00001);
    lxp=length(xp);
    for i=1:m
        if abs(phi(i)) > 10e-7  || lxp > 0
        hors=1;
        end    
    end
    function [PX]=perturbation(A,b,X,F,ksto,d,lm,IN,IB,B,N,gradF)
        [m,n]=size(A);
        tt(ksto,:)=X;
        w(ksto)=F(X);
        Y=step(F,d,X,gradF);
        tt(1,:)=Y;
        w(1)=F(Y);
        [z] = gauss(ksto,n-m);
        for k =2: ksto-1
            zz = z(:,k);
            dlPZ(IN) =lm*zz;
            PN=dlPZ(IN);
            PN=PN';
            BN=-inv(B)*N;
            dlPZ(IB)=BN*PN;
             ttt=step(F,dlPZ,Y,gradF);
             [phi,hors]=cont(A,b,ttt);
             if hors == 1 
                w(k)=w(k-1);
                tt(k,:)=tt(k-1,:);
             end
            if hors == 0 
               tt(k,:)=ttt;
               w(k)=F(ttt);    
            end
        end
        [vm,imin]=min(w);
        PX=tt(imin,:)';
   


