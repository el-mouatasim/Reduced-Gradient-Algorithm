%The program of reduced gradient method for the problem:
% min f(X)
% s.t. AX'=b
%       X >= 0
% A. El Mouatasim 2016
function [FX,X,i,ndN]=RGB(F,gradF,A,b,X0)
[m,n]=size(A);
XX(:,1)=X0';
%disp(F(X0),'F=')
maxiter=150; 
i=1;
 [B,N,IB,IN]=basic(X0,A,1); 
while i<maxiter
    X=XX(:,i)';
    FX=F(X);
 % disp('>>>>>>>>>>--Iteration <  '+string(i)+'  >--<<<<<<<<<<<'); 
   sb=find(X(IB)==0);
   ns=size(sb,1);
   if ns>0
        [B,N,IB,IN]=change(A,B,N,IB,IN,X,sb);
    end
    XB=X(IB); 
    XN=X(IN);  
   g=gradF(X)'; 
   gB=g(IB); 
   gN=g(IN); 
   [d,dB,dN]=direction(B,N,XN,gB,gN,n,m,IB,IN);
   ndN=norm(dN);
    if ndN< 10e-5 
        disp('X is a Kuhn-Tucker point')
        break
    end
    [tt]=step(F,gradF,d,X);
    i=i+1;
    XX(:,i)=tt';
end
function [B,N,IB,IN]=basic(X,A,p)
[m,n]=size(A);
base=[];
in=1:n;
if p==1 
    [Y,in]=sort(X);
end
k1=0;
k2=0;
cb=[];
cn=[];
for j=in
    famille=[base,A(:,j)];
    if rank(famille)>rank(base) && X(j)>0
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
function [B,N,IB,IN]=change(A,B,N,IB,IN,X,sb)
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
   % lmda=[zeros(m,1);Nr'];
   for j=1:n-m
        if 0 <= Nr(j) && XN(j)==0 
            dN(j)=0;
        else 
            dN(j)=-Nr(j);
        end
  end
    dB=IBN*dN';
    d(IB)=dB;
    d(IN)=dN;
function [tt]=step(F,gradF,d,X)
    ip = find(d<0);
    S=-X(ip)./d(ip);
    lmax=min(S);
    if (gradF(X)*d)*(gradF(X+lmax*d)*d)<0 
        lopt=bisection(0,lmax,F,X,d);
    else 
        lopt=lmax;
    end
    lk=min(lopt,lmax);
   tt=X+lk*d;
function l=bisection(a,b,F,X,d)
    prec=0.0001; 
    xa=a; xb=b; xm=0.5*(a+b); l=xm;
    fia=F(X+xa*d'); fib=F(X+xb*d'); fim=F(X+xm*d');
    while abs(xb-xa)>prec 
        xl=0.5*(xa+xm); xr=0.5*(xb+xm);
        fil=F(X+xl*d'); fir=F(X+xr*d');
        fimin=min(fia,fib,fim,fil,fir);
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