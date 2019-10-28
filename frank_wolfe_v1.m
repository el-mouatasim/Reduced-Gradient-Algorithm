function [X,fval,i]=frank_wolfe_v1(f,gradf,X0,e,A,b,Aeq,beq,lb,ub)
%
% A.E.
%
% function function [X,fval,i]=frank_wolfe_v1(f,gradf,X0,e,A,b,Aeq,beq,lb,ub)
%
% order = 3;
%
% Input:    f       cost function
%           X0      starting feasible point
%           e       stopping criteria
%           A       defined in linprog
%           b      
%           Aeq  
%           beq  
%           lb   
%           ub      defined in linprog
%
% Output:   X       optimal point
%           fval    cost at the optimal point
%           i       iterations

X=X0;
%g_f=gradf(X)';
%[p,feval]=linprog(g_f,A,b,Aeq,beq,lb,ub);
%p=p';
i=0;
while(abs(f(X)-e)>10e-6 && i<100)%stopping criteria
g_f=gradf(X)';
[p,feval]=linprog(g_f,A,b,Aeq,beq,lb,ub);
p=p';
d=p-X;
f1=@(lambda) f0(X,lambda,p,f);
[lambda,fval]=fminbnd(f1,0,1);
X = X +lambda*d;
%X=step(f,d,X,gradf);
i=i+1;
end
fval=f(X);
function [tt]=step(F,d,X,gradF)

%    if (gradF(X)*d)*(gradF(X+d)*d)<0 
%         lopt=bisection(1,F,X,d);
%    else 
%         lopt=lmax;
%     end
    lopt=bisection(1,F,X,d);
    lk = min(1,lopt);
    %lk = lopt;
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
function fk=f0(X,lambda,p,f)
    fk=f(X+lambda*(p-X));
    