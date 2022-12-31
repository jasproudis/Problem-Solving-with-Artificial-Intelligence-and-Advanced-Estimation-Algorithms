function[ppss,iterations]=leda(f,q,epsilon)
% Lyapunov Equation Kalman Filter
% Doubling Algorithm
k=1;
a=f';
c=q;
k=k+1;
a1=a*a;
c1=c+a'*c*a;
d=c1-c;
a=a1;
c=c1;
    while norm(d)>epsilon 
        k=k+1;
        a1=a*a;
        c1=c+a'*c*a;
        d=c1-c;
        a=a1;
        c=c1;
    end
ppss=c;
iterations=k;