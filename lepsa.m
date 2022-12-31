function[ppss,iterations]=lepsa(f,q,epsilon)
% Lyapunov Equation Kalman Filter
% Per Step Algorithm
k=0;
p=q;
k=k+1;
p1=q+f*p*f';
d=p1-p;
    while norm(d)>epsilon
        k=k+1;
        p1=q+f*p*f';
        d=p1-p; p=p1;
    end
ppss=p;
iterations=k;