function[ppss,iterations]=reda(f,h,q,r,epsilon) 
% Riccati Equation Kalman Filter 
% Doubling Algorithm 
n=size(q,1); 
mon=eye(n); 
k=1; 
a=f'; 
b=h'*inv(r)*h; 
c=q; 
    k=k+1; 
        a1=a*inv(mon+b*c)*a; 
        b1=b+a*inv(mon+b*c)*b*a'; 
        c1=c+a'*c*inv(mon+b*c)*a; 
        d=c1-c; 
        a=a1; 
        b=b1; 
        c=c1; 
    while norm(d)>epsilon 
        k=k+1; 
            a1=a*inv(mon+b*c)*a; 
            b1=b+a*inv(mon+b*c)*b*a'; 
            c1=c+a'*c*inv(mon+b*c)*a; 
            d=c1-c; 
            a=a1; 
            b=b1; 
            c=c1; 
    end 
ppss=c; 
iterations=k; 
