% Riccati Equation Kalman Filter
% Per Step Algorithm

function[ppss,iterations]=repsa(f,h,q,r,epsilon)
    k=0;
    p=q;
    k=k+1;
    p1=q+f*p*f'-f*p*h'*inv(h*p*h'+r)*h*p*f';
    d=p1-p;
        while norm(d)>epsilon
            k=k+1;
            p1=q+f*p*f'-f*p*h'*inv(h*p*h'+r)*h*p*f';
            d=p1-p; p=p1;
        end
    ppss=p;
    iterations=k;
