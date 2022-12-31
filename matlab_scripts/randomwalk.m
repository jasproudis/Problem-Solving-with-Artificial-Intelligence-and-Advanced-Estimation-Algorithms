function randomwalk(kmax)
% random walk
% rng default;
% model parameters
F=1;
H=1;
Q=1;
R=1;
% initial conditions
x0=0;
P0=0;
% time invariant Kalman filter
% x-axis movement
X=x0;
x=x0;
xe=x0; XE=xe;
pe=P0;
for i=1:kmax
    x=F*x+randn; X=[X x];
    z=H*x+randn;
    xp=F*xe;
    pp=Q+F*pe*F';
    g=pp*H'*inv(H*pp*H'+R);
    xe=(1-g*H)*xp+g*z; XE=[XE xe];
    pe=(1-g*H)*pp;
end;
Xreal=X;
Xest=XE;
% y-axis movement
X=x0;
x=x0;
xe=x0; XE=xe;
pe=P0;
for i=1:kmax
    x=F*x+randn; X=[X x];
    z=H*x+randn;
    xp=F*xe;
    pp=Q+F*pe*F';
    g=pp*H'*inv(H*pp*H'+R);
    xe=(1-g*H)*xp+g*z; XE=[XE xe];
    pe=(1-g*H)*pp;
end;
Yreal=X;
Yest=XE;
dreal=sqrt(Xreal.^2+Yreal.^2);
dest=sqrt(Xest.^2+Yest.^2);
% plots
figure(1);
tk=[1:kmax+1];
plot(tk,dreal,'b',tk,dest,'r--');
xlabel('time');
ylabel('distance');
legend('real','estimated');
figure(2);
plot(Xreal,Yreal,'b',Xest,Yest,'r--');
xlabel('x-axis');
ylabel('y-axis');
legend('real','estimated');
title('position');