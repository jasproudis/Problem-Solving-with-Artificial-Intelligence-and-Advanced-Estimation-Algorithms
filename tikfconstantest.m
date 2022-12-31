function tikfconstantest(c,kmax)
% constant estimation using Time Invariant Kalman Filter

rng default;

% constant
C=c.*ones(1,kmax+1);

% model parameters
F=1;
H=1;
Q=10^(-3);
R=1;

% measurements
Z=[];

for k=0:kmax

    z=c+R*randn; Z=[Z z];

end

% TIKF
% initial conditions
xp0=0;
pp0=1;

xp=xp0;
pp=pp0;
XE=[];
PE=[];

for k=0:kmax

    g=(pp*H)/(H*pp*H+R);
    xe=(1-g*H)*xp+g*Z(k+1);
    pe=(1-g*H)*pp;
    xp=F*xe;
    pp=Q+F*pe*F';
    XE=[XE xe];
    PE=[PE pe];

end
XEKF=XE;
PEKF=PE;

% plots
clf;
figure(1);
timek=[0:kmax];
plot(timek,C,'k',timek,Z,'b',timek,XEKF,'r');
legend('constant','measurements','KF estimation');
xlabel('time k');
ylabel('constant and estimation');

figure(2);
timek=[0:kmax];
plot(timek,PEKF,'r');
xlabel('time k');
ylabel('estimation error variance');
