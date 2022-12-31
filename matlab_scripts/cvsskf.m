function cvsskf(realvx,realvy)
% constant velocity
% Steady State Kalman Filter
% 2 dimensions
% molel
n=4;
dt=0.1;
F=[1 0 dt 0;0 1 0 dt;0 0 1 0;0 0 0 1];
H=eye(n);
sigma1=0.3;
Q=sigma1^2*[0.25*(dt)^4 0 0 0;0 0.25*(dt)^4 0 0;0 0 (dt)^2 0 ;0 0 0 (dt)^2];
sigma2=3; sigma3=0.03;
R=[sigma2^2 0 0 0;0 sigma2^2 0 0;0 0 sigma3^2 0;0 0 0 sigma3^2];
% initial conditions
x0=[0 0 0 0]';
sigma4=4; sigma5=0.4;
pe0=[sigma4^2 0 0 0;0 sigma4^2 0 0;0 0 sigma5^2 0;0 0 0 sigma5^2];
% observability
om=obsv(H,F);
display(om);
rom=rank(om);
display(rom);
% Riccati Equation
PPss=dare(F',H',Q,R);
% Steady State Kalman Filter parameters
Gss=PPss*H'*inv(H*PPss*H'+R);
Ass=(eye(n)-Gss*H)*F;
display(Ass);
display(Gss);
% number of iterations
kmax=1000/dt;
% measurements
sigmav=0.03; sigmap=3;
% rng default;
for k=0:kmax
    vx(k+1)=realvx; zvx(k+1)=vx(k+1)+sigmav*randn;
    vy(k+1)=realvy; zvy(k+1)=vy(k+1)+sigmav*randn;
    px(k+1)=realvx*k*dt; zpx(k+1)=px(k+1)+sigmap*randn;
    py(k+1)=realvy*k*dt; zpy(k+1)=py(k+1)+sigmap*randn;
end;
X=[px; py; vx; vy];
Z=[zpx; zpy; zvx; zvy];
% Steady State Kalman Filter
xe=(eye(n)-Gss*H)*x0+Gss*Z(:,1:1);
XESSKF=xe;
for k=0:kmax-1
    xe=Ass*xe+Gss*Z(:,k+1:k+1);
    XESSKF=[XESSKF xe];
end;
% plots
figure(1);
timet=[0:kmax];
subplot(2,2,1); plot(timet,X(1:1,:),'b',timet,XESSKF(1:1,:),'r--'); xlabel('time'); ylabel('px');
subplot(2,2,2); plot(timet,X(2:2,:),'b',timet,XESSKF(2:2,:),'r--'); xlabel('time'); ylabel('py');
subplot(2,2,3); plot(timet,X(3:3,:),'b',timet,XESSKF(3:3,:),'r--'); xlabel('time'); ylabel('vx');
subplot(2,2,4); plot(timet,X(4:4,:),'b',timet,XESSKF(4:4,:),'r--'); xlabel('time'); ylabel('vy');
% percent Mean Absolute Error
X=X(:,2:kmax+1);
XESSKF=XESSKF(:,2:kmax+1);
pcMAE1=sum(abs(X(1:1,:)-XESSKF(1:1,:))./X(1:1,:))*100/kmax;
pcMAE2=sum(abs(X(2:2,:)-XESSKF(2:2,:))./X(2:2,:))*100/kmax;
pcMAE3=sum(abs(X(3:3,:)-XESSKF(3:3,:))./X(3:3,:))*100/kmax;
pcMAE4=sum(abs(X(4:4,:)-XESSKF(4:4,:))./X(4:4,:))*100/kmax;
pcMAE=[pcMAE1 pcMAE2 pcMAE3 pcMAE4];
display(pcMAE)