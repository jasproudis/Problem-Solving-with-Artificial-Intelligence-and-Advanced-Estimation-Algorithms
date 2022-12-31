function casskf(realax,realay,realaz)
% constant acceleration
% Steady State Kalman Filter
% 3 dimensions
% model
n=6;
dt=0.1;
F=[1 0 0 dt 0 0;0 1 0 0 dt 0;0 0 1 0 0 dt;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
B=[0.5*(dt)^2 0 0;0 0.5*(dt)^2 0;0 0 0.5*(dt)^2;dt 0 0;0 dt 0;0 0 dt];
H=eye(n);
sigma1=0.3;
Q=sigma1^2*[0.25*(dt)^4 0 0 0 0 0;0 0.25*(dt)^4 0 0 0 0;0 0 0.25*(dt)^4 0 0 0;0 0 0 (dt)^2 0 0;0 0 0 0 (dt)^2 0;0 0 0 0 0 (dt)^2];
sigma2=3; sigma3=0.03;
R=[sigma2^2 0 0 0 0 0;0 sigma2^2 0 0 0 0;0 0 sigma2^2 0 0 0;0 0 0 sigma3^2 0 0;0 0 0 0 sigma3^2 0;0 0 0 0 0 sigma1^3];
% initial coditions
x0=[0 0 0 0 0 0]';
sigma4=4; sigma5=0.4;
pe0=[sigma4^2 0 0 0 0 0;0 sigma4^2 0 0 0 0;0 0 sigma4^2 0 0 0;0 0 0 sigma5^2 0 0;0 0 0 0 sigma5^2 0;0 0 0 0 0 sigma5^2];
% observability
om=obsv(H,F);
display(om);
rom=rank(om);
display(rom);
% controllability
cm=ctrb(F,B);
display(cm);
rcm=rank(cm);
display(rcm);
% Riccati Equation
PPss=dare(F',H',Q,R);
% Steady State Kalman Filter parameters
Gss=PPss*H'*inv(H*PPss*H'+R);
Ass=(eye(n)-Gss*H)*F;
Css=(eye(n)-Gss*H)*B;
display(Ass);
display(Css);
display(Gss);
% number of iterations
kmax=1000/dt;
% mesurements
sigmaa=0.3; sigmav=0.03; sigmap=3;
% rng default;
for k=0:kmax
    ax(k+1)=realax; zax(k+1)=ax(k+1)+sigmaa*randn;
    ay(k+1)=realay; zay(k+1)=ay(k+1)+sigmaa*randn;
    az(k+1)=realaz; zaz(k+1)=az(k+1)+sigmaa*randn;
    vx(k+1)=realax*k*dt; zvx(k+1)=vx(k+1)+sigmav*randn;
    vy(k+1)=realay*k*dt; zvy(k+1)=vy(k+1)+sigmav*randn;
    vz(k+1)=realaz*k*dt; zvz(k+1)=vz(k+1)+sigmav*randn;
    px(k+1)=realax*0.5*(k*dt)^2;
    zpx(k+1)=px(k+1)+sigmap*randn;
    py(k+1)=realay*0.5*(k*dt)^2;
    zpy(k+1)=py(k+1)+sigmap*randn;
    pz(k+1)=realaz*0.5*(k*dt)^2;
    zpz(k+1)=pz(k+1)+sigmap*randn;
end;
A=[ax; ay; az];
X=[px; py; pz; vx; vy; vz];
Z=[zpx; zpy; zpz; zvx; zvy; zvz];
% Steady State Kalman Filter
xe=(eye(n)-Gss*H)*x0+Gss*Z(:,1:1);
XESSKF=xe;
U=[realax realay realaz]';
for k=0:kmax-1
xe=Ass*xe+Gss*Z(:,k+1:k+1)+Css*U;
XESSKF=[XESSKF xe];
end;
% plots
figure(1);
timet=[0:kmax];
subplot(3,3,1); plot(timet,A(1:1,:)); xlabel('time'); ylabel('ax');
subplot(3,3,2); plot(timet,A(2:2,:)); xlabel('time'); ylabel('ay');
subplot(3,3,3); plot(timet,A(3:3,:)); xlabel('time'); ylabel('az');
subplot(3,3,4); plot(timet,X(4:4,:),'b',timet,XESSKF(4:4,:),'r--'); xlabel('time'); ylabel('vx');
subplot(3,3,5); plot(timet,X(5:5,:),'b',timet,XESSKF(5:5,:),'r--'); xlabel('time'); ylabel('vy');
subplot(3,3,6); plot(timet,X(6:6,:),'b',timet,XESSKF(6:6,:),'r--'); xlabel('time'); ylabel('vz');
subplot(3,3,7); plot(timet,X(1:1,:),'b',timet,XESSKF(1:1,:),'r--'); xlabel('time'); ylabel('px');
subplot(3,3,8); plot(timet,X(2:2,:),'b',timet,XESSKF(2:2,:),'r--'); xlabel('time'); ylabel('py');
subplot(3,3,9); plot(timet,X(3:3,:),'b',timet,XESSKF(3:3,:),'r--'); xlabel('time'); ylabel('pz');
% percent Mean Absolute Error
X=X(:,2:kmax+1);
XESSKF=XESSKF(:,2:kmax+1);
pcMAE1=sum(abs(X(1:1,:)-XESSKF(1:1,:))./X(1:1,:))*100/kmax;
pcMAE2=sum(abs(X(2:2,:)-XESSKF(2:2,:))./X(2:2,:))*100/kmax;
pcMAE3=sum(abs(X(3:3,:)-XESSKF(3:3,:))./X(3:3,:))*100/kmax;
pcMAE4=sum(abs(X(4:4,:)-XESSKF(4:4,:))./X(4:4,:))*100/kmax;
pcMAE5=sum(abs(X(5:5,:)-XESSKF(5:5,:))./X(5:5,:))*100/kmax;
pcMAE6=sum(abs(X(6:6,:)-XESSKF(6:6,:))./X(6:6,:))*100/kmax;
pcMAE=[pcMAE1 pcMAE2 pcMAE3 pcMAE4 pcMAE5 pcMAE6];
display(pcMAE);