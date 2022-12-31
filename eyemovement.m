function [dreal,destimated]=eyemovement
% eye movement
% steady state Kalman filter
% time
kmax=10000;
kmax=100;
% vevocity
vx=0.4;
vy=0.3;
% real position
POS=[0 0]';
for k=1:kmax
posx=vx*k;
posy=vy*k;
pos=[posx posy]';
POS=[POS pos];
end;
% model
% state = [sx sy vx vy]'
f = [1 0 1 0
0 1 0 1
0 0 1 0
0 0 0 1];
h = [1 0 0 0
0 1 0 0];
q=(10000/6)*[2*eye(2) 3*eye(2)
3*eye(2) 6*eye(2)];
r=102*eye(2);
% steady state Kalman Filter parameters
pp=dare(f',h',q,r);
b=pp*h'*inv(h*pp*h'+r);
a=(eye(4)-b*h)*f;
% steady state Kalman Filter
x0=[0 0 0 0 ]';
p0=[1 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0]/16;
xe=x0;
XE=xe;
for k=1:kmax
% measurements
zx=vx*k+randn*sqrt(102);
zy=vy*k+randn*sqrt(102);
xe=a*xe+b*[zx zy]';
XE=[XE xe];
end;
XESSKF=XE;
% distances
dreal=sqrt(POS(1,length(POS))^2+POS(2,length(POS))^2);
destimated=sqrt(XESSKF(1,length(XESSKF))^2+XESSKF(2,length(XESSKF))^2);
% plots
figure(1);
clf;
timek=[0:kmax];
subplot(2,1,1);
plot(timek,POS(1:1,:),'b',timek,XESSKF(1:1,:),'r');
legend('real','estimation');
xlabel('time k');
ylabel('x position');
subplot(2,1,2);
plot(timek,POS(2:2,:),'b',timek,XESSKF(2:2,:),'r');
legend('real','estimation');
xlabel('time k');
ylabel('y position');
figure(2);
clf;
plot(POS(1:1,:),POS(2:2,:),'b',XESSKF(1:1,:),XESSKF(2:2,:),'r');
legend('real','estimation');
title('position');