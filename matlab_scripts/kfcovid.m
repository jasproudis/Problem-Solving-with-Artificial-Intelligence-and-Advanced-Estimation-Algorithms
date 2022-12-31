function kfcovid
% cases prediction of Covid-19
% steady state Kalman filter

% data cases
Z=[19731 24996 18825 17662 18640 16023 10980 18847 23687 19616];

% steady state Kalman Filter parameters
pp=dare(1,1,1,1);
b = pp/(pp+1);
a=1-b;

% steady state Kalman Filter
x0= 20000;
xe=x0;
XE=[];

for k=1:10
    xe=a*xe+b*Z(k);
    XE=[XE xe];
end;

XESSKF=ceil(XE);

% plots
figure(1);
clf;
timek=[1:10];
plot(timek,Z,'b',timek,XESSKF,'r');
axis([0, 10, 0, 30000]);
legend('real','prediction');
xlabel ('1-10 February 2022');
ylabel ('cases');
title('cases position prediction of Covid-19');
