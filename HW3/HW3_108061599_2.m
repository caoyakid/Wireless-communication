clc;
clear all; 
close all; 
%% when M = 8
M=8;
N=4*M+2;
fm1=0.01;
fm2=0.1;
fm3=0.5
T=1;
B=pi/M;
fn1=zeros(1,N);
fn2=zeros(1,N);
fn3=zeros(1,N);
for a=1:1:N
    fn1(a)=fm1*cos(2*pi*a/N);
    fn2(a)=fm2*cos(2*pi*a/N);
    fn3(a)=fm3*cos(2*pi*a/N);
end
for t=1:1:1001
    for b=1:1:M
    x1(b)=cos(B*b)*cos(2*pi*fn1(b)*t);
    y1(b)=sin(B*b)*cos(2*pi*fn1(b)*t);
    x2(b)=cos(B*b)*cos(2*pi*fn2(b)*t);
    y2(b)=sin(B*b)*cos(2*pi*fn2(b)*t);
    x3(b)=cos(B*b)*cos(2*pi*fn3(b)*t);
    y3(b)=sin(B*b)*cos(2*pi*fn3(b)*t);
    end
% if α=0
gI1(t)=sqrt(2)*(2*sum(x1)+sqrt(2)*cos(2*pi*fm1*t)); 
gQ2(t)=sqrt(2)*(2*sum(y1));
g1(t)=gI1(t)+j*gQ2(t);
g1dB(t)=10*log10(abs(g1(t)));

gI2(t)=sqrt(2)*(2*sum(x2)+sqrt(2)*cos(2*pi*fm2*t)); 
gQ2(t)=sqrt(2)*(2*sum(y2));
g2(t)=gI2(t)+j*gQ2(t);
g2dB(t)=10*log10(abs(g2(t)));

gI3(t)=sqrt(2)*(2*sum(x3)+sqrt(2)*cos(2*pi*fm3*t)); 
gQ3(t)=sqrt(2)*(2*sum(y3));
g3(t)=gI3(t)+j*gQ3(t);
g3dB(t)=10*log10(abs(g3(t)));
end
t=1:1:300;
subplot(2,6,1);
plot(t,g1dB(t));
title('M=8  fm*T=0.01');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,6,2);
plot(t,g2dB(t));
title('M=8  fm*T=0.1');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,6,3);
plot(t,g3dB(t));
title('M=8  fm*T=0.5');
xlabel('t/T');
ylabel('Envelope Level(dB)');

[ACF,lags,bounds]=autocorr(g1,1000);
t=0:0.01:10;
subplot(2,6,7);
plot(t,ACF);
title('M=8  fm*T=0.01');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');
[ACF,lags,bounds]=autocorr(g2,100);
t=0:0.1:10;
subplot(2,6,8);
plot(t,ACF);
title('M=8  fm*T=0.1');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');

[ACF,lags,bounds]=autocorr(g3,100);
t=0:0.1:10;
subplot(2,6,9);
plot(t,ACF);
title('M=8  fm*T=0.5');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');

%% when M = 16
M=16;
N=4*M+2;
fm1=0.01;
fm2=0.1;
fm3=0.5
T=1;
B=pi/M;
fn1=zeros(1,N);
fn2=zeros(1,N);
fn3=zeros(1,N);
for a=1:1:N
    fn1(a)=fm1*cos(2*pi*a/N);
    fn2(a)=fm2*cos(2*pi*a/N);
    fn3(a)=fm3*cos(2*pi*a/N);
end
for t=1:1:1001
    for b=1:1:M
    x1(b)=cos(B*b)*cos(2*pi*fn1(b)*t);
    y1(b)=sin(B*b)*cos(2*pi*fn1(b)*t);
    x2(b)=cos(B*b)*cos(2*pi*fn2(b)*t);
    y2(b)=sin(B*b)*cos(2*pi*fn2(b)*t);
    x3(b)=cos(B*b)*cos(2*pi*fn3(b)*t);
    y3(b)=sin(B*b)*cos(2*pi*fn3(b)*t);
    end
% if α=0
gI1(t)=sqrt(2)*(2*sum(x1)+sqrt(2)*cos(2*pi*fm1*t)); 
gQ2(t)=sqrt(2)*(2*sum(y1));
g1(t)=gI1(t)+j*gQ2(t);
g1dB(t)=10*log10(abs(g1(t)));

gI2(t)=sqrt(2)*(2*sum(x2)+sqrt(2)*cos(2*pi*fm2*t)); 
gQ2(t)=sqrt(2)*(2*sum(y2));
g2(t)=gI2(t)+j*gQ2(t);
g2dB(t)=10*log10(abs(g2(t)));

gI3(t)=sqrt(2)*(2*sum(x3)+sqrt(2)*cos(2*pi*fm3*t)); 
gQ3(t)=sqrt(2)*(2*sum(y3));
g3(t)=gI3(t)+j*gQ3(t);
g3dB(t)=10*log10(abs(g3(t)));
end
t=1:1:300;
subplot(2,6,4);
plot(t,g1dB(t));
title('M=16  fm*T=0.01');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,6,5);
plot(t,g2dB(t));
title('M=16  fm*T=0.1');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,6,6);
plot(t,g3dB(t));
title('M=16  fm*T=0.5');
xlabel('t/T');
ylabel('Envelope Level(dB)');

[ACF,lags,bounds]=autocorr(g1,1000);
t=0:0.01:10;
subplot(2,6,10);
plot(t,ACF);
title('M=16  fm*T=0.01');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');
[ACF,lags,bounds]=autocorr(g2,100);
t=0:0.1:10;
subplot(2,6,11);
plot(t,ACF);
title('M=16  fm*T=0.1');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');

[ACF,lags,bounds]=autocorr(g3,100);
t=0:0.1:10;
subplot(2,6,12);
plot(t,ACF);
title('M=16  fm*T=0.5');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');