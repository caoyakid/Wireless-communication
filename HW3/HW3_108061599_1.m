clc;
clear all; 
close all; 
%% variable
fm_1 = 0.01;
fm_2 = 0.1;
fm_3 = 0.5;

T = 1;
omega_p = 1;
% when fm = 0.01
s_1 = 2-cos(pi*fm_1*T/2)-sqrt((2-cos(pi*fm_1*T/2))^2-1);
sigma_1=((1+s_1)*omega_p)/((1-s_1)*2);
gI1(1)=sqrt(0.01)*randn(1);
gQ1(1)=sqrt(0.01)*randn(1);
g1(1)=gI1(1)+j*gQ1(1);
g1dB(1)=10*log10(abs(g1(1)));
% when fm = 0.1
s_2 = 2-cos(pi*fm_2*T/2)-sqrt((2-cos(pi*fm_2*T/2))^2-1);
sigma_2=((1+s_2)*omega_p)/((1-s_2)*2);
gI2(1)=sqrt(0.1)*randn(1);
gQ2(1)=sqrt(0.1)*randn(1);
g2(1)=gI2(1)+j*gQ2(1);
g2dB(1)=10*log10(abs(g2(1)));
% when fm - 0.5
s_3 = 2-cos(pi*fm_3*T/2)-sqrt((2-cos(pi*fm_3*T/2))^2-1);
sigma_3=((1+s_3)*omega_p)/((1-s_3)*2);
gI3(1)=sqrt(0.5)*randn(1);
gQ3(1)=sqrt(0.5)*randn(1);
g3(1)=gI3(1)+j*gQ3(1);
g3dB(1)=10*log10(abs(g3(1)));

for t=1:1:1001    
wa1=sqrt(sigma_1)*randn(1);
wa2=sqrt(sigma_1)*randn(1);
gI1(t+1)=s_1*gI1(t)+(1-s_1)*wa1;
gQ1(t+1)=s_1*gQ1(t)+(1-s_1)*wa2;
g1(t+1)=gI1(t+1)+j*gQ1(t+1);
g1dB(t+1)=10*log10(abs(g1(t+1)));

wb1=sqrt(sigma_2)*randn(1);
wb2=sqrt(sigma_2)*randn(1);
gI2(t+1)=s_2*gI2(t)+(1-s_2)*wb1;
gQ2(t+1)=s_2*gQ2(t)+(1-s_2)*wb2;
g2(t+1)=gI2(t+1)+j*gQ2(t+1);
g2dB(t+1)=10*log10(abs(g2(t+1)));

wc1=sqrt(sigma_3)*randn(1);
wc2=sqrt(sigma_3)*randn(1);
gI3(t+1)=s_3*gI3(t)+(1-s_3)*wc1;
gQ3(t+1)=s_3*gQ3(t)+(1-s_3)*wc2;
g3(t+1)=gI3(t+1)+j*gQ3(t+1);
g3dB(t+1)=10*log10(abs(g3(t+1)));
end

%% Plot
t=1:1:300;
subplot(2,3,1);
plot(t,g1dB(t));
title('fm*T=0.01');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,3,2);
plot(t,g2dB(t));
title('fm*T=0.1');
xlabel('t/T');
ylabel('Envelope Level(dB)');
subplot(2,3,3);
plot(t,g3dB(t));
title('fm*T=0.5');
xlabel('t/T');
ylabel('Envelope Level(dB)');


[ACF,lags,bounds]=autocorr(g1,1000);
t=0:0.01:10;
subplot(2,3,4);
plot(t,ACF);
title('fm*T=0.01');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');

[ACF,lags,bounds]=autocorr(g2,100);
t=0:0.1:10;
subplot(2,3,5);
plot(t,ACF);
title('fm*T=0.1');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');

[ACF,lags,bounds]=autocorr(g3,10);
t=0:1:10;
subplot(2,3,6);
plot(t,ACF);
title('fm*T=0.5');
xlabel('Time difference(fm*\tau)');
ylabel('Autocorrelation');