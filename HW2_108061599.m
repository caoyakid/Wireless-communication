%% 
clear;
close all;
clc;

%% variable 
c = 3*10^8; %光速
t = 1000000; %用於製造亂數
phase = -pi + 2*pi*rand(1,t); %偏移角度
% a小題部分
v_a = 20000/3600; %速度(轉成公尺/秒
fc_a = 2*10^9; %波長倒數
fm_a = v_a * fc_a/c;  
fD_a = fm_a * cos(phase); %頻偏
sequence_a = round(cos(phase),3); %化簡到小數前3位
prob_a = zeros(1,2001); %創立一個0序列

% b小題部分
v_b = 90000/3600;
fc_b = 26*10^9;
fm_b = v_b * fc_b/c;
fD_b = fm_b * cos(phase);
sequence_b = round(cos(phase),3);
prob_b = zeros(1,2001);

% c小題部分
v_c = 20000/3600 + (70000/3600)*rand(1,t);
fc_c = fc_a;
fm_c = v_c * fc_c/c;
fD_c = fm_c.*cos(phase);
fD_c_nor = fD_c./(max(fD_c)); %經過normalized
sequence_c = round(fD_c_nor,3);
prob_c = zeros(1,2001);

%% Simulation results
for i = 1:2001
    p = (i-1000)/1000;
    prob_a(i) = sum((sequence_a == p));
    prob_b(i) = sum((sequence_b == p));
    prob_c(i) = sum((sequence_c == p));
end

Prob_a = prob_a/t;
Prob_b = prob_b/t;
Prob_c = prob_c/t;

x = -1:0.001:1;

% figure
figure(1);
plot(x,Prob_a);
grid on;
title('PDF of the observed Doppler shift');
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(2);
histogram(fD_a,1000000,'Normalization','probability');
grid on;
xlabel('Doppler Frequency(f)');
ylabel('Probability');

figure(3);
cdfplot(sequence_a);
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(4);
histogram(fD_a,1000000,'Normalization','cdf');
grid on;
ylim([0,1]);
xlabel('Doppler Frequency(f)');
ylabel('Probability');

figure(5);
plot(x,Prob_b);
grid on;
title('PDF of the observed Doppler shift');
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(6);
histogram(fD_b,1000000,'Normalization','probability');
grid on;
xlabel('Doppler Frequency(f)');
ylabel('Probability');

figure(7);
cdfplot(sequence_b);
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(8);
histogram(fD_b,1000000,'Normalization','cdf');
grid on;
ylim([0,1]);
xlabel('Doppler Frequency(f)');
ylabel('Probability');

figure(9);
plot(x,Prob_c);
grid on;
title('PDF of the observed Doppler shift');
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(10);
histogram(fD_c,10000,'Normalization','probability');
grid on;
xlabel('Doppler Frequency(f)');
ylabel('Probability');

figure(11)
cdfplot(sequence_c);
xlabel('Normalized Doppler Frequency(f/fm)');
ylabel('Probability');

figure(12);
histogram(fD_c,10000,'Normalization','cdf');
grid on;
ylim([0,1]);
xlabel('Doppler Frequency(f)');
ylabel('Probability');

%% Theoretical results
% d小題部分
n = 1999;
v_d = v_a;
fc_d = fc_a;
fm_d = fm_a;
w = linspace(-fm_d, fm_d, n);
x_n = -0.999:0.001:0.999;
fD_d = 1./(pi*sqrt(1-x_n.^2)*fm_d);
FD_d = fD_d./(sum(fD_d));

cdf(n) = 0;
cdf(1) = fD_d(2)/(sum(fD_d));
for i = 2:1999
    cdf(i) = cdf(i-1)+fD_d(i)/(sum(fD_d));
end

figure(13);
plot(w,FD_d);
grid on;
title('Theoretical PDF of (a)');
xlabel('Theoretical Doppler Shift of (a)');

figure(14);
plot(w,cdf);
grid on;
title('Theoretical CDF of (a)');
xlabel('Theoretical Doppler Shift of (a)');
ylabel('Theoretical Probability')