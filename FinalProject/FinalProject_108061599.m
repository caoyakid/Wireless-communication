clc;
close all;
clear all;

%%  Parameters setting (look at p.296)
N = input('Please enter the number of antenna elements N:');
delta = input('The antenna separation:');
system = input('SIMO:1 / MISO:2 ');
signal_degree = input('Direction of desiered signal is:');
interference_degree = input('Direction of interference signal is:');
total_length = N*delta;

%% SIMO
if system == 1
%-----Basis-----
omega_r = (0:N-1)/(N*delta); %N antennas : N個入射角 形成1*N matrix
U_r = zeros(N,N); %基底為N*N的矩陣

for L = 0:N-1   %i=0矩陣數值為1 代表最短距離d i=1,2,3 ... 距離為d+(i-1)*delta*omaga_r
    U_r (L+1,:) = exp(-1i*2*pi*L*delta*omega_r)./sqrt(N);
    % 基底(不同天線和不同角度接收的值 % li:表示複數(虛數單位)
end

%-----Correlation-----
omega_rr = -2:0.01:2; %天線彼此間角度差之自相關 x軸範圍-2~2 與cos有關
correlation = abs(sin(pi*total_length*omega_rr)./(N*sin(pi*total_length*omega_rr/N))); % 7.36
figure(1);
plot(omega_rr, correlation,'-b','linewidth',1.5);
xlabel('The value of cosin function');       
ylabel('correlation') ;
title('The correlation between different basis vectors');

%-----Gain Patterns-----
for k=0:(N-1)            
    xx=pi/2:(2*pi)/401:5*pi/2-(2*pi)/401;
    omega_d=cos(xx)-k/(delta*N);           %公式7.63
    corr_d=abs(sin(pi*N*delta*omega_d)./(N*sin(pi*delta*omega_d)));
    figure
    p = polar(xx,corr_d,'b-')
    p.LineWidth = 1.5
    title('The gain pattern of the ULA');
end

%-----The Gain of desired signal-----
signal=zeros(N,1)+1/sqrt(N);

y=zeros(N,1);  %建立一初始選擇的入射角(direcrion of desiered signal)N*1矩陣
for J=0:(N-1)
    y(J+1,1)=exp(-1i*2*pi*J*delta*(cos(pi*signal_degree/180)));
end

y_a=U_r'*y.*signal;   %入射訊號和基底做內積(在哪個天線的投影量最大代表在那根天線的gain最大)

figure
disp('desiered signal gain ')
disp(abs(sum(y_a./y)))
stem(1:N,(abs(y_a)'),'-b','linewidth',1.5)
xlabel('The Number of the antenna');       
ylabel('Gain') ;
title({'The gain of the desired signal for ';'using different radiation/reception beams (SIMO)'});

%-----Signal to interference power ratio (SINR)-----
I=zeros(N,1);
for J=0:(N-1)
    I(J+1,1)=exp(-1i*2*pi*J*delta*(cos(pi*interference_degree/180)-cos(pi*signal_degree/180)));
end

I_a=U_r'*I.*signal;

SINR=sum(y_a./I_a);  %要求功率 需絕對值平方
disp('SINR(dB)')
disp(10*log(abs(SINR)))

% %ingmatrix(SIMO)
%  ingmatrix=zeros(N,1);
%  intrad=2*pi*interference_degree/360;
%  omegarbe=cos(intrad);
%  for x=1:N
%      ingmatrix(x)=exp(-i*2*pi*(x-1)*delta*omegarbe);
%  end
%  interference_gain=U_r*ingmatrix;
%  %SINR(SIMO)
%  SINR=zeros(1,N);
%  signal_totalpower=I_a.*conj(I_a);
%  interference_totalpower=interference_gain.*conj(interference_gain);
%  for x=1:N
%      SINR(x)=signal_totalpower(x)/interference_totalpower(x);
%  end
%  x=1:N;
%  y=10*log10(SINR(x));
%  disp(y)
%  figure;
%  stem(x,y)
%  xlabel('Power of the incoming signal of insert');       
%  ylabel('SINR') ;
%  title('SINR for using different beams');

%-----SINR of multiple input signals (multiple reception directions) with
%diversity combining-----


PP=100000;   %因fading通道具隨機性,所以要多跑幾次SINR取結果平均
SINR_DC=zeros(1,PP);
SINR_EGC=zeros(1,PP);
SINR_MRC=zeros(1,PP);

for D=1:PP
    %N個branch，N種角度(path)-----------------------------------------------
   
    ray_a=(randn(N,1)+1i*randn(N,1))/sqrt(2);  %N*1 
    ray_I=(randn(N,1)+1i*randn(N,1))/sqrt(2);  %因選擇從特定角度入射 所以會有N種path
                                               %就只需產生N個fading channel
    a=signal.*ray_a;
    b=signal.*ray_I;
    %要求功率 需絕對值平方
    
    G_d1=a./exp(1i*angle(ray_a)).*y_a;
    G_I1=b./exp(1i*angle(ray_I)).*I_a;
    
    SINR_EGC(1,D)=abs(sum(G_d1))^2/abs(sum(G_I1))^2;
    
    G_d2=a.*conj(ray_a).*y_a;
    G_I2=b.*conj(ray_I).*I_a;
    
    SINR_MRC(1,D)=abs(sum(G_d2))^2/abs(sum(G_I2))^2;
end
SIR_EGCC=sum(SINR_EGC)/PP;
SIR_MRCC=sum(SINR_MRC)/PP;

disp('Average_SINR_EGC(dB) ')
disp(10*log(abs(SIR_EGCC)))
disp('Average_SINR_MRC(dB) ')
disp(10*log(abs(SIR_MRCC)))

end


%% MISO
if system == 2
%-----Basis-----
omega_r = (0:N-1)/(N*delta); %N antennas : N個入射角 形成1*N matrix
U_t = zeros(N,N); %基底為N*N的矩陣
for L = 0:N-1   %i=0矩陣數值為1 代表最短距離d i=1,2,3 ... 距離為d+(i-1)*delta*omaga_r
    U_t (L+1,:) = exp(-1i*2*pi*L*delta*omega_r)./sqrt(N);
    % 基底(不同天線和不同角度接收的值 % li:表示複數(虛數單位)
end
%-----Correlation-----
omega_rr = -2:0.01:2; %天線彼此間角度差之自相關 x軸範圍-2~2 與cos有關
correlation = abs(sin(pi*total_length*omega_rr)./(N*sin(pi*total_length*omega_rr/N))); % 7.36
figure(1);
plot(omega_rr, correlation,'-b','linewidth',1.5);
xlabel('The value of cosin function');       
ylabel('correlation') ;
title('The correlation between different basis vectors');

%-----Gain Patterns-----
for k=0:(N-1)            
    xx=pi/2:(2*pi)/401:5*pi/2-(2*pi)/401;
    omega_d=cos(xx)-k/(delta*N);           %公式7.63
    corr_d=abs(sin(pi*N*delta*omega_d)./(N*sin(pi*delta*omega_d)));
    figure
    p = polar(xx,corr_d,'b-')
    p.LineWidth = 1.5
    title('The gain pattern of the ULA');
end

%-----The Gain of desired signal-----

y=zeros(N,1);  %建立一初始選擇的入射角(direcrion of desiered signal)N*1矩陣
for J=0:(N-1)
    y(J+1,1)=exp(-1i*2*pi*J*delta*(cos(pi*signal_degree/180)));
end

y_a=U_t'*y;   %入射訊號和基底做內積(在哪個天線的投影量最大代表在那根天線的gain最大)

figure
disp('desiered signal gain ')
disp(abs(sum(y_a./y)))
stem(1:N,(abs(y_a)'),'-b','linewidth',1.5)
xlabel('The Number of the antenna');       
ylabel('Gain') ;
title({'The gain of the desired signal for ';'using different radiation/reception beams (MISO)'});

%-----Signal to interference power ratio (SINR)-----
I=zeros(N,1);
for J=0:(N-1)
    I(J+1,1)=exp(-1i*2*pi*J*delta*cos(pi*interference_degree/180))/sqrt(N);
end

I_a=U_t'*I;

SINR=sum(abs(y_a))/sum(abs(I_a));
disp('SINR(dB) ')
disp(10*log(SINR))

% %imgmatrix(MISO)
% ingmatrix=zeros(N,1);
% intrad=2*pi*interference_degree/360;
% omegarbe=cos(intrad);
%  for x=1:N
%      ingmatrix(x)=exp(-i*2*pi*(x-1)*delta*omegarbe);
%  end
%  interference_gain=U_t*ingmatrix;
%  %SINR(MISO)
%  SINR=zeros(1,N);
%  signal_totalpower=I_a.*conj(I_a);
%  interference_totalpower=interference_gain.*conj(interference_gain);
%  for x=1:N
%      SINR(x)=signal_totalpower(x)/interference_totalpower(x);
%  end
%  x=1:N;
%  y=10*log10(SINR(x));
%  bmatrix=zeros(1,N);
%  bmatrix(1)=sum(y);
%  disp(y)
%  figure;
%  stem(x,bmatrix)
%  xlabel('Power of the incoming signal of insert');       
%  ylabel('SINR') ;
%  title('SINR for using different beams');

end