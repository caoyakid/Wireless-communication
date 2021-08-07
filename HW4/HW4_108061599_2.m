clc;
clear all;
close all;
%% 
N_0 = 1;
E_av = 2;
N = 10e5;
P_e = zeros(4,4,9);

for L = 1:1:4
    for gamma = 1:1:9
        for i = 1:1:N
            r = zeros(4);
            if (randn(1)>=0)
                s=1+1j;
            else
                s=-1-1j;
            end
            
            for LL = 1:1:L
                gain = sqrt((10^(gamma/10)*N_0)/(2*E_av))*randn(1)+1j*sqrt((10^(gamma/10)*N_0)/(2*E_av))*randn(1)+sqrt((10^(gamma/10)*N_0)/(E_av));
                noise = sqrt(N_0/2)*randn(1)+1j*sqrt(N_0/2)*randn(1);
                received=gain*s+noise;
                if(abs(conj(gain)*received)>=abs(r(1)))
                    r(1)=conj(gain)*received;
                end
                r(2)=r(2)+conj(gain)*received;
                r(3)=r(3)+conj(gain)*received/abs(gain);
                r(4)=r(4)+received;
            end
            
            for k=1:1:4
                if (real(r(k))+imag(r(k)))>=0
                    d=1+1j;
                else
                    d=-1-1j;
                end
                if(d~=s)
                    P_e(k,L,gamma)=P_e(k,L,gamma)+1;
                end
            end
        end
    end
end             

P_e=P_e/N;
P1=shiftdim(P_e(1,:,:));
P2=shiftdim(P_e(2,:,:));
P3=shiftdim(P_e(3,:,:));
P4=shiftdim(P_e(4,:,:));
                
gamma=1:1:9;
%% 
figure();
subplot(4,1,1);
semilogy(gamma, P1(1,:),'-o',gamma, P1(2,:),'-^',gamma, P1(3,:),'-v',gamma, P1(4,:),'-s','linewidth',2);
title('Ricean fading (Selective Combining)');
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');
ylabel('BER')
subplot(4,1,2);
semilogy(gamma, P2(1,:),'-o',gamma, P2(2,:),'-^',gamma, P2(3,:),'-v',gamma, P2(4,:),'-s','linewidth',2);
title('Ricean fading (Maximal Ratio Combining)');
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');
ylabel('BER')
subplot(4,1,3);
semilogy(gamma, P3(1,:),'-o',gamma, P3(2,:),'-^',gamma, P3(3,:),'-v',gamma, P3(4,:),'-s','linewidth',2);
title('Ricean fading (Equal Gain Combining)');
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');
ylabel('BER')
subplot(4,1,4);
semilogy(gamma, P4(1,:),'-o',gamma, P4(2,:),'-^',gamma, P4(3,:),'-v',gamma, P4(4,:),'-s','linewidth',2);
title('Ricean fading (Direct Combining)');
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');
ylabel('BER')                