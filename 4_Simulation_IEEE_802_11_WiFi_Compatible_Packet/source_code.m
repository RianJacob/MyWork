clear all
close all


j = sqrt(-1);

Sk = sqrt(1/2)*[0, 0, 1+j, 0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, -1-j,...
                0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, 0, 0, 0, 0, ...
                -1-j, 0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, 1+j, 0,...
                0, 0, 1+j, 0,0,0, 1+j, 0, 0];
            
 df = 312.5e3;  % frequency separation in Hz
 fs = 20e6;     % sampling frequency in Hz;
 Ts = 1/fs;     %sampling time in sec.
 N = 64;
 
 
 x_STF = zeros(1,160);
 
 Sk_index = -26:1:26;
 
 
 for n=0:159
     xx = 0;
     for k=1:length(Sk)
         
         xx = xx + Sk(k)*exp(j*2*pi*Sk_index(k)*df*n*Ts);
     end
     x_STF(n+1) = xx/sqrt(12);
 end
 

 
 Lk = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, ...
       -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, ...
       -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
 
 Lk_index = -26:1:26;
 x_LTF = zeros(1,160);
 
 for n=0:159
     summ = 0;
     for k=1:length(Lk)
         summ = summ + Lk(k)*exp(j*2*pi*Lk_index(k)*df*n*Ts);
     end
     
     x_LTF(n+1) = summ/sqrt(52);
 end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Q1 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_STF = x_STF(1:16);
figure(1)
plot(0:15,abs(z_STF))
xlabel('n')
ylabel('|z-STF|')
title('magnitude of z-STF')
grid on

figure(2)
plot(0:15,angle(z_STF))
xlabel('n')
ylabel('angle of z-STF')
title('angle of z-STF')
grid on

x = [x_STF x_LTF];

snr = 20;  % snr in dB

% y = x+w(n), since h=[1] which is delta function 
y = awgn(x,snr,'measured');   

[r,n_index] = xcorr(conj(y),z_STF);



SNR = 0:2:20;    %dB

%this matrix is used to hold the indices where exactly 
% we have made detection
peak_index = zeros(1,10);

% array is used to evaluate probability of detection
array = zeros(1,length(SNR));

for i=1:length(SNR)
    counter = 0;
    for m=1:1e4
    y = awgn(x,SNR(i),'measured');
    r = xcorr(conj(y),z_STF);
    
    %threshold 
    th = mean(abs(r)) + 1.08*std(abs(r));
    
    for k=1:length(r)
    if (r(k)>th)
        r(k) = 1;
    else 
        r(k)=0;
    end
    end
    
    %index used to go through peak_index matrix
    s=1;
    for q=1:length(r)
        if (r(q)==1)
            peak_index(s)=q;
            s=s+1;
        end
    end
    a = diff(peak_index);
    if (sum(a)>=9*16)
        counter=counter+1;
    end
   
    end
    array(i) = counter;
end


p = 1-array./(1e4);
figure(3)
semilogy(SNR,p)
xlabel('SNR (dB)')
ylabel('Detection Probability')
title('Probability of Packet Detection')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Q2 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_f = 30e3;    %frequency offset (only in Q.2)


% this loop is doing the multiplication 
%of x_STF by exp(j*2Pi*delta_f*n*Ts) 
for n=0:159
    x_STF(n+1) = x_STF(n+1)*exp(j*2*pi*delta_f*n*Ts);
end


M = 1e4;  % number of Monte-Carlo Simulation

% Matrix used to hold delta_f_estimate
delta_f_est = zeros(length(SNR),M);

%Y_STF matrix is used to hold values of y_STF at each value of SNR
Y_STF = zeros(length(SNR),160);

%Matrix is used to hold square error (delta_f_estimate - delta_f).^2
square_error = zeros(length(SNR),M);

%Matrix is used to hold RMSE values to be ploted 
df_error = zeros(length(SNR));

for i=1:length(SNR)
    for m=1:M
        y = awgn(x_STF,SNR(i),'measured');
        Y_STF(i,:) = y;
        y_STF = Y_STF(i,:);
        EST = 0;
        for n=0:143
            EST = EST + angle(conj(y_STF(n+1+16))*y_STF(n+1));
        end
        delta_f_est(i,m) = EST/(-144*32*pi*Ts);
        square_error(i,m) = (delta_f_est(i,m) - delta_f)^2;
    end
    
    df_error(i) = sqrt(mean(square_error(i,:)));
       
        
end

figure(4)
plot(SNR,df_error)
xlabel('SNR')
ylabel('RMSE')
title('root mean square error estimate of frequency')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Q3 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M3 = 1e4;           %Monte_carlo Simulation 
h = [1 0.5];       %impulse response of channel
Hk = fft(h,N);     %taking fft of H with N=64
xh = conv(x,h);    %convolution of h and x

% Matrix to hold estimated values of Lk_10_est for M3=1e4 Monte-Carlo
Lk_10_est = zeros(length(SNR),M3);
% Matrix to hold estimated values of Lk_22_est for M3=1e4 Monte-Carlo
Lk_22_est = zeros(length(SNR),M3);

% Matrix to hold estimated values of Hk_10_est for M3=1e4 Monte-Carlo
Hk_10_est = zeros(length(SNR),M3);
% Matrix to hold estimated values of Hk_22_est for M3=1e4 Monte-Carlo
Hk_22_est = zeros(length(SNR),M3);

% Matrix to hold square error values of Lk_10_est for M3=1e4 Monte-Carlo
Hk_10_error = zeros(length(SNR),M3);
% Matrix to hold square error values of Lk_20_est for M3=1e4 Monte-Carlo
Hk_22_error = zeros(length(SNR),M3);

% Matrix to hold mean of square error values of Lk_10_est 
%for M3=1e4 Monte-Carlo
Hk_10_error_mean = zeros(1,length(SNR));
% Matrix to hold mean of square error values of Lk_22_est 
%for M3=1e4 Monte-Carlo
Hk_22_error_mean = zeros(1,length(SNR));


for i=1:length(SNR)
    
    for m=1:M3
        y3 = awgn(xh,SNR(i),'measured');
        y3_LTF = y3(161:320);
        vn = y3_LTF(33:96);     %y_LTF samples from 33 to 96
        
        temp_10 = 0;
        temp_22 = 0;
        for n=0:63
             temp_10 = temp_10 + vn(n+1)*exp(-j*2*pi*10*df*n*Ts);
             temp_22 = temp_22 + vn(n+1)*exp(-j*2*pi*22*df*n*Ts);
        end
        
        Lk_10_est(i,m) = temp_10/8;
        Lk_22_est(i,m) = temp_22/8;
        
        Hk_10_est(i,m) = Lk_10_est(i,m)/Lk(26+1+10);
        Hk_22_est(i,m) = Lk_22_est(i,m)/Lk(26+1+22);
        
        Hk_10_error(i,m) = ((abs(Hk_10_est(i,m) - Hk(10)))^2)/(abs(Hk(11)))^2;
        Hk_22_error(i,m) = ((abs(Hk_22_est(i,m) - Hk(22)))^2)/(abs(Hk(23)))^2;
        
    end
    
    Hk_10_error_mean(i) = mean(Hk_10_error(i,:));
    Hk_22_error_mean(i) = mean(Hk_22_error(i,:));
end

figure(5)
plot(SNR,Hk_10_error_mean)
hold on 
plot(SNR,Hk_22_error_mean)
legend('k=10','k=22')
xlabel('SNR (dB)')
ylabel('mean of normalized square estimation error');
title('mean of normalized square estimation error of Hk')
grid on






