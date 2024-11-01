
close all
clear all
j = sqrt(-1);

Sk = sqrt(1/2)*[0, 0, 1+j, 0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, -1-j,...
                0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, 0, 0, 0, 0, ...
                -1-j, 0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, 1+j, 0,...
                0, 0, 1+j, 0,0,0, 1+j, 0, 0];
            
  
 fs = 20e6;     % sampling frequency in Hz;
 Ts = 1/fs;     % sampling time in sec.
 N = 64;        % number of subcarriers 
 df = fs/N;     % frequency separation in Hz
 
 x_STF = zeros(1,160);
 
 Sk_index = -26:1:26;
 
 
x1_STF = zeros(1,160);
 for n=1:160
     % the value of index n starts at zero inside the complex exponential
     x_STF(n) = sum(Sk.*exp(j*2*pi*Sk_index*df*(n-1)*Ts))/sqrt(12);
 end

%{
%The two plots below are for the x_STF
subplot(2,1,1)
plot(0:159,abs(x_STF),'LineWidth',1.5)
xlabel('n')
ylabel('|x-STF|')
title('magnitude of x-STF')
grid on
 
subplot(2,1,2)
plot(0:159,angle(x_STF), 'LineWidth',1.5)
xlabel('n')
ylabel('angle(x-STF)')
title('angle(x-STF)')
grid on
%}
 
 
 Lk = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, ...
       -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, ...
       -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
 
 Lk_index = -26:1:26;
 x_LTF = zeros(1,160);
 
 for n=1:160 
     % the value of index n starts at zero inside the complex exponential
     x_LTF(n) = sum(Lk.*exp(j*2*pi*Lk_index*df*(n-1)*Ts))/sqrt(52);
 end
 
%{
%The two plots below are for the x_LTF
figure
subplot(2,1,1)
plot(0:159,abs(x_LTF), 'LineWidth',1.5)
xlabel('n')
ylabel('|x-LTF|')
title('magnitude of x-LTF')
grid on

subplot(2,1,2)
plot(0:159,angle(x_LTF), 'LineWidth',1.5)
xlabel('n')
ylabel('angle of x-LTF')
title('angle of x-LTF')
grid on
%}


%%%%%%%%%%
%%%%Q1%%%%
%%%%%%%%%%

%z_STF represents one period of x_STF
z_STF = x_STF(1:16);
figure
plot(0:15,abs(z_STF), 'LineWidth',1.5)
xlabel('n')
ylabel('|z-STF|')
title('magnitude of z-STF')
grid on

figure
plot(0:15,angle(z_STF), 'LineWidth',1.5)
xlabel('n')
ylabel('angle of z-STF')
title('angle of z-STF')
grid on

x = [x_STF, x_LTF];

snr_20 = 20;  % snr in dB

% y = x+w(n), since h=[1] which is delta function 
noise_pwr = rms(x)^2/(10^(snr_20/10));
y = x + sqrt(noise_pwr/2)*(randn(size(x))+j*randn(size(x)));

%r is correlation function between y and z_STF used for packet detection 
r = zeros(1,length(x_STF));
m=1:16;
for n=1:145
    r(n) = sum(conj(y(n-1+m)).*z_STF(m));
end


% This figure is needed to come up with the threshold value used for packet
% detection
figure
plot(0:length(r)-1,abs(r), 'LineWidth',1.5)
xlabel('n');
ylabel('r(n)');
title('r(n) at SNR=20 dB')
grid on

SNR = 0:2:20;    %dB
SNR_abs = 10.^(SNR/10);


% array is used to evaluate probability of detection
P_detection = zeros(1,length(SNR));

% u is used as index to access z_STF below in computation of correlation
% between received signal y and z_STF
u=1:length(z_STF);

% max number of iterations 
max_iter = 1e4;

for i=1:length(SNR)
    counter = 0;
    for m=1:max_iter
        noise_power = rms(x)^2/SNR_abs(i);
        y = x+sqrt(noise_power/2)*(randn(size(x))+j*randn(size(x)));

        %correlation between y and z_STF
        for n=1:145
            r(n) = sum(conj(y(n-1+u)).*z_STF(u));
        end
        %threshold (determined based on plot of r(n) at SNR=20 dB
        %threshold is SNR dependent 
        th = mean(abs(r)) + 2*std(abs(r));
        
        % need to find at least 9 peaks to declare packet detection
        a = diff(find(r>th));
        if (sum(a)==9*16)
            counter=counter+1;
        end
   
    end
    P_detection(i) = counter/max_iter;
end


figure
semilogy(SNR,P_detection, 'LineWidth',1.5)
xlabel('SNR (dB)')
ylabel('Probability of Packet Detection ')
title('Probability of Packet Detection')
grid on



%%%%%%%%%%
%%%%Q2%%%%
%%%%%%%%%%

delta_f = 30e3;    %frequency offset (only in Q.2)


%short training field under the effect of frequency offset
n=0:159;
x_STF_f_offset = x_STF.*exp(j*2*pi*delta_f*n*Ts);

M = 1e4;  % number of Monte-Carlo Simulation


%Matrix is used to hold square error (delta_f_estimate - delta_f).^2
square_error = zeros(1,M);

%Matrix is used to hold RMSE values to be ploted 
df_error = zeros(1,length(SNR));

% time index used for correlation below for frequency offset calculation
n=0:143;

for i=1:length(SNR)
    noise_power = rms(x_STF_f_offset)^2/(SNR_abs(i));

    for m=1:M
        noise = sqrt(noise_power/2)*(randn(size(x_STF_f_offset))+j*randn(size(x_STF_f_offset)));
        y = x_STF_f_offset + noise;
        
       
        delta_f_est = sum((angle(conj(y(n+1+16)).*y(n+1))))/(-144*32*pi*Ts);
        square_error(m) = (delta_f_est - delta_f)^2;
    end
    df_error(i) = sqrt(sum(square_error)/M);      
end

figure
%scale by 1e3 to convert to kHz
plot(SNR,df_error/1e3, 'LineWidth',1.5)
xlabel('SNR')
ylabel('RMSE')
title('root mean square error estimate of frequency (kHz)')
grid on




%%%%%%%%%%
%%%%Q3%%%%
%%%%%%%%%%


 
h = [1 0.5];       %impulse response of channel
Hk = fft(h,N);     %taking fft of H with N=64
xh = conv(x,h);    %convolution of h and x


% Matrix to hold estimated values of Hk_10_est for M3=1e4 Monte-Carlo
Hk_10_est = zeros(length(SNR),max_iter);
% Matrix to hold estimated values of Hk_22_est for M3=1e4 Monte-Carlo
Hk_22_est = zeros(length(SNR),max_iter);

% Matrix to hold square error values of Lk_10_est for M3=1e4 Monte-Carlo
Hk_10_error = zeros(length(SNR),max_iter);
% Matrix to hold square error values of Lk_20_est for M3=1e4 Monte-Carlo
Hk_22_error = zeros(length(SNR),max_iter);

% Matrix to hold mean of square error values of Lk_10_est 
%for M3=1e4 Monte-Carlo
Hk_10_error_mean = zeros(1,length(SNR));
% Matrix to hold mean of square error values of Lk_22_est 
%for M3=1e4 Monte-Carlo
Hk_22_error_mean = zeros(1,length(SNR));

%time index used below for Lk_10_est and Lk_22_est
n=0:63;
for i=1:length(SNR)
    noise_power1 = rms(xh)^2/SNR_abs(i);
    
    for m=1:max_iter
        noise = sqrt(noise_power1/2)*(randn(size(xh)) + j*randn(size(xh)));
        y3 = xh+noise;
        y3_LTF = y3(161:320);
        vn = y3_LTF(33:96);     %y_LTF samples from 33 to 96
        
        % n=0:63
        Lk_10_est = sum(vn.*exp(-j*2*pi*10*df*n*Ts))/8;
        Lk_22_est = sum(vn.*exp(-j*2*pi*22*df*n*Ts))/8;
        
        %Hk_10_est(i,m) = Lk_10_est(i,m)/Lk(26+1+10);
        %Hk_22_est(i,m) = Lk_22_est(i,m)/Lk(26+1+22);

        Hk_10_est(i,m) = Lk_10_est/Lk(26+1+10);
        Hk_22_est(i,m) = Lk_22_est/Lk(26+1+22);
        
        Hk_10_error(i,m) = ((abs(Hk_10_est(i,m) - Hk(11)))^2)/(abs(Hk(11)))^2;
        Hk_22_error(i,m) = ((abs(Hk_22_est(i,m) - Hk(23)))^2)/(abs(Hk(23)))^2;
        
    end
    
    Hk_10_error_mean(i) = mean(Hk_10_error(i,:));
    Hk_22_error_mean(i) = mean(Hk_22_error(i,:));
end

figure
plot(SNR,Hk_10_error_mean, 'LineWidth',1.5)
hold on 
plot(SNR,Hk_22_error_mean, 'LineWidth',1.5)
legend('k=10','k=22')
xlabel('SNR (dB)')
ylabel('mean of normalized square estimation error');
title('mean of normalized square estimation error of Hk')
grid on


num_sym = 12*64;  % we need 12 OFDM block
map = [-1 1];     % map used to create 4-QAM signal

%Generate 4-QAM signal
tx_symb = 1/sqrt(2)*(map(randi([1 2],max_iter,num_sym)) + j*map(randi([1 2],max_iter,num_sym)));


%qam4_sym = unique(tx_symb);  % 4 unique symbols of 4-QAM signal 
sym_error_k10 = zeros(max_iter,length(SNR));
sym_error_k22 = zeros(max_iter,length(SNR));

for iter = 1:max_iter
    x_data = zeros(1,960);
    for sym_idx = 1:12    %index used to go through 12 OFDM blocks
        useful_sym = ifft(tx_symb((sym_idx-1)*64+1:sym_idx*64));
        useful_sym = useful_sym/rms(useful_sym);
        
        %create cyclic prefix, the first 16 samples of each OFDM block
        x_data((sym_idx-1)*80+1:(sym_idx-1)*80+16) = useful_sym(end-15:end);
        
        %append the data part of OFDM block
        x_data((sym_idx-1)*80+17:(sym_idx)*80) = useful_sym;
    end
    
    y_data = conv(x_data,h);
    y_data = y_data(1:960);
    
    %Add noise
    for i=1:length(SNR)
        noise_power = 1/(10^(SNR(i))/10);
        y_data_noisy =  y_data + 1/sqrt(2)*sqrt(noise_power)*(randn(1,length(y_data))+j*randn(1,length(y_data)));
        
        % Rx
        num_error_k10 = 0;
        num_error_k22 = 0;
        
        for sym = 1:12
            useful_sym_rx = fft(y_data_noisy((sym-1)*80+17 : sym*80));
            useful_sym_rx = useful_sym_rx/rms(useful_sym_rx);
            
            %Channel equalization 
            
            useful_sym_rx_10 = useful_sym_rx(11)/Hk_10_est(i,iter);
            useful_sym_rx_22 = useful_sym_rx(23)/Hk_22_est(i,iter);
            
            
            %detection 
            useful_sym_detected_10 = 1/sqrt(2)*(sign(real(useful_sym_rx_10)) + j*sign(imag(useful_sym_rx_10)));
            useful_sym_detected_22 = 1/sqrt(2)*(sign(real(useful_sym_rx_22)) + j*sign(imag(useful_sym_rx_22)));
           
            
            %error 
            num_error_k10 = num_error_k10 + sum(useful_sym_detected_10 ~=tx_symb((sym-1)*64+11));
            num_error_k22 = num_error_k22 + sum(useful_sym_detected_22 ~=tx_symb((sym-1)*64+23));
        end
        
        sym_error_k10(iter,i) = num_error_k10;
        sym_error_k22(iter,i) = num_error_k22;
    end
end


figure
plot(0:N-1,abs(Hk), 'LineWidth',1.5)
xlabel('k')
ylabel('|FFT(h)|')
title('64 point |FFT(h)|')
grid on


figure
semilogy(SNR,sum(sym_error_k10(:,:))./(12*M),'LineWidth',1.5)
hold on 
semilogy(SNR,sum(sym_error_k22(:,:))./(12*M),'LineWidth',1.5)
legend('k=10', 'k=22')
xlabel('SNR (dB)')
ylabel('symbol error rate')
title('Probability of symbol error rate')
grid on