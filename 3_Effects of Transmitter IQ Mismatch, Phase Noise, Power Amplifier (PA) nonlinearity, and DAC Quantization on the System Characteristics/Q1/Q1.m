
close all
clear all
j = sqrt(-1);
N = 1e3;

%%%%%%%%%%%%%%%%%%%
%%%% Q1.A (a) %%%%%
%%%%%%%%%%%%%%%%%%%

%%%%%%4-QAM
K=1/2; %gain imbalance 
data4_I = randi([1 2],1,N);
data4_Q = randi([1 2],1,N);

x4_ref = [1+j, -1+j, -1-j, 1-j];

map_4 = [-1 1];

xi_4 = map_4(data4_I);
xq_4 = map_4(data4_Q);

 
x1_4 = xi_4 + j*K*xq_4;

%SNR of 20 dB
SNR_20 = 20;
SNR_20_abs = 10^(SNR_20/10);
y4 = x1_4 + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N)); 


y4plot = scatterplot(y4,1,0,'g.');
hold on
scatterplot(x4_ref,1,0,'r*',y4plot)
title('gain imbalance effects K=1/2')



%%%%%%%16-QAM
data16_I = randi([1 4],1,N);
data16_Q = randi([1 4],1,N);
map_16 = [-3 -1 1 3];

xi_16 = map_16(data16_I);
xq_16 = map_16(data16_Q);

x16 = xi_16 + j*K*xq_16;
y16 = x16 + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N)); 

x16_ref = [map_16(1)+j*map_16, map_16(2)+j*map_16, map_16(3)+j*map_16, map_16(4)+j*map_16 ];

y16plot = scatterplot(y16,1,0,'g.');
hold on
scatterplot(x16_ref,1,0,'r*',y16plot)
title('gain imbalance effects K=1/2')

%%%%%%%%%%%%%%%%%%%
%%%% Q1.A (b) %%%%%
%%%%%%%%%%%%%%%%%%%

K1=1; %gain imbalance 
psi_I = pi/8; %phase imbalance

%%%%%%% 4_QAM
xi_4_tld = xi_4;
xq_4_tld = -K1*sin(psi_I)*xi_4 + K1*cos(psi_I)*xq_4;
x4_tld = xi_4_tld + j*xq_4_tld;
y4_tld = x4_tld + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N)); 

y4_tld_plot = scatterplot(y4_tld,1,0,'g.');
hold on 
scatterplot(x4_ref,1,0,'r*',y4_tld_plot)
title('phase imbalance effects \psi_I=\pi/8')

%%%%%%%% 16-QAM
xi_16_tld = xi_16;
xq_16_tld = -K1*sin(psi_I)*xi_16 + K1*cos(psi_I)*xq_16;
x16_tld = xi_16_tld + j*xq_16_tld;
y16_tld = x16_tld + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N)); 
y16_tld_plot = scatterplot(y16_tld,1,0,'g.');
hold on 
scatterplot(x16_ref,1,0,'r*',y16_tld_plot)
title('phase imbalance effects \psi_I=\pi/8')

%%%%%%%%%%%%%%%%%%%
%%%% Q1.A (c) %%%%%
%%%%%%%%%%%%%%%%%%%



Eb_No = 0:2:20;  %range of Eb/No in dB scale
Eb_No_abs = 10.^(Eb_No./10);

BER_4 = zeros(1,length(Eb_No));
BER_16 = zeros(1,length(Eb_No));

norm16 = 10;  % average power for 16-QAM symbol is 10 (normalization factor)

for i=1:length(Eb_No)
    x2_4 = map_4(randi([1 2],1,N)) + j*map_4(randi([1 2],1,N));

    y2_4 = 1/sqrt(2)*x2_4 + 1/sqrt(4*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y2_4_d = sign(real(y2_4)) + j*sign(imag(y2_4));
    BER_4(i) = (sum(real(y2_4_d) ~= real(x2_4)) + sum(imag(y2_4_d) ~= imag(x2_4)))/(2*N);
    
    x2_16 = map_16(randi([1 4],1,N)) + j*map_16(randi([1 4],1,N));
    y2_16 = x2_16 + sqrt(norm16)/sqrt(8*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y2_16_d = QAM16_demod(y2_16);
    BER_16(i) = (sum(real(y2_16_d) ~= real(x2_16)) + sum(imag(y2_16_d) ~= imag(x2_16)))/(4*N);
  
end
figure(5)
semilogy(Eb_No,BER_4,'LineWidth',1.5)
hold on
semilogy(Eb_No,BER_16,'LineWidth',1.5)
xlabel('Eb/No (dB)')
ylabel('BER')
title('BER without IQ imbalance')
grid on
legend('4-QAM', '16-QAM')


%%%%%%%%%BER K=1/2, psi_I = pi/8

BER_4_IQ = zeros(1,length(Eb_No));
BER_16_IQ = zeros(1,length(Eb_No));

for i=1:length(Eb_No)
    %4-QAM
    x3_4_i = map_4(randi([1 2],1,N));
    x3_4_q = map_4(randi([1 2],1,N));
    x3_4_i_tld = x3_4_i;
    x3_4_q_tld = -K*sin(psi_I)*x3_4_i + K*cos(psi_I)*x3_4_q;
    x3_4_tld = x3_4_i_tld + j*x3_4_q_tld;
    y3_4_tld = 1/sqrt(2)*x3_4_tld + 1/sqrt(4*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y3_4_tld_d = sign(real(y3_4_tld)) + j*sign(imag(y3_4_tld));
    BER_4_IQ(i) = (sum(real(y3_4_tld_d) ~= (x3_4_i)) + sum(imag(y3_4_tld_d) ~= (x3_4_q)))/(2*N);
    
    %16-QAM
    x3_16_i = map_16(randi([1 4],1,N));
    x3_16_q = map_16(randi([1 4],1,N));
    x3_16_i_tld = x3_16_i;
    x3_16_q_tld = -K*sin(psi_I)*x3_16_i + K*cos(psi_I)*x3_16_q;
    x3_16_tld = x3_16_i_tld + j*x3_16_q_tld;
    y3_16_tld = x3_16_tld + sqrt(norm16)/sqrt(8*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y3_16_tld_d = QAM16_demod(y3_16_tld);
    BER_16_IQ(i) = (sum(real( y3_16_tld_d) ~= (x3_16_i)) + sum(imag( y3_16_tld_d) ~= (x3_16_q)))/(4*N);

end
figure(6)
semilogy(Eb_No,BER_4_IQ,'LineWidth',1.5)
hold on 
semilogy(Eb_No,BER_16_IQ,'LineWidth',1.5)
xlabel('EB/No (dB)')
ylabel('BER')
title('BER with IQ imbalance')
grid on
legend('4-QAM', '16-QAM')



%%%%%%%%%%%%%%%%%%%
%%%% Q1.B (a) %%%%%
%%%%%%%%%%%%%%%%%%%
                            
df = 20;      %delta_f in Hz
Ts = 1e-6;    %symbol time
%Effects of phase Noise (psi_N4, psi_N16)

%create received signal (transmitted signa + AWGN)
y4_pn = (xi_4 + j*xq_4) + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N));
y16_pn =(xi_16 + j*xq_16) + 1/sqrt(SNR_20_abs)*(randn(1,N) + j*randn(1,N));

%Phase noise
psi_N4 = zeros(1,N);
psi_N16 = zeros(1,N);
%Gaussian noise
vn_4 = sqrt(4*pi*df*Ts)*randn(1,N);
vn_16 = sqrt(4*pi*df*Ts)*randn(1,N);

%Phase noise defined as [psi(n+1) = ps(n)+v(n)], ps(0)=0
for i=1:N-1
    psi_N4(i+1) = psi_N4(i) + vn_4(i);
end

for i=1:N-1
    psi_N16(i+1) = psi_N16(i) + vn_16(i);
end

%multiply the received signal by phase noise (random complex exponential)
y4_tld_pn = y4_pn.*exp(j*psi_N4);
y16_tld_pn = y16_pn.*exp(j*psi_N16);

plot4 = scatterplot(y4_tld_pn,1,0,'g.');
hold on 
scatterplot(x4_ref,1,0,'r*',plot4)
title('phase noise effects')

plot16 = scatterplot(y16_tld_pn,1,0,'g.');
hold on 
scatterplot(x16_ref,1,0,'r*',plot16)
title('phase noise effects')


BER_4_pn = zeros(1,length(Eb_No));
BER_16_pn = zeros(1,length(Eb_No));
for i=1:length(Eb_No)

    x4_4 = map_4(randi([1 2],1,N)) + j*map_4(randi([1 2],1,N));
    y4_pn = 1/sqrt(2)*x4_4 + 1/sqrt(4*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y4_pn_tld = y4_pn.*exp(j*psi_N4);
    y4_pn_tld_d = sign(real(y4_pn_tld)) + j*sign(imag(y4_pn_tld));
    BER_4_pn(i) = (sum(real(y4_pn_tld_d) ~= real(x4_4)) + sum(imag(y4_pn_tld_d) ~= imag(x4_4)))/(2*N);
    
    x4_16 = map_16(randi([1 4],1,N)) + j*map_16(randi([1 4],1,N));
    y16_pn =  x4_16 + sqrt(norm16)/sqrt(8*Eb_No_abs(i))*(randn(1,N) + j*randn(1,N));
    y16_pn_tld = y16_pn.*exp(j*psi_N16);
    y16_pn_tld_d = QAM16_demod(y16_pn_tld);
    BER_16_pn(i) = (sum(real(y16_pn_tld_d) ~= real(x4_16)) + sum(imag(y16_pn_tld_d) ~= imag(x4_16)))/(4*N);
end

figure(9)
semilogy(Eb_No,BER_4_pn,'LineWidth',1.5)
hold on 
semilogy(Eb_No,BER_16_pn,'LineWidth',1.5)
xlabel('Eb/No(dB)')
ylabel('BER')
title('BER with phase noise')
grid on
legend('4-QAM', '16-QAM')






