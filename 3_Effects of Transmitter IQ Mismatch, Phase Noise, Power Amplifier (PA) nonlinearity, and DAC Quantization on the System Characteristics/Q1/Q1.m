clear all
close all

j = sqrt(-1);
N = 1e3;
                            %%%%%%%%%%%%%%%%%%%
                            %%%% Q1.A (a) %%%%%
                            %%%%%%%%%%%%%%%%%%%

%%%%%%4-QAM
K=1/2;
data4_I = randi([1 2],1,N);
data4_Q = randi([1 2],1,N);

map_4 = [-1 1];
xi_4 = zeros(1,N);
xq_4 = zeros(1,N);

for i=1:N
    xi_4(i) = map_4(data4_I(i));
    xq_4(i) = map_4(data4_Q(i));
end
 
x4 = xi_4 + j*K*xq_4;
y4 = awgn(x4,20,'measured');

xmod4 =xi_4 + j*xq_4;

y4plot = scatterplot(y4,1,0,'g.');
hold on
scatterplot(xmod4,1,0,'r*',y4plot)
title('gain imbalance effects K=1/2')

%%%%%%%16-QAM
data16_I = randi([1 4],1,N);
data16_Q = randi([1 4],1,N);
map_16 = [-3 -1 1 3];
xi_16 = zeros(1,N);
xq_16 = zeros(1,N);
for i=1:N
    xi_16(i) = map_16(data16_I(i));
    xq_16(i) = map_16(data16_Q(i));
end
x16 = xi_16 + j*K*xq_16;
y16 = awgn(x16,20,'measured');

xmod16 = xi_16 + j*xq_16;

y16plot = scatterplot(y16,1,0,'g.');
hold on
scatterplot(xmod16,1,0,'r*',y16plot)
title('gain imbalance effects K=1/2')

                            %%%%%%%%%%%%%%%%%%%
                            %%%% Q1.A (b) %%%%%
                            %%%%%%%%%%%%%%%%%%%
KK=1;
psi_I = pi/8;

%%%%%%% 4_QAM
xi_4_tld = xi_4;
xq_4_tld = -KK*sin(psi_I)*xi_4 + KK*cos(psi_I)*xq_4;
x4_tld = xi_4_tld + j*xq_4_tld;
y4_tld = awgn(x4_tld,20,'measured');

y4_tld_plot = scatterplot(y4_tld,1,0,'g.');
hold on 
scatterplot(xmod4,1,0,'r*',y4_tld_plot)
title('phase imbalance effects \psi_I=\pi/8')

%%%%%%%% 16-QAM
xi_16_tld = xi_16;
xq_16_tld = -KK*sin(psi_I)*xi_16 + KK*cos(psi_I)*xq_16;
x16_tld = xi_16_tld + j*xq_16_tld;
y16_tld = awgn(x16_tld,20,'measured');
y16_tld_plot = scatterplot(y16_tld,1,0,'g.');
hold on 
scatterplot(xmod16,1,0,'r*',y16_tld_plot)
title('phase imbalance effects \psi_I=\pi/8')

                            %%%%%%%%%%%%%%%%%%%
                            %%%% Q1.A (c) %%%%%
                            %%%%%%%%%%%%%%%%%%%
SNR_dB = 0:2:20;
BER_4 = zeros(1,length(SNR_dB));
BER_16 = zeros(1,length(SNR_dB));

xx4 = qamdemod(xmod4,4,'bin');
xx16 = qamdemod(xmod16,16,'bin');
for i=1:length(SNR_dB)
    yy4 = awgn(xmod4,SNR_dB(i),'measured');
    yy4_demod = qamdemod(yy4,4,'bin');
    [nn,BER_4(i)] = biterr(xx4,yy4_demod);
    
    yy16 = awgn(xmod16,SNR_dB(i),'measured');
    yy16_demod = qamdemod(yy16,16,'bin');
    [nm, BER_16(i)] = biterr(xx16,yy16_demod);
end
figure(5)
semilogy(SNR_dB,BER_4,'LineWidth',1.5)
xlabel('SNR')
ylabel('BER')
title('BER of 4-QAM without IQ imbalance')
grid on
figure(6)
semilogy(SNR_dB,BER_16,'LineWidth',1.5)
xlabel('SNR')
ylabel('BER')
title('BER of 16-QAM without IQ imbalance')
grid on

%%%%%%%%%BER K=1/2, psi_I = pi/8
%%%%4-QAM
xxxi_4_tld = xi_4;
xxxq_4_tld = -K*sin(psi_I)*xi_4 + K*cos(psi_I)*xq_4;
xxx4_tld = xxxi_4_tld + j*xxxq_4_tld;

%%%%16-QAM
xxxi_16_tld = xi_16;
xxxq_16_tld = -K*sin(psi_I)*xi_16 + K*cos(psi_I)*xq_16;
xxx16_tld = xxxi_16_tld + j*xxxq_16_tld;


BER_4_IQ = zeros(1,length(SNR_dB));
BER_16_IQ = zeros(1,length(SNR_dB));

xxx4_tld_demod = qamdemod(xxx4_tld,4,'bin');
xxx16_tld_demod = qamdemod(xxx16_tld,16,'bin');

for i=1:length(SNR_dB)
    yyy4_tld = awgn(xxx4_tld,SNR_dB(i),'measured');
    yyy4_tld_demod = qamdemod(yyy4_tld,4,'bin');
    [mn,BER_4_IQ(i)] = biterr(xxx4_tld_demod,yyy4_tld_demod);
    
    yyy16_tld = awgn(xxx16_tld,SNR_dB(i),'measured');
    yyy16_tld_demod = qamdemod(yyy16_tld,16,'bin');
    [mm,BER_16_IQ(i)] = biterr(xxx16_tld_demod,yyy16_tld_demod);
end
figure(7)
semilogy(SNR_dB,BER_4_IQ,'LineWidth',1.5)
hold on 
semilogy(SNR_dB,BER_16_IQ,'LineWidth',1.5)
xlabel('SNR')
ylabel('BER')
title('BER with IQ imbalance')
grid on
legend('4-QAM', '16-QAM')



                            %%%%%%%%%%%%%%%%%%%
                            %%%% Q1.B (a) %%%%%
                            %%%%%%%%%%%%%%%%%%%
                            
df = 20;      %delta_f in Hz
Ts = 1e-6;    %symbol time

x4_pn = xi_4 + j*xq_4;
x16_pn = xi_16 + j*xq_16;

yy4_pn = awgn(x4_pn,20,'measured');
yy16_pn = awgn(x16_pn,20,'measured');

psi_N = zeros(1,N);
vn = sqrt(4*pi*df*Ts)*randn(1,N);
for i=1:N-1
    psi_N(i+1) = psi_N(i) + vn(i);
end

yy4_tld_pn = yy4_pn.*exp(j*psi_N);
yy16_tld_pn = yy16_pn.*exp(j*psi_N);

plot4 = scatterplot(yy4_tld_pn,1,0,'g.');
hold on 
scatterplot(xmod4,1,0,'r*',plot4)
title('phase noise effects')

plot16 = scatterplot(yy16_tld_pn,1,0,'g.');
hold on 
scatterplot(xmod16,1,0,'r*',plot16)
title('phase noise effects')


BER_4_pn = zeros(1,length(SNR_dB));
BER_16_pn = zeros(1,length(SNR_dB));
for i=1:length(SNR_dB)
    y4_pn = awgn(x4_pn,SNR_dB(i),'measured');
    y16_pn = awgn(x16_pn,SNR_dB(i),'measured');
    
    y4_pn_tld = y4_pn.*exp(j*psi_N);
    y16_pn_tld = y16_pn.*exp(j*psi_N);
    
    y4_pn_tld_demod = qamdemod(y4_pn_tld,4,'bin');
    y16_pn_tld_demod = qamdemod(y16_pn_tld,16,'bin');
    
    [uv,BER_4_pn(i)] = biterr(xx4,y4_pn_tld_demod);
    [vv,BER_16_pn(i)] = biterr(xx16,y16_pn_tld_demod);
end
figure(10)
semilogy(SNR_dB,BER_4_pn,'LineWidth',1.5)
hold on 
semilogy(SNR_dB,BER_16_pn,'LineWidth',1.5)
xlabel('SNR')
ylabel('BER')
title('BER with phase noise')
grid on
legend('4-QAM', '16-QAM')