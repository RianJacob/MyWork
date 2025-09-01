% Adaptive Equalizer 
close all
clear all

j = sqrt(-1);

% SNR (dB) at the output of the channel
SNR = 30;
equalizer_length = 35;
Delta = 15;

mu = 0.4;
epsilon = 1e-6;

% length of training data 
training = 500;
% length of decision directed data
decision_directed = 5000;
% total number of data
N = training + decision_directed;

training_data = zeros(1,training+Delta);
% vector that includes both training data and decision directed data
s = zeros(1,N);

% Channel
channel = [0.5, 1.2, 1.5, -1];


%256-QAM, factor to average the power
k_factor = sqrt(85);

%Noise variance 
% during decision-directed
sigma_v_dd = sqrt((k_factor^2 * norm(channel)^2)/10^(SNR/10));
% during training 
sigma_v_tt = sqrt((norm(channel)^2)/10^(SNR/10));
%Noise samples
v = zeros(1,N);
v(1:training) = sigma_v_tt * 1/sqrt(2)*(randn(1,training) + j*randn(1,training));
v(training+1:N) = sigma_v_dd * 1/sqrt(2) * (randn(1, decision_directed) + j*randn(1, decision_directed));

%training data
QPSK_map = 1/sqrt(2)*[1,-1];
s(1:training) = QPSK_map(randi([1,2],1,training)) + j*QPSK_map(randi([1,2],1,training));
% QPSK training Data, the first Delta samples are zero
training_data(Delta+1:training+Delta) = s(1:training);

%Decision-directed data
QAM_256_map = -15:2:15;
s(training+1:N) = QAM_256_map(randi([1,16],1,decision_directed)) + j* QAM_256_map(randi([1,16],1,decision_directed));

%filter the data through the channel
y = filter(channel,1,s);
% additive white Gaussian noise
r = y + v;

%initializing the equilizier 
% equlizier coefficients (column vector)
w = zeros(equalizer_length,1);
%regressor vector (u)
u = zeros(1,equalizer_length);
%error vector
e = zeros(1,N);
%number of errors
num_errors = 0;

%Adaptive equaliziation 

%Training mode
% note that this will access the first delta (15) symbol of
% actual data (decision-directed data) 
for i=1:training+Delta
    u = [r(i) u(1:equalizer_length-1)];

    d = training_data(i);
    e(i) = d-u*w;
    w = w + mu*u'*e(i)/(u*u'+epsilon);
end

%Decision-directed mode
s_hat = zeros(1,N);
s_check = zeros(1,N);
for i=training+Delta+1:N
    u = [r(i) u(1:equalizer_length-1)];
    s_hat(i) = u*w;
    s_check(i) = QAM_256_demod(s_hat(i));

    d = s_check(i);
    e(i) = d - u*w;

    w = w +mu*e(i)*u'/(u*u'+epsilon);

    if (s_check(i) ~= s(i-Delta))
        num_errors = num_errors + 1;
    end
end

subplot(221)
plot(real(training_data(Delta+1:training)), imag(training_data(Delta+1:training)),'x',LineWidth=1.5)
grid 
xlabel('Re')
ylabel('Im')
title('Training Sequence')

subplot(222)
plot(real(s(training+Delta+1:N)), imag(s(training+Delta+1:N)), 'x', LineWidth=1.5)
xlim([-20 20])
ylim([-20 20])
grid
xlabel('Re')
ylabel('Im')
title('Transmitted Sequence')

subplot(223)
plot(real(r(training+Delta+1:N)), imag(r(training+Delta+1:N)),'.')
grid
xlabel('Re')
ylabel('Im')
title('Received Sequence')

subplot(224)
plot(real(s_hat(training+Delta+1:N)), imag(s_hat(training+Delta+1:N)), '.')
grid 
xlabel('Re')
ylabel('Im')
title('Equalizer Output')


