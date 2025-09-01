% Adaptive Equalizer 
close all
clear all

j = sqrt(-1);
load channel 
c = channel;

% SNR (dB) at the output of the channel
SNR = 40;
equalizer_length = 55;
Delta = 30;

%epsilon NLMS
mu = 0.4;
epsilon = 1e-6;

%RSL
lambda = 0.995;
%define alpha as reciprocal of lambda for convenience later in claculations
alpha = 1/lambda;

% length of training data 
% epislon NLMS uses 2000 symbols, RLS uses 100 symbol
training = [2000, 100];
% length of decision directed data
decision_directed = 5000;
% total number of data
N = training + decision_directed;

% 4-QAM, average power factor
k_factor = sqrt(2);


% algo = 1 (epsilon NLMS), algo=2 (RSL)
for algo=1:2

    training_data = zeros(1,training(algo)+Delta);
    % vector that includes both training data and decision directed data
    s = zeros(1,N(algo));
    
        
    
    %Noise variance     
    % during decision-directed
    sigma_v_dd = sqrt((k_factor^2 * norm(c)^2)/10^(SNR/10));
    % during training 
    sigma_v_tt = sqrt((norm(c)^2)/10^(SNR/10));
    %Noise samples
    v = zeros(1,N(algo));
    v(1:training(algo)) = sigma_v_tt * 1/sqrt(2)*(randn(1,training(algo)) + j*randn(1,training(algo)));
    v(training(algo)+1:N(algo)) = sigma_v_dd * 1/sqrt(2) * (randn(1, decision_directed) + j*randn(1, decision_directed));

    %training data
    QPSK_map = 1/sqrt(2)*[1,-1];
    s(1:training(algo)) = QPSK_map(randi([1,2],1,training(algo))) + j*QPSK_map(randi([1,2],1,training(algo)));
    % QPSK training Data, the first Delta samples are zero
    training_data(Delta+1:training(algo)+Delta) = s(1:training(algo));

    %Decision-directed data
    QAM_4_map = -1:2:1;
    s(training(algo)+1:N(algo)) = QAM_4_map(randi([1,2],1,decision_directed)) + j*QAM_4_map(randi([1,2],1,decision_directed));

    %filter the data through the channel
    y = filter(c,1,s);
    % additive white Gaussian noise
    r = y + v;

    %initializing the equilizier 
    % equlizier coefficients (column vector)
    w = zeros(equalizer_length,1);
    %regressor vector (u)
    u = zeros(1,equalizer_length);
    %error vector
    e = zeros(1,N(algo));
    %number of errors
    num_errors = 0;
    %RLS 
    P = eye(equalizer_length)*1/epsilon;

    %Decision-directed mode preallocation vectors
    s_hat = zeros(1,N(algo));
    s_check = zeros(1,N(algo));

    %Adaptive equaliziation 
        
    %epsilon NLMS (algo = 1)
    if algo==1

        %Training mode
        % note that this will access the first delta (15) symbol of
        % actual data (decision-directed data) 
        for i=1:training(algo)+Delta
            u = [r(i) u(1:equalizer_length-1)];
    
            d = training_data(i);
            e(i) = d-u*w;
            w = w + mu*u'*e(i)/(u*u'+epsilon);
        end
    
        %decision directed mode
        for i=training(algo)+Delta+1:N(algo)
            u = [r(i) u(1:equalizer_length-1)];
            s_hat(i) = u*w;
            s_check(i) = QAM_4_demod(s_hat(i));
    
            d = s_check(i);
            e(i) = d - u*w;
        
            w = w +mu*e(i)*u'/(u*u'+epsilon);
        
            if (s_check(i) ~= s(i-Delta))
                num_errors = num_errors + 1;
            end
        end
    
    %RLS (algo = 2)    
    else
        %Training mode
        for i=1:training(algo)+Delta
            u = [r(i) u(1:equalizer_length-1)];
            
            P = alpha * (P - alpha*P*u'*u*P/(1+alpha*u*P*u'));
            d = training_data(i);
            e(i) = d-u*w;
            w = w + P*u'*e(i);
        end

        %Decision-directed mode
        for i=training(algo)+Delta+1:N(algo)
            u = [r(i) u(1:equalizer_length-1)];
            s_hat(i) = u*w;
            s_check(i) = QAM_4_demod(s_hat(i));
            d = s_check(i);
            e(i) = d - u*w;
            P = alpha * (P - alpha*P*u'*u*P/(1+alpha*u*P*u'));
            w = w + P*u'*e(i);
            if (s_check(i) ~= s(i-Delta))
                num_errors = num_errors + 1;
            end
        end 
    end
        
         if algo == 1 
            subplot(131)
            plot(real(training_data(Delta+1:end)), imag(training_data(Delta+1:end)),'x',LineWidth=1.5)
            xlabel('Re')
            ylabel('Im')
            title('Training Sequence')
            grid on
            xlim([-2 2])
            ylim([-2 2])
            subplot(132)
            plot(real(s_hat(training(algo)+Delta+1:N(algo))), imag(s_hat(training(algo)+Delta+1:N(algo))), '.')
            xlabel('Re')
            ylabel('Im')
            title('Equalizer output, \epsilonNLMS with 2000 symbols')
            grid on
            axis tight
         else
             subplot(133)
             plot(real(s_hat(training(algo)+Delta+1:N(algo))), imag(s_hat(training(algo)+Delta+1:N(algo))), '.')
             xlabel('Re')
             ylabel('Im')
             title('Equalizer output, RLS with 100 training symbols')
             grid on
             axis tight
         end
end

% convolution of impulse response of the channel and equalizer (RSL)
comb = conv(c,w);

%frequency response of the channel
[H,ww] = freqz(c,1,256);
%frequency response of equalizer (RLS)
[W,~] = freqz(w,1,256);
% frequency response of combination (convolution) of channel and equalizer 
[COMB,~] = freqz(comb,1,256);
%%
figure(2)
subplot(231)
plot(0:length(c)-1, c)
xlabel('Tab index')
ylabel('Amplitude')
title('Impulse res. of channel')
xlim([0 length(c)-1])
grid on

subplot(232)
plot(0:length(w)-1, real(w))
xlabel('Tab index')
ylabel('Amplitude')
title('Impulse res. of equalizer (RLS algorithm)')
xlim([0 length(w)-1])
grid on

subplot(233)
plot(0:length(comb)-1,real(comb))
xlabel('Tab index')
ylabel('Amplitude')
title('Impulse res. of combination')
xlim([0 length(comb)-1])
grid on 

subplot(234)
plot(ww, 20*log10(abs(H)))
xlabel('\omega (rad/sample)')
ylabel('Amplitude (dB)')
title('Freq. response of channel')
grid on 

subplot(235)
plot(ww,20*log10(abs(W)))
xlabel('\omega (rad/sample)')
ylabel('Amplitude (dB)')
title('Freq. response of Equalizer')
grid on 

subplot(236)
plot(ww, 20*log10(abs(COMB)))
xlabel('\omega (rad/sample)')
ylabel('Amplitude (dB)')
title('Freq. response of combination')
grid on

       



    




