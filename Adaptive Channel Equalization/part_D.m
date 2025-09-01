% Adaptive Equalizer 
close all
clear all

j = sqrt(-1);

% SNR (dB) at the output of the channel (dB scale)
SNR = 5:1:30;
%equalizer 
equalizer_length = 35;
Delta = 15;

mu = 0.4;
epsilon = 1e-6;

%length of training symbols
training = 500;
%length of decision-directed data
decision_directed = 5000;
% total number of data
N = training + decision_directed;

%training data
training_data = zeros(1,training + Delta);
% vector that includes both training data and decision directed data
s = zeros(1,N);

%training data
QPSK_map = 1/sqrt(2)*[1,-1];
s(1:training) = QPSK_map(randi([1,2],1,training)) + j*QPSK_map(randi([1,2],1,training));
% QPSK training Data, the first Delta samples are zero
training_data(Delta+1:training+Delta) = s(1:training);



% Channel
channel = [0.5, 1.2, 1.5, -1];

% QAM constellation size
%Probability of Symbol error
Ps = zeros(4,length(SNR));


for k=1:4
   
    %k = 1 (4-QAM)
    if k==1
        QAM_map = -1:2:1;
        k_factor = sqrt(2);
        s(training+1:N) = QAM_map(randi([1 2],1,decision_directed)) + j*QAM_map(randi([1 2],1,decision_directed));

    %k = 2 (16-QAM)
    elseif k==2
        QAM_map = -3:2:3;
        k_factor = sqrt(10);
        s(training+1:N) = QAM_map(randi([1 4],1,decision_directed)) + j*QAM_map(randi([1 4],1,decision_directed));

    %k = 3 (64-QAM)
    elseif k==3
        QAM_map = -7:2:7;
        k_factor = sqrt(42);
        s(training+1:N) = QAM_map(randi([1 8],1,decision_directed)) + j*QAM_map(randi([1 8],1,decision_directed));

    else
        %k = 4 (256-QAM)
        QAM_map = -15:2:15;
        k_factor = sqrt(170);
        s(training+1:N) = QAM_map(randi([1 16],1,decision_directed)) + j*QAM_map(randi([1 16],1,decision_directed));
    end

    for i=1:length(SNR)

        %Noise variance 
        % during decision-directed
        sigma_v_dd = sqrt((k_factor^2 * norm(channel)^2)/10^(SNR(i)/10));
        % during training 
        sigma_v_tt = sqrt(norm(channel)^2/10^(SNR(i)/10));
        %Noise samples
        v = zeros(1,N);
        % scale by sqrt(2) is needed for normalization because the noise is complex
        v(1:training) = sigma_v_tt * 1/sqrt(2)*(randn(1,training) + j*randn(1,training));
        v(training+1:N) = sigma_v_dd * 1/sqrt(2) * (randn(1, decision_directed) + j*randn(1, decision_directed));
        
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

        %Adaptive Equalization

        %Training mode
        % note that this will access the first delta (15) symbol of
        % actual data (decision-directed data) 
        for n=1:training+Delta
            u = [r(n) u(1:equalizer_length-1)];
            d = training_data(n);
            e(n) = d-u*w;

            w = w + mu*u'*e(n)/(u*u'+epsilon);
        end

        %Decision-directed mode
        s_hat = zeros(1,N);
        s_check = zeros(1,N);
        
        %4-QAM
        if k==1 

            for n=training+Delta+1:N
                u = [r(n) u(1:equalizer_length-1)];
                s_hat(n) = u*w;
                s_check(n) = QAM_4_demod(s_hat(n));

                d = s_check(n);
                e(n) = d-u*w;
                w = w +mu*e(n)*u'/(u*u'+epsilon);
                num_errors = num_errors + (s_check(n) ~= s(n-Delta));
            end
            
        
        %16-QAM
        elseif k==2

            for n=training+Delta+1:N
                u = [r(n) u(1:equalizer_length-1)];
                s_hat(n) = u*w;
                s_check(n) = QAM_16_demod(s_hat(n));

                d = s_check(n);
                e(n) = d-u*w;
                w = w +mu*e(n)*u'/(u*u'+epsilon);
                num_errors = num_errors + (s_check(n) ~= s(n-Delta));  
            end

        elseif k==3

            for n=training+Delta+1:N
                u = [r(n) u(1:equalizer_length-1)];
                s_hat(n) = u*w;
                s_check(n) = QAM_64_demod(s_hat(n));

                d = s_check(n);
                e(n) = d-u*w;
                w = w +mu*e(n)*u'/(u*u'+epsilon);
                num_errors = num_errors + (s_check(n) ~= s(n-Delta));
            end
                 
        else
            for n=training+Delta+1:N
                u = [r(n) u(1:equalizer_length-1)];
                s_hat(n) = u*w;
                s_check(n) = QAM_256_demod(s_hat(n));

                d = s_check(n);
                e(n) = d-u*w;
                w = w + mu*e(n)*u'/(u*u'+epsilon);
                num_errors = num_errors + (s_check(n) ~= s(n-Delta));
            end
        end        
        
        
        Ps(k,i) = num_errors/decision_directed;
    end
end


semilogy(SNR,Ps(1,:),'-o', SNR, Ps(2,:), '-^', SNR, Ps(3,:), '-*', SNR, Ps(4,:),'-s')

xlabel('SNR(dB)')
title("Symbol error rate with feedforward equalizer")
legend('4-QAM', '16-QAM', '64-QAM', '256-QAM')
grid on

   









