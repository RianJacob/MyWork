% Adaptive Equalizer 
close all
clear all

j = sqrt(-1);

% SNR at the output of the channel (dB scale)
SNR = 5:1:30;
%equalizer 
%number of feedforward taps
Nf = 10;
%number of feedback taps
Q = 2;
Delta = 7;

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
        sigma_v_tt = sqrt((norm(channel)^2)/10^(SNR(i)/10));
        %Noise samples
        v = zeros(1,N);
        v(1:training) = sigma_v_tt * 1/sqrt(2)*(randn(1,training) + j*randn(1,training));
        v(training+1:N) = sigma_v_dd * 1/sqrt(2) * (randn(1, decision_directed) + j*randn(1, decision_directed));

        %filter the data through the channel
        y = filter(channel,1,s);
        % additive white Gaussian noise
        r = y + v;

        %initializing the equilizier 
        % equlizier coefficients - feedforward (column vector)
        wf = zeros(Nf,1);
        %equlizer coefficients - feedback filter (column vector)
        wb = zeros(Q,1);
        w = [wb;wf];
        %regressor vector (uf: feedforward, ub:feedback)
        uf = zeros(1,Nf);
        ub = zeros(1,Q);
        %error vector
        e = zeros(1,N);
        %number of errors
        num_errors = 0;

        %Adaptive Equalization

        %Training mode
        % note that this will access the first delta (15) symbol of
        % actual data (decision-directed data) 
        s_hat = zeros(1,N);
        for n=1:training+Delta
            uf = [r(n) uf(1:Nf-1)];
            s_hat(n) = uf*wf + ub*wb;
            d = training_data(n);
            e(n) = d-s_hat(n);
            % u
            u = [ub, uf];
            w = w + mu*u'*e(n)/(u*u'+epsilon);
            wb = w(1:Q);
            wf = w(Q+1:Q+Nf);
        end

        %Decision-directed mode
        s_check = zeros(1,N);
        
        %4-QAM
        if k==1 

            %Decision-directed mode
            s_check = zeros(1,N);
            for n=training+Delta+1:N
                uf = [r(n) uf(1:Nf-1)];
                s_hat(n) = uf*wf + ub*wb;
                s_check(n) = QAM_4_demod(s_hat(n));
            
                d = s_check(n);
                e(n) = d - s_hat(n);
                u = [ub, uf];
                w = w + mu*e(n)*u'/(u*u'+epsilon);
                wb = w(1:Q);
                wf = w(Q+1:Q+Nf);
            
                if (s_check(n) ~= s(n-Delta))
                    num_errors = num_errors + 1;
                end
            end

        %16-QAM
        elseif k==2

            for n=training+Delta+1:N
                uf = [r(n) uf(1:Nf-1)];
                s_hat(n) = uf*wf + ub*wb;
                s_check(n) = QAM_16_demod(s_hat(n));
            
                d = s_check(n);
                e(n) = d - s_hat(n);
                u = [ub, uf];
                w = w + mu*e(n)*u'/(u*u'+epsilon);
                wb = w(1:Q);
                wf = w(Q+1:Q+Nf);
            
                if (s_check(n) ~= s(n-Delta))
                    num_errors = num_errors + 1;
                end
            end

        elseif k==3

            for n=training+Delta+1:N
                uf = [r(n) uf(1:Nf-1)];
                s_hat(n) = uf*wf + ub*wb;
                s_check(n) = QAM_64_demod(s_hat(n));
            
                d = s_check(n);
                e(n) = d - s_hat(n);
                u = [ub, uf];
                w = w + mu*e(n)*u'/(u*u'+epsilon);
                wb = w(1:Q);
                wf = w(Q+1:Q+Nf);
            
                if (s_check(n) ~= s(n-Delta))
                    num_errors = num_errors + 1;
                end
            end
        else
            for n=training+Delta+1:N
                uf = [r(n) uf(1:Nf-1)];
                s_hat(n) = uf*wf + ub*wb;
                s_check(n) = QAM_256_demod(s_hat(n));
            
                d = s_check(n);
                e(n) = d - s_hat(n);
                u = [ub, uf];
                w = w + mu*e(n)*u'/(u*u'+epsilon);
                wb = w(1:Q);
                wf = w(Q+1:Q+Nf);
            
                if (s_check(n) ~= s(n-Delta))
                    num_errors = num_errors + 1;
                end
            end
        end

        Ps(k,i) = num_errors/decision_directed;
    end
end


semilogy(SNR,Ps(1,:),'-*', SNR, Ps(2,:), '-o', SNR, Ps(3,:), '-diamond', SNR, Ps(4,:),'-s')

xlabel('SNR(dB)')
title("Symbol error rate with Decision Feedback Equalizer")
legend('4-QAM', '16-QAM', '64-QAM', '256-QAM')
grid on

   









