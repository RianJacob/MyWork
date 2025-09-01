% Adaptive Equalizer 
close all
clear all

j = sqrt(-1);

% SNR (dB) at the output of the channel
SNR = 30;
equalizer_length = 35;
Delta = 15;

%epsilon NLMS
mu = 0.4;
epsilon = 1e-6;

%LMS
mu1 = 0.001;

% length of training data 
training = [150, 300, 500];
% length of decision directed data
decision_directed = 5000;
% total number of data
N = training + decision_directed;

% index used to access all the subplots
index = 1;
% algo = 1 (LMS), algo=2 (epsilon NLMS)
for algo=1:2

    for k=1:3
    
        training_data = zeros(1,training(k)+Delta);
        % vector that includes both training data and decision directed data
        s = zeros(1,N(k));
    
        % Channel
        channel = [0.5, 1.2, 1.5, -1];
    
    
        %16-QAM, factor to average the power
        k_factor = sqrt(10);

        %Noise variance     
        % during decision-directed
        sigma_v_dd = sqrt((k_factor^2 * norm(channel)^2)/10^(SNR/10));
        % during training 
        sigma_v_tt = sqrt((norm(channel)^2)/10^(SNR/10));
        %Noise samples
        v = zeros(1,N(k));
        v(1:training(k)) = sigma_v_tt * 1/sqrt(2)*(randn(1,training(k)) + j*randn(1,training(k)));
        v(training(k)+1:N(k)) = sigma_v_dd * 1/sqrt(2) * (randn(1, decision_directed) + j*randn(1, decision_directed));
    
        %training data
        QPSK_map = 1/sqrt(2)*[1,-1];
        s(1:training(k)) = QPSK_map(randi([1,2],1,training(k))) + j*QPSK_map(randi([1,2],1,training(k)));
        % QPSK training Data, the first Delta samples are zero
        training_data(Delta+1:training(k)+Delta) = s(1:training(k));

        %Decision-directed data
        QAM_16_map = -3:2:3;
        s(training(k)+1:N(k)) = QAM_16_map(randi([1,4],1,decision_directed)) + j* QAM_16_map(randi([1,4],1,decision_directed));
    
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
        e = zeros(1,N(k));
        %number of errors
        num_errors = 0;

         %Decision-directed mode
         s_hat = zeros(1,N(k));
         s_check = zeros(1,N(k));
    
        %Adaptive equaliziation 
        
        %LMS (algo = 1)
        if algo==1

            % Training mode
            for i=1:training(k)+Delta
                u = [r(i) u(1:equalizer_length-1)];
                d = training_data(i);
                e(i) = d-u*w;
                w = w + mu1*u'*e(i);
            end

            %decision-directed mode
            for i=training(k)+Delta+1:N(k)
                u = [r(i) u(1:equalizer_length-1)];
                s_hat(i) = u*w;
                s_check(i) = QAM_16_demod(s_hat(i));
                d = s_check(i);
                e(i) = d-u*w;
                w = w + mu1*u'*e(i);

                if (s_check(i) ~= s(i-Delta))
                    num_errors = num_errors + 1;
                end
            end
        %epsilon NLMS (algo = 2)    
        else
            %Training mode
            % note that this will access the first delta (15) symbol of
            % actual data (decision-directed data) 
            for i=1:training(k)+Delta
                u = [r(i) u(1:equalizer_length-1)];
        
                d = training_data(i);
                e(i) = d-u*w;
                w = w + mu*u'*e(i)/(u*u'+epsilon);
            end
    
           
            for i=training(k)+Delta+1:N(k)
                u = [r(i) u(1:equalizer_length-1)];
                s_hat(i) = u*w;
                s_check(i) = QAM_16_demod(s_hat(i));
        
                d = s_check(i);
                e(i) = d - u*w;
        
                w = w +mu*e(i)*u'/(u*u'+epsilon);
        
                if (s_check(i) ~= s(i-Delta))
                    num_errors = num_errors + 1;
                end
            end
        end

        sgtitle('Equalizier output')
        subplot(2,3,index)

        if algo ==1 && k==1
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            xlabel('Re')
            ylabel('Im')
            title('LMS, 150 training symbols')
        elseif algo ==1 && k==2
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            xlabel('Re')
            ylabel('Im')
            title('LMS, 300 training symbols')
        elseif algo ==1 && k==3
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            xlabel('Re')
            ylabel('Im')
            title('LMS, 500 training symbols')

        elseif algo ==2 && k==1
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            grid 
            xlabel('Re')
            ylabel('Im')
            title('\epsilonNLMS, 150 training symbols')
        elseif algo ==2 && k==2
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            grid 
            xlabel('Re')
            ylabel('Im')
            title('\epsilonNLMS, 300 training symbols')
        elseif algo ==2 && k==3
            plot(real(s_hat(training(k)+Delta+1:N(k))), imag(s_hat(training(k)+Delta+1:N(k))), '.')
            grid 
            xlabel('Re')
            ylabel('Im')
            title('\epsilonNLMS, 500 training symbols')

        end


        index = index + 1;
    end


end

