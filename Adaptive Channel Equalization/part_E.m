% Decision Feedback Equalizer (DFE)
close all
clear all

j = sqrt(-1);

% SNR (dB) at the output of the channel
SNR = 30;
%feedforward taps 
Nf = [10,20];
%feedback taps
Q = 2;
Delta = [7,10];

mu = 0.4;
epsilon = 1e-6;


% length of training data 
training = 500;
% length of decision directed data
decision_directed = 5000;
% total number of data
N = training + decision_directed;

% iteration 1: feedfooward tabs = 10, Delta = 7
% iteration 2: feedforward tabs = 20, Delta = 10
% index used to create subplots
index = 1;
for iter=1:2
    training_data = zeros(1,training+Delta(iter));
    % vector that includes both training data and decision directed data
    s = zeros(1,N);

    % Channel
    channel = [0.5, 1.2, 1.5, -1];
    
    
    %64-QAM, factor to average the power
    k_factor = sqrt(42);
    
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
    training_data(Delta(iter)+1:training+Delta(iter)) = s(1:training);
    
    %Decision-directed data
    QAM_64_map = -7:2:7;
    s(training+1:N) = QAM_64_map(randi([1,8],1,decision_directed)) + j* QAM_64_map(randi([1,8],1,decision_directed));
    
    %filter the data through the channel
    y = filter(channel,1,s);
    % additive white Gaussian noise
    r = y + v;

    %initializing the equilizier 
    % equlizier coefficients - feedforward (column vector)
    wf = zeros(Nf(iter),1);
    %equlizer coefficients - feedback filter (column vector)
    wb = zeros(Q,1);
    w = [wb;wf];
    %regressor vector (uf: feedforward, ub:feedback)
    uf = zeros(1,Nf(iter));
    ub = zeros(1,Q);
    %error vector
    e = zeros(1,N);
    %number of errors
    num_errors = 0;
    
    %Adaptive equaliziation 
    
    %Training mode

    s_hat = zeros(1,N);
    for i=1:training+Delta(iter)
        uf = [r(i) uf(1:Nf(iter)-1)];
        s_hat(i) = uf*wf + ub*wb;
        d = training_data(i);
        e(i) = d-s_hat(i);
        % u
        u = [ub, uf];
        w = w + mu*u'*e(i)/(u*u'+epsilon);
        wb = w(1:Q);
        wf = w(Q+1:Q+Nf(iter));
    
    end
    
    %Decision-directed mode
    s_check = zeros(1,N);
    for i=training+Delta(iter)+1:N
        uf = [r(i) uf(1:Nf(iter)-1)];
        s_hat(i) = uf*wf + ub*wb;
        s_check(i) = QAM_64_demod(s_hat(i));
    
        d = s_check(i);
        e(i) = d - s_hat(i);
        u = [ub, uf];
        w = w + mu*e(i)*u'/(u*u'+epsilon);
        wb = w(1:Q);
        wf = w(Q+1:Q+Nf(iter));
    
        if (s_check(i) ~= s(i-Delta))
            num_errors = num_errors + 1;
        end
    end

if iter==1

    subplot(2,2,index)
    plot(real(r(training+Delta(iter)+1:N)), imag(r(training+Delta(iter)+1:N)),'.')
    grid
    xlabel('Re')
    ylabel('Im')
    title('Received Sequence')
    
    subplot(2,2,index+1)
    plot(real(s_hat(training+Delta(iter)+1:N)), imag(s_hat(training+Delta(iter)+1:N)), '.')
    grid 
    xlabel('Re')
    ylabel('Im')
    title('Decision feedback Equalizer (10 feedforward tabs, 2 feedback taps)')
    index =  index+2;
else
    subplot(2,2,index)
    plot(real(r(training+Delta(iter)+1:N)), imag(r(training+Delta(iter)+1:N)),'.')
    grid
    xlabel('Re')
    ylabel('Im')
    title('Received Sequence')

    subplot(2,2,index+1)
    plot(real(s_hat(training+Delta(iter)+1:N)), imag(s_hat(training+Delta(iter)+1:N)), '.')
    grid 
    xlabel('Re')
    ylabel('Im')
    title('Decision feedback Equalizer (20 feedforward tabs, 2 feedback taps)')
end
end


