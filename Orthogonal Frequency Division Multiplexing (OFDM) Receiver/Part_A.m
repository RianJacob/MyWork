clear all
N = 64;      % OFDM Symbol Size
P = 16;      % size of Cyclic Prefix

M = 8;       % Channel Length
L = 8;       % Training tones in symbol

p_ind = 4:8:60;

% QPSK Modulation scheme
map = 1/sqrt(2) * [-1,1];
j = sqrt(-1);

F = fft(eye(N))/sqrt(N);

x = map(randi([1,2],N,L)) + j*map(randi([1,2],N,L));
x(p_ind,:) = 1/sqrt(2)*(ones(length(p_ind),L) + j*ones(length(p_ind),L));

x_bar = F'*x;   %transform domain

h = (randn(1,M) + j*randn(1,M))/sqrt(2*M);
h_app = [h zeros(1,N-M)];

H = zeros(64,64);
H(1,:) = h_app;
for i=2:N
    H(i,:) = circshift(H(i-1,:),1);
end


SNR = 10;     % dB
var = 0.1;

v_bar = sqrt(var/2) * (randn(N,L) + j*randn(N,L));

y_bar = H*x_bar + v_bar;

y = F*y_bar;

Qa = F(1:M,1:N);
Qb = F(1:M,p_ind);

S = diag(x(:,1));

h_hata = 1/sqrt(N)*inv(Qa*Qa')*Qa*(conj(x(:,1)).*y(:,1));

h_hatb = 1/sqrt(N)*inv(Qb*Qb')*Qb*(conj(x(p_ind,:)).*y(p_ind,:));

h_b = zeros(M,1);
for i=1:M
    h_b(i) = mean(h_hatb(i,:));
end
    
subplot(311)
stem(abs(h))
xlabel('Tap index')
ylabel('Amplitude')
title('Exact channel taps')
grid on

subplot(312)
stem(abs(h_hata))
xlabel('Tap index')
ylabel('Amplitude')
title('Taps estimated by using 64 tones in one symbol')
grid on

subplot(313)
stem(abs(h_b))
xlabel('Tap index')
ylabel('Amplitude')
title('Taps estimated by using 64 tones over 8 symbols')
grid on


h_bar = F'*([h.' ;zeros(N-M,1)]);
figure(2)
subplot(311)
stem(abs(h_bar))
xlabel('Tone index')
title('H (transform domain)')
grid on

h_bar_a = F'*([h_hata;zeros(N-M,1)]);
subplot(312)
stem(abs(h_bar_a))
xlabel('Tone index')
title('Ha')
grid on

h_bar_b = F'*([h_b; zeros(N-M,1)]);
subplot(313)
stem(abs(h_bar_b))
xlabel('Tone index')
title('Hb')
grid on












