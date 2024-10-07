clear all
close all

%Time varying Channel Estimation uisng Kalman Filter

A = [0.99 0; 0 0.999];
Q = 0.001 * eye(2);
N = 102;

% Initialization
h_hat = zeros(2,1);
M = 100*eye(2);


% channel coefficients 
h_0 = zeros(1,N-1);
h_1 = zeros(1,N-1);

Mn = zeros(2,N-1);
Kn = zeros(2,N-1);
H = zeros(2,N-1);


p1 = zeros(1,5);
p2 = ones(1,5);
p = [p1 p2];
V = repmat(p,1,10);
V = [V,1 1];


% Time Varying channel
h1 = exp(-rand([1, N-1]));
h1 = sort(h1);
h1 = h1(end:-1:1);
H(1,:) = h1;
b = 0.85:0.02:0.95;
H(2,:) = b(randi([1 6],1,N-1));

% noiseless output
y = zeros(1,N-1);
% noisyoutput
x1 = zeros(1,N-1);

for n=1:N-1

    hn_1_hat = A*h_hat;
    Mn_1 = A*M*A' + Q;
    

    %generate pilot signal 
    v = [V(n);V(n+1)];
    % generate h (time varying channel)
    h = [H(1,n);H(2,n)];
 
    % generate observation 
    x = v'*h + sqrt(0.1)*randn;
    K = Mn_1 * v./(0.1 + v'*Mn_1*v);
    h_hat = hn_1_hat + K * (x - v'*hn_1_hat);
    M = (eye(2) - K*v')*Mn_1;

    y(n) = v'*h;
    x1(n) = x;


    h_0(n) = h_hat(1);
    h_1(n) = h_hat(2);

    Kn(1,n) = K(1);
    Kn(2,n) = K(2);

    Mn(1,n) = M(1,1);
    Mn(2,n) = M(2,2);
end

figure(1)
subplot(211)
plot(0:N-2,H(1,:))
xlabel('Sample number, n')
ylabel('Tap weight, h_n[0]')
grid on
axis([0 100 0 2])

subplot(212)

plot(0:N-2,H(2,:))
xlabel('Sample number, n')
ylabel('Tap weight, h_n[1]')
grid on
axis([ 0 100 0 2])

figure(2)
subplot(311)
plot(0:N-2,V(1:end-1))
xlabel('Sample number, n')
ylabel('Channel input, v[n]')
grid on
axis([0 100 -1 3])

subplot(312)
plot(0:N-2,y)
xlabel('Sample number, n')
ylabel('Noiseless channel output, y[n]')
grid on
axis([0 100 -1 3])

subplot(313)
plot(0:N-2,x1)
xlabel('Sample number, n')
ylabel('Channel output, x[n]')
grid on
axis([0 100 -1 3])   

figure(3)
subplot(211)
plot(0:N-2,H(1,:))
hold on 
plot(0:N-2,h_0)
xlabel('Sample number, n')
ylabel('True and estimated tap weight, h_n[0]')
grid on
axis([0 100 0 2])
legend('True', 'Estimate')

subplot(212)
plot(0:N-2,H(2,:))
hold on 
plot(0:N-2,h_1)
xlabel('Sample number, n')
ylabel('True and estimated tap weight, h_n[1]')
grid on
axis([0 100 0 2])
legend('True', 'Estimate')

figure(4)
subplot(211)
plot(0:N-2,Kn(1,:))
xlabel('Sample number, n')
ylabel('Kalamn gain, K_1[n]')
grid on
axis([0 100 -0.6 1])

subplot(212)
plot(0:N-2,Kn(2,:))
xlabel('Sample number, n')
ylabel('Kalamn gain, K_2[n]')
grid on
axis([0 100 -0.6 1])

figure(5)
subplot(211)
plot(0:N-2,Mn(1,:))
xlabel('Sample number, n')
ylabel('Minimum MSE, M_{11}[n]')
grid on
axis([0 100 0 0.2])

subplot(212)
plot(0:N-2,Mn(2,:))
xlabel('Sample number, n')
ylabel('Minimum MSE, M_{11}[n]')
grid on
axis([0 100 0 0.2])
