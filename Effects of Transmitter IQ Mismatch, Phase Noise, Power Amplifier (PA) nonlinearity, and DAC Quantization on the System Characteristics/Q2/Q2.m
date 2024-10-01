j = sqrt(-1);
num_syms = 1e3;


%%% generate symbols for UE-1 and UE-2
syms = sqrt(2)*(randi([0,1],2, num_syms)-0.5 + 1j*((randi([0,1], 2, num_syms)-0.5)));

syms(2,:) = 0*syms(2,:);



%%% upsmaple and pulse shape (use this only for to see impact of PA and DAC
%%% non_linearity on frequ. domains)

%upsample
syms_up = transpose(upsample(transpose(syms),4));

%pulse shaping filter
g = rcosdesign(0.5, 8, 4, 'normal');
%   g = g/max(g);
syms_tilde_1 = conv(syms_up(1,:),g);
syms_tilde_2 = conv(syms_up(2,:),g);

syms_tilde = [syms_tilde_1; syms_tilde_2];

%number of BS antennas
M8 = 8;  
phi_1 = 30; %angle to UE-1 in degrees
h1_M8 = exp(-1j*pi*sind(phi_1)*(0:M8-1));  %channel vector to UE-1

phi_2 = 40; %angle to US-2 in degrees 
h2_M8 = exp(-1j*pi*sind(phi_2)*(0:M8-1));

%Number of BS antennas
M32 = 32;   
h1_M32 = exp(-1*j*pi*sind(phi_1)*(0:M32-1));

h2_M32 = exp(-1j*pi*sind(phi_2)*(0:M32-1));





%%% ZF pre-coding

H_M8 = [h1_M8; h2_M8];
H_M32 = [h1_M32; h2_M32];

P8 = H_M8'*inv(H_M8*H_M8');
P32 = H_M32'*inv(H_M32*H_M32'); 

x8 = P8*syms_tilde;
x32 = P32*syms_tilde;

%%% DAC uniform quantization 

%%%%%%%% b=12, M=8
b12 = 12;
delta_M8_b12 = 2*(1/M8)/2^b12;  %delta, separation distance between boundary
I_M8_b12 = -1/M8:delta_M8_b12:1/M8;  %Interval: [-1/M,1/M]

code_M8_b12 = zeros(1,length(I_M8_b12)-1);
%code_M32_b12 is simply a vector that includes the quantization levels
for i=1:length(code_M8_b12)
    code_M8_b12(i) = I_M8_b12(i)+delta_M8_b12/2;
end
x8_b12_real = real(x8);
x8_b12_imag = imag(x8);

y8_b12_real = qunt(x8_b12_real,code_M8_b12); %quantize the real part
y8_b12_imag = qunt(x8_b12_imag,code_M8_b12); %quantize the imaginary part
y8_b12 = y8_b12_real + j*y8_b12_imag;        %output of quantizer

%%%%%%%%% b=4, M=8
b4=4;
delta_M8_b4 = 2*(1/M8)/2^b4;
I_M8_b4 = -1/M8:delta_M8_b4:1/M8;
code_M8_b4 = zeros(1,length(I_M8_b4)-1);
%code_M32_b12 is simply a vector that includes the quantization levels
for i=1:length(code_M8_b4)
    code_M8_b4(i) = I_M8_b4(i) + delta_M8_b4/2;
end
x8_b4_real = real(x8);
x8_b4_imag = imag(x8);

y8_b4_real = qunt(x8_b4_real,code_M8_b4);  %quantize the real part
y8_b4_imag =  qunt(x8_b4_imag,code_M8_b4); %quantize the imaginary part
y8_b4 = y8_b4_real + j*y8_b4_imag;

%%%%%%% b=12, M=32
delta_M32_b12 = 2*(1/M32)/2^b12;
I_M32_b12 = -1/M32:delta_M32_b12:1/M32;
code_M32_b12 = zeros(1,length(I_M32_b12)-1);
%code_M32_b12 is simply a vector that includes the quantization levels
for i=1:length(code_M32_b12)
    code_M32_b12(i) = I_M32_b12(i) + delta_M32_b12/2;
end
x32_b12_real = real(x32);
x32_b12_imag = imag(x32);

y32_b12_real = qunt(x32_b12_real,code_M32_b12);
y32_b12_imag = qunt(x32_b12_imag,code_M32_b12);
y32_b12 = y32_b12_real + j*y32_b12_imag;

%%%%%% b=4, M=32
delta_M32_b4 = 2*(1/M32)/2^b4;
I_M32_b4 = -1/M32:delta_M32_b4:1/M32;
code_M32_b4 = zeros(1,length(I_M32_b4)-1); 
%code_M32_b12 is simply a vector that includes the quantization levels
for i=1:length(code_M32_b4)
    code_M32_b4(i) = I_M32_b4(i) + delta_M32_b4/2;
end
x32_b4_real = real(x32);
x32_b4_imag = imag(x32);

y32_b4_real = qunt(x32_b4_real,code_M32_b4);
y32_b4_imag = qunt(x32_b4_imag,code_M32_b4);
y32_b4 = y32_b4_real + j*y32_b4_imag;





%%%%Q2.A
%%%(a)   M=8, b=12,  Ideal PA
z8_b12 = y8_b12;
figure(1); 
pwelch(z8_b12(1,:))
title('PSD, M=8, b=12')

%%%(b)   M=8, b=4,   Ideal PA
z8_b4 = y8_b4;
figure(2); pwelch(z8_b4(1,:))
title('PSD, M=8, b=4')

%%%(c)  M=32, b=4,   Ideal PA
z32_b12 = y32_b12;
figure(3); pwelch(z32_b12(1,:))
title('PSD, M=32, b=12')

%%%(d)  M=32, b=12,  Ideal PA
z32_b4 = y32_b4;
figure(4); pwelch(z32_b4(1,:))
title('PSD, M=32, b=4')


%%% Q2(B)
%%% angular response of UE-1 signal with perfect transmission (no DAC
%%% quantization and PA non-linearity)

phi_range = -90:90;
ang_response8_b12 = zeros(1, length(phi_range));

%%% (a)  M=8, b=12
for phi = phi_range
    a8 = exp(-1*j*pi*sind(phi)*(0:M8-1));
    ang_response8_b12(phi==phi_range) = rms(a8*z8_b12).^2;
end
figure(5)
plot(phi_range, 10*log10(ang_response8_b12),'LineWidth',2);
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB), M=8, b=12');



%%% (b)  M=8, b=4
ang_response8_b4 = zeros(1, length(phi_range));
for phi = phi_range
    a8 = exp(-1*j*pi*sind(phi)*(0:M8-1));
    ang_response8_b4(phi==phi_range) = rms(a8*z8_b4).^2;
end
figure(6)
plot(phi_range, 10*log10(ang_response8_b4),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=8, b=4');

%%%(c)  M=32, b=12
ang_response32_b12 = zeros(1, length(phi_range));
for phi=phi_range
    a32 = exp(-1*j*pi*sind(phi)*(0:M32-1));
    ang_response32_b12(phi==phi_range) = rms(a32*z32_b12).^2;
end
figure(7)
plot(phi_range, 10*log10(ang_response32_b12),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=32, b=12');

%%% (d)  M=32, b=4
ang_response32_b4 = zeros(1, length(phi_range));
for phi=phi_range
    a32 = exp(-1*j*pi*sind(phi)*(0:M32-1));
    ang_response32_b4(phi==phi_range) = rms(a32*z32_b4).^2;
end
figure(8)
plot(phi_range, 10*log10(ang_response32_b4),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=32, b=4');


%%% Q2(C)
%%%(a)  M=8, ideal PA
z8_a = x8;
figure(9)
pwelch(z8_a(1,:));
title('PSD, M=8, ideal PA')

%%% (b)   M=8, non ideal PA
z8_b = x8 + (-133)*x8.*(abs(x8)).^2;
figure(10)
pwelch(z8_b(1,:));
title('PSD, M=8, non ideal PA')

%%% (c) M=32, ideal PA
z32_c = x32;
figure(11)
pwelch(z32_c(1,:));
title('PSD, M=32, ideal PA')

%%% (d) M=32, non ideal PA
z32_d = x32 + (-133)*x32.*(abs(x32)).^2;
figure(12)
pwelch(z32_d(1,:));
title('PSD, M=32, non ideal PA')

%%% Q2 (D)
%%% (a) M=8, ideal PA
ang_response8_a = zeros(1, length(phi_range));
for phi=phi_range
    a8 = exp(-1*j*pi*sind(phi)*(0:M8-1));
    ang_response8_a(phi==phi_range) = rms((a8*z8_a).^2);
end
figure(13)
plot(phi_range, 10*log10(ang_response8_a),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=8, ideal PA');

%%% (b) M=8, non ideal PA
ang_response8_b = zeros(1, length(phi_range));
for phi=phi_range
    a8 = exp(-1*j*pi*sind(phi)*(0:M8-1));
    ang_response8_b(phi==phi_range) = rms(a8*z8_b).^2;
end
figure(14)
plot(phi_range, 10*log10(ang_response8_b),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=8, non ideal PA');

%%% (c), M=32, ideal PA
ang_response32_c = zeros(1, length(phi_range));
for phi=phi_range
    a32 = exp(-1*j*pi*sind(phi)*(0:M32-1));
    ang_response32_c(phi==phi_range) = rms(a32*z32_c).^2;
end
figure(15)
plot(phi_range, 10*log10(ang_response32_c),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=32, ideal PA');

%%% (d), M=32, non ideal PA
ang_response32_d = zeros(1, length(phi_range));
for phi=phi_range
    a32 = exp(-1*j*pi*sind(phi)*(0:M32-1));
    ang_response32_d(phi==phi_range) = rms(a32*z32_d).^2;
end
figure(16)
plot(phi_range, 10*log10(ang_response32_d),'LineWidth',2)
grid on
ylabel('g(\phi) (dB)');
xlabel('\phi (degree)');
title('g(\phi) (dB),M=32, non ideal PA');




