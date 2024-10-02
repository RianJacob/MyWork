clear all
close all

N = 64;
M = 8;
L = 8;
n = 16;  % OFDM Symbols
j = sqrt(-1);

F = fft(eye(N))/sqrt(N);

p_ind = 4:8:60;
p_d = setdiff(1:N,p_ind);

SNR = [0:2:24,25,26];

n_c = 100;  % 100 iterations for channel

map = 1/sqrt(2)*[-1 1];

error_v = zeros(1,length(SNR));
error_va = zeros(1,length(SNR));
error_vb = zeros(1,length(SNR));
for ch=1:n_c
    h = (randn(1,M) + j*randn(1,M))/sqrt(2*M);
    h_app = [h zeros(1,N-M)];
    H = zeros(N,N);
    H(1,:) = h_app;
    
    for k=2:N
        H(k,:) = circshift(h_app,k-1);
    end
    
    
    for i = 1:length(SNR)
        n_error = 0;
        n_error_a = 0;
        n_error_b = 0;
        x = map(randi([1 2],N,n)) + j*map(randi([1 2],N,n));
        x_bar = F'*x;
        
        if (ch==1 && SNR(i)==25)
            x_25 = x;  % transmitted Sequence
            
        end
            
        
        xa = map(randi([1 2],N,n)) + j*map(randi([1 2],N,n));
        xa(:,1) = ones(N,1) * (1+j);   % used for training 
        xa_bar = F'*xa;
        
        xb = zeros(N,n);  % Scheme b
        xb(p_ind,:) = ones(length(p_ind),n)*(1+j);  % used for training
        
        xb(p_d,:) = map(randi([1 2],length(p_d),n)) + j*map(randi([1 2],length(p_d),n));
        
        xb_bar = F'*xb;
        
        
        var = 10^(-SNR(i)/10);
        v_bar = sqrt(var/2)*(randn(N,n) + j*randn(N,n));
    
        y_bar = H*x_bar+ v_bar;
        y = F*y_bar;       % DFT
        
        ya_bar = H*xa_bar + v_bar;
        ya = F*ya_bar;
        
        if(ch==1 && SNR(i)==25)
            ya_25 = ya;
        end
        
        yb_bar = H*xb_bar+v_bar;
        yb = F*yb_bar;
    
        Lamda = diag(F'*h_app(:));   % assuming exact knowledge of the channel
        x_hat = inv(Lamda)*y/sqrt(N);

        if (ch==1 && SNR(i)==25)
            x_hat25 = x_hat;
        end
        

        
        

        x_hat = slicer_4(x_hat);
       
        
        %Scheme a
        Qa = F(1:M,1:N);
        h_hata = 1/sqrt(N)*inv(Qa*Qa')*Qa*(conj(xa(:,1)).*ya(:,1));
        Lamda_a = diag(F'*[h_hata; zeros(N-L,1)]);
        x_hata = inv(Lamda_a)*ya/sqrt(N);
        if (ch==1 && SNR(i)==25)
            x_hata_25 = x_hata;
        end
        
        x_hata = slicer_4(x_hata);
        
        % Scheme b
        Qb = F(1:M,p_ind);
        h_hatb = 1/sqrt(N)*inv(Qb*Qb')*Qb*(conj(xb(p_ind,:)).*yb(p_ind,:));
        h_b = zeros(N,1);
        for index=1:M
            h_b(index) = mean(h_hatb(index,:));
        end
        Lamda_b = diag(F'*h_b);
        x_hatb = inv(Lamda_b)*yb/sqrt(N);
        x_hatb = slicer_4(x_hatb);
        
        for m=1:N
            for k=2:n-1     % ignore first OFDM symbol to make fair comparison with other schemes (same number of bits)
                if (real(x_hat(m,k))~= real(x(m,k)))
                    n_error = n_error+1;
                end
                if(imag(x_hat(m,k)) ~= imag(x(m,k)))
                    n_error = n_error+1;
                end
            end
        end
        
        for ar=1:N
            for ac=2:n-1
                if (real(x_hata(ar,ac)) ~= real(xa(ar,ac)))
                    n_error_a = n_error_a + 1;
                end
                if (imag(x_hata(ar,ac)) ~=imag(xa(ar,ac)))
                    n_error_a = n_error_a + 1;
                end
            end
        end
        
        for br=p_d
            for bc=1:n
                if (real(x_hatb(br,bc)) ~= real(xb(br,bc)))
                    n_error_b = n_error_b + 1;
                end
                if (imag(x_hatb(br,bc)) ~= imag(xb(br,bc)))
                    n_error_b = n_error_b + 1;
                end
            end
        end
        
        error_v(i) = error_v(i)+n_error;
        error_va(i) = error_va(i)+n_error_a;
        error_vb(i) = error_vb(i)+n_error_b;
                
        
    end
    
end
Ber = error_v/(2*N*(n-2)*n_c);  
Bera = error_va/(2*N*(n-2)*n_c);
Berb = error_vb/(2*length(p_d)*n*n_c);
    



semilogy(SNR(1:11), Ber(1:11),'-^')
hold on 
semilogy(SNR(1:11), Bera(1:11),'-o')
semilogy(SNR(1:11), Berb(1:11),'-s')


grid on 
xlabel('SNR (dB)')
ylabel('BER')
title('BER vs SNR for two training schemes')
legend('Exact', 'Scheme a', 'Scheme b')

figure(2)
subplot(221)
plot(x_25,'b+')
axis([-1,1,-1,1]*2)
title('transmitted signal')
grid on
subplot(222)
plot(ya_25,'b+')
axis([-1,1,-1,1]*3)
title('Received signal')
grid on
subplot(223)
plot(x_hat25,'b+')
axis([-1,1,-1,1]*2)
title('Estimated Sequence using exact channel')
grid on
subplot(224)
plot(x_hata_25,'b+')
axis([-1,1,-1,1]*2)
title('Estimated Sequence using Scheme (a)')
grid on

    
    
    
    
    
    
    
    
    
    

    
    
    