function x_hat = QAM_16_demod(x)
j = sqrt(-1);

a = sign_fun(real(x));
b = sign_fun(imag(x));

x_real = abs(real(x));
x_imag = abs(imag(x));

if x_real<=2
    x_hat_real = 1;
else
    x_hat_real = 3;
end

if x_imag<=2
    x_hat_imag = 1;
else
    x_hat_imag = 3;
end

x_hat = a*x_hat_real + j*b*x_hat_imag;

end