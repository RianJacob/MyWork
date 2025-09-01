function x_hat = QAM_64_demod(x)
j = sqrt(-1);

a = sign_fun(real(x));
b = sign_fun(imag(x));

x_real = abs(real(x));
x_imag = abs(imag(x));

if x_real<=2
    x_hat_real = 1;
elseif x_real <= 4
    x_hat_real = 3;
elseif x_real <= 6
    x_hat_real = 5;
else
    x_hat_real = 7;
end


if x_imag<=2
    x_hat_imag = 1;
elseif x_imag <= 4
    x_hat_imag = 3;
elseif x_imag <= 6
    x_hat_imag = 5;
else
    x_hat_imag = 7;
end

x_hat = a*x_hat_real + j*b*x_hat_imag;

end