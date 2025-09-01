function x_hat = QAM_4_demod(x)
j = sqrt(-1);

a = sign_fun(real(x));
b = sign_fun(imag(x));

x_hat = a + j*b;

end