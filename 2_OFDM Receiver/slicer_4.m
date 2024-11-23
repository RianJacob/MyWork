function [y] = slicer_4(u)
a = sign(real(u));
b = sign(imag(u));

j = sqrt(-1);
y = 1/sqrt(2)*(a+j*b);
end

