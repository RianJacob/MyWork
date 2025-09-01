function y = QAM16_demod(x)
j = sqrt(-1);
y = zeros(1,length(x));

for i=1:length(x)

    a = real(x(i));
    b = imag(x(i));

    sign_a = sign(a);
    sign_b = sign(b);

    if (abs(a)<=2)
        a=1;
    else
        a=3;
    end


    if (abs(b)<=2)
        b=1;
    else
        b=3;
    end


    y(i) = a*sign_a + j*b*sign_b;
end

end

