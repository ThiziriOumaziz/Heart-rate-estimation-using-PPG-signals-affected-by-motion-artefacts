function z = Smoother(y,lambda) 
    y= transpose(y);
    m = length(y);
    I = eye(m);
    D = diff(I);
    z = (I + lambda * (D') * D)\y;
end