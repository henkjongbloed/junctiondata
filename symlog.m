function y = symlog(x,c)
if nargin==1
    c = 0;
end
y = sign(x).*(log10(1+abs(x)./(10.^c))) ;
end