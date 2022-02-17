function t = symlog(x)
C=2;
t=sign(x).*log(1+abs(x)./10.^C);
end