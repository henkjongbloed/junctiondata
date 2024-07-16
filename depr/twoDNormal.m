function f = twoDNormal(Y,Z,mY,mZ,sY,sZ)


f = 1/(2*pi*sY*sZ)*exp(-1/2*(((Y-mY)/sY).^2 + ((Z-mZ)/sZ).^2));
f = f./max(max(f));

end