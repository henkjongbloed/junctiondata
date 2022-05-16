% TestGeom


B = 100;
H = 25;
n = 10000;
x = rand([n,2])*diag([B,H]);



plot(zb(x(:,1), B, H))


function zb = zb(y, B, H)
    zb = nan(size(y));
    y = sort(y);
    zb(y<=B/2) = -H/B.*y(y<=B/2) + H/2;
    zb(y>B/2) = 4*H/(3*B^2).*(y(y>B/2).^2 - B.*y(y>B/2)) + H/3;
end

function eta = eta(t,a)
    eta = a.*sin(t);
end