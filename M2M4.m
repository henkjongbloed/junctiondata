function f = M2M4(A, P, t)

T2 = 12.42; O2 = 2*pi/T2;
T4 = T2/2;  O4 = 2*pi/T4;

f = A(1).*cos(O2.*t + P(1)) + A(2)*cos(O4.*t + P(2));

end