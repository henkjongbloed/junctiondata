function M = assembleM(t, dr, Ntide, Nspace)

N = length(t);

Ns = [1, 4, 10];
Nt = [1, 3, 5];
Np = 3*Nt(Ntide+1)*Ns(Nspace+1);

M = zeros(N, Np/3); %model matrix

TM2 = 12.42;

if Ntide == 0
    c = @(t) 1;
elseif Ntide ==1
    c = @(t) [1 cos(2*pi/TM2*t) sin(2*pi/TM2*t)];
else
    c = @(t) [1 cos(2*pi/TM2*t) sin(2*pi/TM2*t) cos(4*pi/TM2*t) sin(4*pi/TM2*t)];
end

% elseif lp == 4 %  first order approximation
%     r = [ones([lr,1]), dr'];
% elseif lp == 10  %  second order approximation
%     r = [ones([lr,1]), dr', dr(1,:)'.^2/2, dr(1,:)'.*dr(2,:)', dr(1,:)'.*dr(3,:)', dr(2,:)'.^2/2, dr(2,:)'.*dr(3,:)', dr(3,:)'.^2];
% end

if Nspace == 0
    for i = 1:N
        M(i,:) = [c(t(i))];
    end
end

if Nspace == 1
    for i = 1:N
        M(i,:) = [c(t(i)), dr(1,i).*c(t(i)), dr(2,i).*c(t(i)), dr(3,i).*c(t(i))];
    end
end

if Nspace == 2
    for i = 1:N
        M(i,:) = [c(t(i)), dr(1,i).*c(t(i)), dr(2,i).*c(t(i)), dr(3,i).*c(t(i)), dr(1,i)'.^2/2.*c(t(i)), dr(1,i)'.*dr(2,i)'.*c(t(i)), dr(1,i)'.*dr(3,i)'.*c(t(i)), dr(2,i)'.^2/2.*c(t(i)), dr(2,i)'.*dr(3,i)'.*c(t(i)), dr(3,i)'.^2/2.*c(t(i))];
    end
end

% M3 = zeros(N);
% M3(1:size)

end