function [amp, phas] = getAmpPhase(X, t)
% X is N by K matrix, containing temporally varying columns and K different
% datasets.
% t must be a column vector of times at which samples were taken

% par is a 4-by-K matrix containing columns [A2; A4; Phi2; Phi4].

T2 = 12.42; O2 = 2*pi/T2;
T4 = T2/2;  O4 = 2*pi/T4;

%We now employ the normal equations. Model:
% u = A2*cos(O2*t + P2) + A4*cos(O4*t + P4);

A = [cos(O2*t) sin(O2*t) cos(O4*t) sin(O4*t)]; %M2 and M4 tide only.

p = (A'*A)\A'*X;

amp(1,:) = sqrt(p(1,:).^2+p(2,:).^2); % A2
amp(2,:) = sqrt(p(3,:).^2+p(4,:).^2); % A4

phas(1,:) = atan2(-p(2,:), p(1,:)); % Phi2
phas(2,:) = atan2(-p(4,:), p(3,:)); % Phi4

end


