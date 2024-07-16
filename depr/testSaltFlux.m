% Test Salt Flux Decompositions

% Dim1 = y, Dim2 = z, Dim3 = t

% Create artifical data

Bh = 100;
H = 20;

y = -Bh:1:Bh;
z = 0:-1:-H;

[Y,Z] = meshgrid(y,z);

sY = 50; mY = 0;
sZ = 10;  mZ = 0;  

Ucross = twoDNormal(Y,Z,mY,mZ,sY,sZ);
Scross = twoDNormal(Y,Z,mY,mZ-H,3*sY,sZ);

h = 0:12;
phi = pi/2;
tidU = sin(2*pi/12*h + 0);
tidS = sin(2*pi/12*h + phi);
A = zeros(1,1,13);
A(:,:,:) = Bh*(H+tidS);

Q = 1000;
Ulat = twoDNormal(Y,Z,Bh/3,-H/3,sY,sZ);
Sver = twoDNormal(Y,Z,mY,mZ-H,3*sY,sZ);
U = zeros(size(Ucross,1),size(Ucross,2),length(h));
S = zeros(size(Scross,1),size(Scross,2),length(h));
for hr = 1:length(h)
    U(:,:,hr) = Q/(A(hr)) + Ulat + tidU(hr).*Ucross;
    S(:,:,hr) = 10*(1+ Sver + tidS(hr).*Scross);
end

figure(1)
hr = 1;
subplot(3,1,1)
caxis([min(U,[],'all'), max(U,[],'all')])
pU = pcolor(Y,Z,U(:,:,hr));
colorbar;
subplot(3,1,2)
caxis([min(S,[],'all'), max(S,[],'all')])
pcolor(Y,Z,S(:,:,hr));
colorbar;
%pause(1)
for hr = 2:length(h)
    pause(1)
    subplot(3,1,1)
    caxis([min(U,[],'all'), max(U,[],'all')])
    pcolor(Y,Z,U(:,:,hr));
    colorbar;
    subplot(3,1,2)
    caxis([min(S,[],'all'), max(S,[],'all')])
    pcolor(Y,Z,S(:,:,hr));
    colorbar;
    pause(1)
end
%S = U;


[C,T] = crosssectionalAvg(U,S,A);

SF = log10(cell2mat(T));
SF(SF<0) = 0;

subplot(3,1,3)
bar(SF)
title('Kranenburg 2016 (log)')
