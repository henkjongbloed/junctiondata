% Fourdimensional regression: Estimate (u(t,x,y,z), v(t,x,y,z), w(t,x,y,z))
% under constraint of mass conservation using basis functions

%% Generate artificial data

% r = [t, x, y, z] Coordinates of measured velocities, t in hours, rest in
% meters
n = 2;
t = linspace(0, 13.5, n);
r = [t; 2*rand([1,n])-1; sin(2*pi*t); rand([1,n])-1];

%% Forward model
p = 1:6;
u = fmu(r,p);
s = fms(r);

%% Simulate ADCP measurements (only measure radial component of velocity)
rb = [t; 0*rand([1,n]); sin(2*pi*t)+0.01*randn(size(t)); 0*rand([1,n])];
Q = zeros([3,n]);
b = zeros([1,n]);
sig = 0*.01;
for i = 1:n
    Q(:,i) = (r(2:end,i)-rb(2:end,i))./norm(r(2:end,i)-rb(2:end,i));
    b(i) = Q(:,i)'*u(:,i) + sig*randn;
end

%% Inverse model

% Caveat: We know the forward model (in reality not the case)

pest = Q'\b';
er = norm(pest-p);

%% Investigate performance

