% GeometryTesting.m
% Get intuition about coordinate transforms in the ADCP toolbox.

% Ship sails over a line in (x,y,z) space

width = 100;
n = 1000;
p_a = linspace(0,width,n)'.*[2, -3, 0] + 0*[randn(n,2), zeros(n,1)];
p_abar = mean(p_a,1)';
[V,D] = eig(p_a'*p_a);

[~,mi] = find(D==max(max(D)));

t = V(:,mi)/norm(V(:,mi));
k = [0;0;-1];
o = [1; - t(1)/t(2)*1; 0]; o = o/norm(o);

dot(t,k)
dot(t,o)
dot(k,o)

C = [o';t';k'];
C0 = [o';t';0 0 0];

p = [-100; -150; -8];
m = C*p - C0*p_abar; %Formula (9) of Vermeulen 2014
%plot(p_a)

