% UTest

% Specify a flow profile in 3D and compute corresponding decompositions

function [C, T] = crosssectionalAvg(U, S, CS)

% Input: U(x = x_0, y, z, t) 3D double on rectangular underlying grid.
% Output: Kranenburg 2016

%U = Ua + Ub + Uc + Ud

d = size(U);

%[Ubar, Uhat, Ubarhat] = avg(U);

Ubar = bar(U);
Uhat = hat(U);
Ubarhat = bar(hat(U));

disp(adiff(bar(Ubar),Ubar))

Ua = Ubarhat; %cross sectionally and tidally averaged
Ub = Ubar - Ua; %cross sectionally averaged, tidally varying
Uc = Uhat - Ua; %cross sectionally varying, tidally averaged
Ud = U - Ua - Ub - Uc; %cross sectionally varying, tidally varying

%[Sbar, Shat, Sbarhat] = avg(S);
Sbar = bar(S);
Shat = hat(S);
Sbarhat = bar(hat(S));
Sa = Sbarhat;
Sb = Sbar - Sa;
Sc = Shat - Sa;
Sd = S - Sa - Sb - Sc;

A = hat(CS);
a = CS - A;

T{1} = A.*Ua.*Sa;
T{2} = Sa.*hat(Ub.*a);
T{3} = Sa.*hat(bar(Uc.*(A + a)));
T{4} = Sa.*hat(bar(Ud.*(A + a)));
T{5} = Ua.*hat(Sb.*a);
T{6} = hat(Ub.*Sb.*a);
T{7} = hat(Sb.*bar(Uc).*a);
T{8} = hat(Sb.*bar(Ud).*(A+a));
T{9} = Ua.*hat(bar(Sc.*(A+a)));
T{10} = hat(Ub.*bar(Sc).*a);
T{11} = A.*hat(Uc.*Sc);
T{12} = bar(Sc.*hat(Ud.*(A+a)));
T{13} = Ua.*(hat(bar(Sd.*(A+a))));
T{14} = hat(bar(Sd).*Ub.*(A+a));
T{15} = bar(Uc.*hat(Sd.*(A+a)));
T{16} = hat(bar(Sd.*Ud.*(A+a)));

for i=1:16
    T{i} = mean(T{i},'all');
end

C.Ua = Ua; C.Ub = Ub; C.Uc = Uc; C.Ud = Ud; 
C.Sa = Sa; C.Sb = Sb; C.Sc = Sc; C.Sd = Sd; 
C.A = A; C.a = a;

end

function hat = hat(x)
d = size(x);
if sum(d(1:2)) == 2 %x is already cross sectionally averaged
    hat = repmat(mean(x, 3), [1,1,d(3)]);
else
    hat = repmat(mean(x, 3), [1,1,d(3)]);
end
end

function bar = bar(x)
d = size(x);
if sum(d(1:2)) == 2 %x is already cross sectionally averaged
    bar = repmat(mean(x, [1,2]), [d(1), d(2), 1]);
else
    bar = repmat(mean(x, [1,2]), [d(1), d(2), 1]);
end
end


function adiff = adiff(x,y)
adiff = max(abs(x - y),[],'all');
end