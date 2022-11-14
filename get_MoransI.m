function [I_loc, I] = get_MoransI(U, W, varargin)

if varargin{1} == 1
    p = U.pars;
    %p = ones(size(U.pars));
elseif varargin{1}==0
    p = U.pars0;
elseif varargin{1}==-1
    p = U.pars - U.pars0;
end


m2 = sum((p - mean(p, 1)).^2, 1)./size(p, 1);
m22 = (p - mean(p, 1))./m2;

I_loc = m22.*(W*(p-mean(p, 1)));%./sum(dat.IM,"all");
I = mean(I_loc,1);
end