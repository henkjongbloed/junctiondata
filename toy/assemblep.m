function p = assemblep(t,s)


% if nargin==1
%     p = varargin{1};
% elseif nargin==2
%     p = zeros([12,1]);
%     p([1, 5, 9]) = varargin{1};
%     p([2:4, 6:8, 10:12]) = varargin{2};
% elseif nargin==3
%     p = zeros([30,1]);
%     p([1, 11, 21]) = varargin{1};
%     p([2:4, 12:14, 22:24]) = varargin{2};
%     p([5:10, 15:20, 25:30]) = varargin{3};
% end


Ns = [1, 4, 10];
Nt = [1, 3, 5];

p = rand([3*Ns(s+1)*Nt(t+1),1])-1/2;
p = ones([3*Ns(s+1)*Nt(t+1),1]);
p(1) = 1;


end