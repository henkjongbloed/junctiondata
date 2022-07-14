function Mb = assembleMb(M, r, rb)


Mb = zeros(size(M,1), 3*size(M,2));
Q = zeros(size(r));
for i = 1:size(r,2)
    Q(:,i) = (r(:,i)-rb(:,i))./norm(r(:,i)-rb(:,i));
    Mb(i,:) = [Q(1,i)*M(i,:), Q(2,i)*M(i,:), Q(3,i)*M(i,:)];
end