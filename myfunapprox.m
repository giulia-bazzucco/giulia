function f_approx = myfunapprox(A,x)
s = size(A);
t = cell(length(s),1);
for i = 1:length(s)
    t{i} = chebyshevT(0:s(i)-1,x(i));
end
Y = tensor(A, s);
% for i = 1:s(3)
%     Y(:,:,i)=A(:,:,i);
% end
array = 1:length(s);
f_approx = ttm(Y,t,array);
end