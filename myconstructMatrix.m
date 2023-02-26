function M = myconstructMatrix (C,k,n,f)

%T2 = reshape(permute(T2,[2,1,3]),n(2),r(1)*r(3));

Cp = C;
cheb = @(i,n) -cos((i-1).*pi/(n-1));
Cp{k} = cheb(1:n(k), n(k)); %all of the cheb points

N = length(C);

T = myconstructTensor (f,Cp);
T = permute(T , [k, 1:k-1, k+1:N]); % I have to permute the columns like above
M = reshape(T, n(k), []); %flattening

% for i = 1 : prod(length(C{i}))
%     ind = sub2ind(i, length(C{i}));
%     
%     for j = 1 : n(i)
%         [linind,] = ind2sub(j, ind);
%     end
% end
% 
% T{linind} = M(j, );

end
