%% Adaptive Cross Approximation (ACA)
% see Algorithm 1 of the Functional Tucker approx using Chebyschev
% interpolation

% iterative application of an approximation of Gaussian elimination with complete pivoting

function[Mc, Mr, Mt, rInd, cInd] = myACA (M, tol, maxIter)
Mc = [];
Mr = [];
Mt = [];
rInd = [];
cInd = [];
Moriginal = M;

for iter=1:maxIter
    
    [err,I2] = max(abs(M(:)));
    if isempty(err) || err < tol
        Mc = Moriginal(:,cInd);
        Mr = Moriginal(rInd,:)';
        Mt = Moriginal(rInd,cInd);
        return
    end
    
     [I,J] = ind2sub(size(M), I2); %ind2sub -> convert linear indices to subscripts
    rInd = [rInd, I];
%     cInd = [cInd, J];
    if J < size(M,2)
      cInd = [cInd, J];
    end

    M = M-M(:,J)*M(I,:)./M(I,J);
end
Mc = Moriginal(:,cInd);
Mr = Moriginal(rInd,:)';
Mt = Moriginal(rInd,cInd);
end





% yields a low-rank approximation in terms of function fibers and interpolating these fibers
    
    