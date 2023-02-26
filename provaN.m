function mycfN = provaN (f, N)
%% codice più pulito rispetto a ttryN.m

% Initialize

n = linspace(17,17,N); %coarseResolution
m = n;                   %fineResolution
r = linspace(6,6,N);             %rank
tol = 1e-6;
err_vec = Inf;
cheb = @(i,n) -cos((i-1).*pi/(n-1)); %punti di chebyshev
reffun = @(n) floor(sqrt(2)^(floor(2*log2(n)) + 1)) + 1;
pref = chebfunpref(); %chebfunpref = class for managing Chebfun construction-time preferences


%% main loop
% happy = 0;
% while ~happy
% while err_vec > tol
%% Phase 1
   happyPhase1 = 0;
   while ~happyPhase1 
       
       J = cell(N,1);
       for i=1:N
       J{i} = initializeIndexRandomly(r(i), n(i));
       end
       
       C = cell(N,1);
       for i=1:N
           C{i} = cheb(J{i}, n(i));
       end
       
        
        U = cell(N,1);
        I2 = cell(N,1);
        M = cell(N,1);
        for i=1:N
            Ci = C;
            Ci{i} = -cos(((1:n(i))-1).*pi/(n(i)-1));
            T = myconstructTensor(f,Ci);
            M{i} = myconstructMatrix(C,i,n,f);
            [U{i}, ~, ~, I,I2{i}] = myACA(M{i}, tol, n(i));  % I and I2?
            r(i) = size(I,2);
            J{i} = I;
            C{i} = cheb(J{i}, n(i));
        end
        
            
            refineFlag = 0;
            for i=1:N
                while r(i)*2*sqrt(2) > n(i)
                    n(i) = reffun(n(i)); %reffun --> see the initialization
                    refineFlag = 1;
                end
            end
            if refineFlag == 0
                happyPhase1 = 1;
                break
            end
       
   end
%    max_index = max(cellfun(@length,I2));
%    disp(max_index);
   %% Phase 2

for i=1:N
if size(J{i},i) == 0 
        mycfN.U = chebfun(zeros([n(i),1]), pref);
        mycfN.V = chebfun(zeros([n(2),1]), pref);
        mycfN.W = chebfun(zeros([n(3),1]), pref);
        mycfN.C = 0;

    else
m = n;

resolved = zeros(N,1); %variabili booleane che indicano se la funzione è stata risolta per ogni variabile
Uf = cell(N,1);
        for l = 1:N
            Uf{l} = U{l};
            ct2 = mycreateCT2(Uf{l});
            resolved(i) = happinessCheck(ct2, [], ct2.coeffs, [], pref);
            if ~resolved(l)
                m(l) = 2*m(l)-1;
            end
        end
        
        % Add function evaluations and check again
        while any(~resolved) 
            
           for p=1:N
               if resolved(p)
                    J{p} = 2*J{p}-1;
               end            
           end
% 


%             M = cell(N,1);
            
            for q=1:N
                if ~resolved(q)
%                     II{q}=ind2sub([size(J{q},2)],I2{q});
                    M{q} = myconstructMatrix(C,q,m,f);
%                     [U{q}, ~, ~, ~,I2{q}] = myACA(M{q}, tol, m(q));
                    Uf{q} = M{i}(:,I2{q});
%                     try
%                         U{q} = M{i}(:,I2{q});
%                     catch
%                         % gestione dell'errore, ad esempio stampare un messaggio di avviso
%                         disp('Attenzione: errore nell''assegnazione di U{q}.')
%                     end
%                    UU = cellfun(@(idx) cat(2, M(:, idx{:})), I2, 'UniformOutput', false);
%                     for d = 1:length(I2)
%                         UU{d} = M(:,I2{d});
                end
            end
% 
% %                     for d = 1:length(transpose(I2))
% %                         if numel(I2{d}) > 1
% %                             U{d} = M(:,I2{d});
% %                         else
% %                             U{d} = M(:,I2{d}(1));
% %                         end
% %                     end
%                    
%                 end
%             end
%                for q=1:N
%                    if ~resolved(q)
%                     M = myconstructMatrix(C,q,m,f);
%                     max_row = size(M, 1);
%                     for q = 1:length(I2)
%                         if any(I2{q} > max_row)
%                             error('Indice I2 fuori range per M');
%                         end
%                     end

%                     I2q = I2{q};
%                     if size(M,2) ~= length(I2q)
%                         error('Dimension mismatch between M and I2q');
%                     end
%                     U{q} = M(:,I2q);
%                     end
                
                resolved = zeros(N,1);
                for j = 1:N
                    Uf{j} = U{j};
                    ct2 = mycreateCT2(Uf{j});
                    resolved(j) = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                        if ~resolved(j)
                            m(j) = 2*m(j)-1;
                        end
                end
        end
   
        
            
            
        
    %% Phase 3
        
        % Compute factor matrices
        for k = 1:N
            [Q{k},~] = qr(Uf{k},0);
            K{k} = myDEIM(Q{k});    %I?
            U{k} = Q{k}/Q{k}(K{k},:);
        end
        
        % Store U(x)
         mycfN.U = cell(N,1);
            for i = 1:N
               mycfN.U{i} = chebfun(U{i}, pref);
            end
                
        % Compute the core
        C = cell(N,1);
            for i=1:N
                C{i} = cheb(K{i}, n(i));       
            end
        mycfN.C = myconstructTensor(f,C);
        
%         disp(size(mycfN.C))
       G = mycfN.C;
%        disp(G)
     
       for i = 1:N
            disp(size(mycfN.U{i}))
       end 
end
end
%% Calcolo dell'errore di approssimazione
coeffs = cell(N,1);
for i= 1:N 
    coeffs{i} = chebcoeffs(full(mycfN.U{i}));
end

% A = tensor(@zeros, [17 17 17]);
A = mycfN.C;

for i=1:N
A = tmprod(A, {coeffs{i}}, i);
end
% disp(A)
mycfN.A = A;

f_approx_vec = zeros(1, length(size(A)));
fval_vec = zeros(1, length(size(A)));
err_vec = zeros(1, length(size(A)));

for i = 1:10
    x = rand(1,length(size(A)));
    f_approx = myfunapprox(A,x);
    f_approx_vec(i) = f_approx;
    args = num2cell(x(1:length(size(A))));
    fval = f(args{:});%fval=f(x(1),x(2),x(3));
    fval_vec(i) = fval;
    err_vec(i) = norm(fval_vec(i) - f_approx_vec(i), "fro"); %abs(fval_vec - f_approx_vec)
end
% err_vec = norm(fval_vec - f_approx_vec, "fro"); %abs(fval_vec - f_approx_vec)
% Stampa i risultati
disp("f_approx_vec:");
disp(f_approx_vec);

disp("fval_vec:");
disp(fval_vec);

disp("err_vec:");
disp(err_vec);

% s = size(A);
% for i = 1:length(s)
%     x = ones(1,length(s)); %rand
%     t{i} = chebyshevT(0:s(i)-1,x(i));
% end
% disp(t)
% Y = tensor(A, s);
% % for i = 1:s(3)
% %     Y(:,:,i)=A(:,:,i);
% % end
% 
% fapprox = ttm(Y,t,[1,2,3]);
% mycfN.fapprox=fapprox;
% fval= f(1,1,1);
% err=norm(fval-fapprox, "fro");
% disp(err)

%   CALCOLO ERRORE
% if err_vec > tol
% %     if ~happy
%     
%     for i=1:N
%         n(i) = floor(sqrt(2)^(floor(2*log2(n(i))) + 1)) + 1;
%     end
% end
% end
%         if any(r(i)>1)
%             if r(i)<3
%                  r(i) = max(6,2*r(i));
%             end
%         end
%         r(i) = max(r(i),3);
%     end
%         break
%     elseif happy == 0
%         break
%     end
% end
end
