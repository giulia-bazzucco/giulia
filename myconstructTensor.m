function T = myconstructTensor (f,C)
%scrivi cosa dovrebbe essere f, C e output
% if vectorize == 1 
 k = numel(C); %number of elements in an array or subscripted array expression
    vec = cell(k,1); %create cell array
    n = cellfun(@numel, C); %apply a function to each cell of a cell array.
    for j = 1:k
        linind = 1:n(j);
        vec{j} = C{j}(linind);
    end
    [vec{:}] = ndgrid(vec{:});
    T = f(vec{:});
% else 
% 
% % Genera le griglie per le variabili
% grids = cell(1,N);
% for v = 1:N
%     % Definisci l'intervallo di campionamento per la variabile v-esima
%     % Esempio: intervallo da 0 a 1 con passo di 0.1
%     start = 0;
%     stop = 1;
%     step = 0.1;
%     grids{v} = start:step:stop;
% end
% 
% % Inizializza l'array T di output
% sizes = cellfun(@length, grids);
% T = zeros(sizes);
% 
% % Calcola i valori della funzione ff su tutti i punti della griglia
% input = cell(1,N);
% for i = 1:sizes(1)
%     input{1} = grids{1}(i);
%     for j = 1:sizes(2)
%         input{2} = grids{2}(j);
%         for k = 1:sizes(3)
%             input{3} = grids{3}(k);
%             % Esegui la chiamata alla funzione ff sui valori della griglia
%             T(i,j,k) = ff(input{:});
%         end
%     end
% end

end



%ALTRO MODO
%N = length(C);
% N=length(C);
% n1 = cellfun(@length,C);
% n = transpose(n1);
% disp(n);
% T = zeros(n);
%T = arrayfun(f,C);
%T = reshape(cell2mat(C),n);

% for i = 1 : prod(n)
%     linind1 = ind2sub(n, i);
%     linind = transpose(linind1);
% %     vec = cellfun(@(x) x(linind(:)),C);
% %     T(linind{:}) = f(vec);
% %     vec = cellfun(@(x) x(linind(i)),C);
% %     T(linind) = f(vec);
% 
% end
% end




%     N = length(C);
%     n1 = cellfun(@length, C);
%     n = transpose(n1);
%     T = reshape(cell2mat(C), n);
%     T = arrayfun(f, T);
%     
%     for i = 1 : prod(n1)
%         linind = ind2sub(n1, i);
%         vec = cellfun(@(x) x(linind(i)),C);
%         T(linind{:}) = f(vec);
%     end
% end



%
% N = length(C);
% n1 = cellfun(@length,C);
% n = transpose(n1);
% disp(n);
% T = zeros(n);
% %T = arrayfun(f,C);
% %T = reshape(cell2mat(C),n);
% 
% for i = 1 : prod(n)
%     linind = ind2sub(n, i);
%     
% %     vec = cellfun(@(x) x(linind(i)),C);
% %     T(linind{:}) = f(vec(1),vec(2),vec(3));
% 
%     vec = cellfun(@(x) x(linind(i)),C);
%     T(linind) = f(vec);
% end
% end



% N = length(C);
% n1 = cellfun(@length,C);
% n = transpose(n1);
% disp(n);
% T = zeros(n);
% %T = arrayfun(f,C);
% %T = reshape(cell2mat(C),n);
% 
% for i = 1 : prod(n)
%     linind = ind2sub(n, i);
%     
% %     vec = cellfun(@(x) x(linind(i)),C);
% %     T(linind{:}) = f(vec(1),vec(2),vec(3));
% 
%     vec = cellfun(@(x) x(linind(i)),C);
%     T(linind) = f(vec);
% end

%PRIMO MODO
% N = length(C);
% n = zeros(N,1);
% for i=1:N
%     n(i)=length(C{i});
% end
% 
% T = zeros(n); %T = arrayfun(f,C);
% 
% 
% for i = 1 : prod(n)
%     linind = ind2sub(n, i);
%     vec = zeros(N,1);
%     for j = 1 : N
%         vec(j) = C{j}(linind(j));
%     end
%     T(linind) = f(vec);
% end

%MIO MODO
% N = length(C);
% n = zeros(N,1);
% for i=1:N
%     n(i)=length(C{i});
% end
% %cellfun(length,J)
% %cheb = @(i,n) -cos((i-1).*pi/(n-1));
% % for i=1:N
% %      n(i)=length(C{i});
% 
% for i = 1:N
% T = zeros(n(i));
% end
% %T = zeros(n(i)); T = zeros(n); T = zeros(N)
%     for l = 1 : prod(n)
%       linind = ind2sub(n, l);
%       vec = zeros(N,1);
%      for j = 1 : N
%         vec(j) = C{j}(linind(j)); %cellfun(length,J)
%      end
%         T(linind) = f(vec);
%      end
% end
% 
