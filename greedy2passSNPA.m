function [J,H] = greedy2pass(M,r,normalize) 

% Successive Nonnegative Projection Algorithm (variant with f(.) = ||.||^2)
%
% *** Description ***
% At each step of the algorithm, the column of M maximizing ||.||_2 is 
% extracted, and M is updated with the residual of the projection of its 
% columns onto the convex hull of the columns extracted so far. 
% 
% See N. Gillis, Successive Nonnegative Projection Algorithm for Robust 
% Nonnegative Blind Source Separation, arXiv, 2013. 
%  
%
% [J,H] = SNPA(M,r,normalize) 
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W is full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
%
% ****** Output ******
% J        : index set of the extracted columns. 
% H        : optimal weights, that is, H argmin_{X >= 0} ||M-M(:,K)X||_F

[m,n] = size(M); 
maxitn = 100; 

if nargin <= 2, normalize = 0; end
if normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags(((sum(M)+1e-16).^(-1))', 0, n, n); M = M*D; 
end

% normM = sum(M.^2); 
% nM = max(normM); 

i = 1; 
% Perform r recursion steps (unless the relative approximation error is 
% smaller than 10^-9)
R = M;

% M
% J = [0];

M_normed = M*diag(1./sqrt(sum(M.^2)));
% Perform r recursion steps (unless the relative approximation error is 
% smaller than 10^-9)
while i <= r %&& max(normM)/nM > 1e-9

    % disp(['Greedy 1st Pass iter = ' int2str(i) '/' int2str(r) '...']);
    [a,b] = max(sum(max(R'*M_normed,0)));
    J(i) = b;
    
    % Update residual 
    if i == 1
        % Initialization using 10 iterations of coordinate descent
        H = nnlsHALSupdt(M,M(:,J),[],10); 
        % Fast gradient method for min_{y in Delta} ||M(:,i)-M(:,J)y||
        H = FGMfcnls(M,M(:,J),[],maxitn); 
    else
        H(:,J(i)) = 0; 
        h = zeros(1,n); h(J(i)) = 1; 
        H = [H; h]; 
        H = FGMfcnls(M,M(:,J),H,maxitn); 
    end
    R = M - M(:,J)*H; 
    
%     % Update norms
%     normM = sum(R.^2);
    
%     b
%     best
%     actual_dec = oldNorm - sum(normM)
%     best_dec
    
    i = i + 1; 
end

i = i-1; % which is = r if loop doesn't exit early
if i > 1
    for jj = 1:r

        % disp(['Greedy 2nd Pass iter = ' int2str(jj) '/' int2str(r) '...']);
        
        J = J(2:i);
        H(1:(end-1), :) = nnlsHALSupdt(M,M(:,J),H(2:end, :),maxitn);
        R = M - M(:,J)*H(1:(end-1), :);

        [a,b] = max(sum(max(R'*M_normed,0)));
        J(i) = b;


        H(:,J(i)) = 0;
        h = zeros(1,n); h(J(i)) = 1;

        H(end, :) = h;
        %H = FGMfcnls(M,M(:,J),H,maxitn);
    end
end
H = nnlsHALSupdt(M,M(:,J),H,maxitn);
end % of function FastSepNMF
