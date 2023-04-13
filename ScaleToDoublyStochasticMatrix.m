function [B, C, err]= ScaleToDoublyStochasticMatrix(A, tol)
% Scale to doubly stochastic matrix with the alternative normalizing
% method. 
% Step 2 in Alg. 1 of [1].
%
% Input:
%           A:          non-negtive matrix, symmetric and zero diagonal
% Output:
%           B:          doubly stochasitc matrix A scaled to, B = CAC
%           C:          diagnal matrix, B = CAC
%           err:        sum(B-1).^2 + sum(B'-1).^2
%
% [1] Li He and Hong Zhang, Doubly Stochastic Distance Clustering, to
% appear on IEEE Transactions on Circuits and Systems for Video Technology
% (TCSVT).
%
% hel@sustech.edu.cn
% 2023.4.13

B = A;
if nargin==1
    tol = 1e-4;
end
maxiter = 1000; % maximum number of iterations
counter = 0;
err = tol*100;

a = ones(size(B,1),1);

% main loop
while err>tol && counter<maxiter
    % normalize by row-sum
    for i=1:size(B,1)
        k = sum(B(i,:));
        a(i) = a(i)/k;
        B(i,:) = B(i,:)/k;
        B(:,i) = B(:,i)/k;
    end
    % normalize by column-sum
    for j=1:size(B,2)
        k = sum(B(:,j));
        a(j) = a(j)/k;
        B(:,j) = B(:,j)/k;
        B(j,:) = B(j,:)/k;
    end
    
    % error: How close B is to a doubly stochastic matrix
    err = sum( (sum(B)-1).^2 + (sum(B')-1).^2 )^.5;

    counter = counter+1;
end

C = diag(a);