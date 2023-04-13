function [cen, s] = SolveCentroidLinear(pt,C)
% Solve the centroid of the unknown sphere given diagonal matrix C and 
% input data pt according to Lemma 3, 8 and Theorem 4, 6, 7 of [2].  
% Step 3 in Alg. 1 of [1].
%
% Input:
%           pt:             n*m, n data described in m variables
%           C:              n*n, a diagonal scaling matrix that scales the 
%                           EDM A of pt to a doubly stochastic matrix B, 
%                           B=CAC
% Output:
%           cen:            centroid of the sphere of [1]
%           s:              E_p=sA, where E_p is the EDM of the inverse 
%                           stereographic projections on the sphere
%
% [1] Li He and Hong Zhang, Doubly Stochastic Distance Clustering, to
% appear on IEEE Transactions on Circuits and Systems for Video Technology
% (TCSVT).
%
% [2] C. R. Johnson, R. D. Masson, and M. W. Trosset. On the diagonal 
% scaling of euclidean distance matrices to doubly stochastic matrices. 
% Linear Algebra and its Applications, vol. 397, pp. 253¨C264, 2005.
%
% hel@sustech.edu.cn
% 2023.4.13

[n,~] = size(pt);

r_B = sqrt(.5/n);

P = [pt r_B./diag(C)];
H = [2*P -ones(n,1)];

b = sum(P.^2,2)-(r_B./diag(C)).^2;
if issparse(H)
    Y = pinv_sparse(H'*H)*H'*b;
else
    Y = pinv(H'*H)*H'*b;
end


% v = Y(end);
Z = Y(1:end-2);

% 2r_Bs^1/2
r_Bs = Y(end-1);

% center z
cen = Z';

% s
sqrts = r_Bs/2/r_B;
s = sqrts^2;

% radius r=r_B*s^1/2, and add radius as the (d+1) dimension of the center
cen = [cen r_B*sqrts];