function spt = inverseStereographicProjection(pt, radius, center)
% Inverse stereographic projection used in [1]. Project points on the d-dim 
% plane to a (d+1)-dim sphere by Lemma 5 of [2].
%
% Taking d=2 for example. Given a sphere S with center (x,y,r) and the 
% radius r, and 2D data pt on the x=0 plane. Let L be the line linking the 
% north pole of S and pt, and the intersecting points of L and S is spt.
%
% Input:
%       pt:         points on the d-dim plane
%       radius:     radius of sphere, shpere is S(center, radius)
%       center:     center of shpere, shpere is S(center, radius)
% Output:
%       spt:        points on (d+1)-dim sphere
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

[n,dim] = size(pt);
if nargin<2
    radius = (.5/n)^.5;
    center = zeros(1,dim);
elseif nargin<3
    center = zeros(1,dim);
end
% since Lemma 5 requires the sphere lies on [0,0,...,0,r], so move
% points and let the center be the original
cenOffset = -center;
cenOffset(dim+1) = 0;

% move data with offset, now, center of sphere is [0,...0,r]
A = pt+repmat(cenOffset(1:dim),n,1);

% Lemma 5
t  =1./( 4*radius^2+sum(A.^2,2) );  % denominator of Lemma 5

spt = 4*radius^2*A.*repmat(t,1,dim);    % dim: 1~d
spt(:,dim+1) = 2*radius*sum(A.^2,2).*t; % dim: d+1

spt = spt-repmat(cenOffset,n,1);