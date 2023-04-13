function demo
% Demostration of Doubly Stochastic Distance Clustering in [1].
%
% [1] Li He and Hong Zhang, Doubly Stochastic Distance Clustering, to
% appear on IEEE Transactions on Circuits and Systems for Video Technology
% (TCSVT).
%
% hel@sustech.edu.cn
% 2023.4.13
%% 1. Initialization
clc
clear
close all

load Iris.mat; % data, labels
numClasses = size(unique(labels),1);

% Euclidean Distance Matrix
A = pdist2(data,data).^2;

%% 2. Doubly Stochastic Distance Clustering
% Sinkhorn scaling, step 2 of Alg. 1 in [1]
[~, C, ~] = ScaleToDoublyStochasticMatrix(A);

% Solve the center of the sphere, step 2 of Alg. 1 [1]
[cen, ~] = SolveCentroidLinear(data,C);

% Inverse Stereographic Projections P, step 4 of Alg. 1 in [1]
P = inverseStereographicProjection(data, cen(end), cen(1:end-1));

%% 3. Final clustering
lab = kmeans(P,numClasses);
