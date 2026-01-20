%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using the spectral clustering algorithm of 
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points 
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------
function groups = SpectralClustering(S,n)
%S为相似度矩阵，n为类的个数
warning off;
N = size(S,1);
% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
DN = diag( 1./sqrt(sum(S)+eps) );
LapN = speye(N) - DN * S * DN;
[uN,sN,vN] = svd(LapN);
P = vN(:,N-n+1:N);
for i = 1:N
   P(i,:) = P(i,:) ./ norm(P(i,:)+eps);
end
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans
groups = kmeans(P,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
end

