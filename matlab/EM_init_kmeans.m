function [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates)
%
% This function initializes the parameters of a Gaussian Mixture Model 
% (GMM) by using k-means clustering algorithm.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:     D x N array representing N datapoints of D dimensions.
%   o nbStates: Number K of GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       D x K array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the 
%               K GMM components.
% Comments ---------------------------------------------------------------
%   o This function uses the 'kmeans' function from the MATLAB Statistics 
%     toolbox. If you are using a version of the 'netlab' toolbox that also
%     uses a function named 'kmeans', please rename the netlab function to
%     'kmeans_netlab.m' to avoid conflicts. 

[nbVar, nbData] = size(Data);

%Use of the 'kmeans' function from the MATLAB Statistics toolbox
%Empty action changed to singleton - DHG
%[Data_id, Centers] = kmeans(Data', nbStates,'EmptyAction','singleton'); 
%For octave, need diff kmeans
[Centers, pic,Data_id]=myKmeans(Data',nbStates);
Mu = Centers';
for i=1:nbStates
  idtmp = find(Data_id==i);
  Priors(i) = length(idtmp);
  Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
  %Add a tiny variance to avoid numerical instability
  Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
end
Priors = Priors ./ sum(Priors);


function [centroid, pointsInCluster, assignment]= myKmeans(data, nbCluster)
% usage
% function[centroid, pointsInCluster, assignment]=
% myKmeans(data, nbCluster)
%
% Output:
% centroid: matrix in each row are the Coordinates of a centroid
% pointsInCluster: row vector with the nbDatapoints belonging to
% the centroid
% assignment: row Vector with clusterAssignment of the dataRows
%
% Input:
% data in rows
% nbCluster : nb of centroids to determine
%
% (c) by Christian Herta ( www.christianherta.de )
%
%Modified to run in matlab also by Dang
data_dim = length(data(1,:));
nbData   = length(data(:,1));


% init the centroids randomly
data_min = min(data);
data_max = max(data);
data_diff = data_max - data_min ;
% every row is a centroid
centroid = ones(nbCluster, data_dim) .* rand(nbCluster, data_dim);
for i=1 : 1 : length(centroid(:,1))
  centroid( i , : ) =   centroid( i , : )  .* data_diff;
  centroid( i , : ) =   centroid( i , : )  + data_min;
end
% end init centroids



% no stopping at start
pos_diff = 1.;

% main loop until
while pos_diff > 0.0

  % E-Step
  assignment = [];
  % assign each datapoint to the closest centroid
  
  nbdata = length(data(:,1));
  
  min_diff = data - repmat(centroid(1,:),nbdata,1); 
  min_diff = sum(min_diff.^2,2);
  curAssignment = ones(nbdata,1);
  for c=2:nbCluster
    diff2c = data - repmat(centroid(c,:),nbdata,1);
    diff2c = sum(diff2c.^2,2);
    idx = find(min_diff >= diff2c);
    curAssignment(idx) = c;
    min_diff(idx) = diff2c(idx);
  end
  
  assignment = curAssignment;

  % for the stoppingCriterion
  oldPositions = centroid;

  % M-Step
  % recalculate the positions of the centroids
  centroid = zeros(nbCluster, data_dim);
  pointsInCluster = zeros(nbCluster, 1);

  for c = 1:nbCluster
    idx = find(assignment == c);
    centroid(c,:) = mean(data(idx,:));
    if length(idx) == 0
      centroid( c , : ) = (rand( 1, data_dim) .* data_diff) + ...
	  data_min;
    end
  end

  %stoppingCriterion
  pos_diff = sum (sum( (centroid - oldPositions).^2 ) );

end
