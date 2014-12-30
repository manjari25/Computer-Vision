function [ w ] = computeWeight( q_u,p_u,info )
% Takes 2 probability histograms and the info matrix as inputs
% Returns a weight value for each pixel in info
for i=1:size(info,2)
    j1 = ceil(info(3,i)/16);
    j2 = ceil(info(4,i)/16);
    j3 = ceil(info(5,i)/16);
    w(1,i) = sqrt(q_u(j1,j2,j3)/p_u(j1,j2,j3));
end

