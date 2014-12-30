function [ q ] = computeHist( info,center)
% Takes the info matrix and center coordinates of the patch
% and returns a histogram distribution cube (16x16x16) of the 
% patch under consideration
q_u = zeros(16,16,16);
countR = zeros(1,16);
countG = zeros(1,16);
countB = zeros(1,16);
% Model q_u for n channel
for i=1:size(info,2)
    j1 = ceil(info(3,i)/16);
    j2 = ceil(info(4,i)/16);
    j3 = ceil(info(5,i)/16);
    countR(1,j1) = countR(1,j1)+1;
    countG(1,j2) = countG(1,j2)+1;
    countB(1,j3) = countB(1,j3)+1;
    x = [info(1,i),info(2,i)];
    % Unrounded center
    t = ((center-x)./25);
    t1 = (sqrt(t(1,1)^2+t(1,2)^2))^2;
    q_u(j1,j2,j3) = q_u(j1,j2,j3)+functionK(t1);
end
% Normalize
q = q_u./(sum(sum(sum(q_u))));
