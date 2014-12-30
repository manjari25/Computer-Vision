% Manjari Akella
% CSE5524 - HW6
% 10/07/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
clear all;
close all;
clc;
% Covariance matrix of model
modelCovMatrix = [47.917 0 -146.636 -141.572 -123.269;
0 408.250 68.487 69.828 53.479;
-146.636 68.487 2654.285 2621.672 2440.381;
-141.572 69.828 2621.672 2597.818 2435.368;
-123.269 53.479 2440.381 2435.368 2404.923];
% Read image
I = double(imread('given_pics/target.jpg'));
r=1;
cfk=zeros(5,5);
while ((r+69)<=240)
    c=1;
    while ((c+23)<=320)
        % Window of 70x24        
        Im = I(r:r+69,c:c+23,1:3);
        n=1;
        % Compute feature vector for window
        for i=1:size(Im,1)
            for j=1:size(Im,2)
                fk(1,n) = j;
                fk(2,n) = i;
                fk(3,n) = Im(i,j,1);
                fk(4,n) = Im(i,j,2);
                fk(5,n) = Im(i,j,3);
                n=n+1;
            end
        end
        % Compute covariance          
        candCovMatrix  = cov(fk');
        % Compute score with model matrix
        score = compare(modelCovMatrix,candCovMatrix);
        % Store score value
        scoreMatrix(r,c) = score;
        c=c+1;
    end
    r=r+1;
end
bestMatch = min(min(scoreMatrix));
[m,ind] = min(scoreMatrix(:));
[x,y] = ind2sub(size(scoreMatrix),ind);
figure('Name','Q1: Score Matrix Surface Plot','NumberTitle','off'),surf(scoreMatrix);
figure('Name','Q1: Best Match Window','NumberTitle','off'),imagesc(uint8(I));
axis('image');
hold on;
rectangle('Position',[y,x,24,70],'LineWidth',2,'EdgeColor',[0.498039, 1, 0]);
hold off;
fprintf('Best score = %f',bestMatch);
fprintf('\n(r,c) of best match = (%d,%d)',x,y);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2
clear all;
close all;
clc;
% Read first frame image
I1 = double(imread('given_pics/img1.jpg'));
center = [150,175];
figure('Name','Q2: Frame 1, Target Model','NumberTitle','off'),imagesc(uint8(I1));
axis('image');
hold on;
ellipse(25,25,0,150,175,'b');
hold off;
[info] = circRegion25(I1,center);
% Compute histogram for R,G,B channels
[q_u] = computeHist(info,center); 
% New frame image
I2 = double(imread('given_pics/img2.jpg'));
center1  = center;
figure('Name','Q2: Frame 2, Target Candidate','NumberTitle','off');
imagesc(uint8(I2));
axis('image');
hold on;
ellipse(25,25,0,center(1,1),center(1,2),'b');
num=0;
di = zeros(25,1);
for i=1:25
    % Compute weights
    [info] = circRegion25(I2,round(center1));
    [p_y] = computeHist(info,center1); 
    [w] = computeWeight(q_u,p_y,info);
    
    % Find next best location
    num(1,1) = sum(info(1,:).*w);
    num(1,2) = sum(info(2,:).*w);
    den = sum(w);
    new_center = num./den;
    % Distance between consecutive iterations
    di(i,1) = sqrt((center1(1,1)-new_center(1,1))^2 + (center1(1,2)-new_center(1,2))^2);
    % Set center1 = new_center
    center1 = new_center;
    % Plot the new patch    
    h = ellipse(25,25,0,round(center1(1,1)),round(center1(1,2)),'r');
    pause(0.5);
    % If not the last iteration, delete circle plot
    if(i~=25)
        delete(h);
    end
end
hold off;
% Final position
fprintf('Final Position(x,y)=(%.10f,%.10f)',center1(1,1),center1(1,2));
% fprintf('Final Position(x,y)=(%f,%f)',center1(1,1),center1(1,2));

% Distance between last 2 iterations
fprintf('\nDistance between last two iterations(x,y)=%.10f',di(size(di,1),1));
% fprintf('\nDistance between last two iterations(x,y)=%f',di(size(di,1),1));

%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
