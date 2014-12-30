% Manjari Akella
% CSE5524 - HW5
% 09/30/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
% Load the data
clc;
clear; 
close all;
load 'given_data/eigdata.txt';
X = eigdata;
figure('Name','Q1: Data','NumberTitle','off'),subplot(2,1,1);
plot(X(:,1),X(:,2),'b.');
axis('equal');
title('Raw Data');
% mean-subtract data
m = mean(X);
Y = X - ones(size(X,1),1)*m;
subplot(2,1,2);
plot(Y(:,1),Y(:,2),'r.');
axis('equal');
title('Mean Subtracted Data');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2
CY = cov(Y);
c=9;
[U,V] = eig(CY);
disp('Eigen Values');
disp(V);
disp('Eigen Vectors');
disp(U);
% Minor Axis 
minaxis = sqrt(c*V(1,1)).*U(:,1);
% Major Axis 
majaxis = sqrt(c*V(2,2)).*U(:,2);
% Plot data
figure('Name','Q2: Data with 3*sigma length of axis','NumberTitle','off'),plot(Y(:,1),Y(:,2),'r.');
axis('equal');
hold on;
line([0,minaxis(1,1)],[0,minaxis(2,1)]);
line([0,majaxis(1,1)],[0,majaxis(2,1)]);
ellipse(sqrt(c*V(1,1)),sqrt(c*V(2,2)),-atan((majaxis(2,1)/majaxis(1,1))),0,0,'g');

% 1*sigma to 5*sigma plots
figure('Name','Q2: Data with n*sigma length of axis','NumberTitle','off'),plot(Y(:,1),Y(:,2),'r.');
axis('equal');
hold on;
for i=1:5
    c=i^2;
    % Minor Axis
    minaxis = sqrt(c*V(1,1)).*U(:,1);
    % Major Axis
    majaxis = sqrt(c*V(2,2)).*U(:,2);
    ellipse(sqrt(c*V(1,1)),sqrt(c*V(2,2)),-atan((majaxis(2,1)/majaxis(1,1))),0,0,'b');
end
hold off;
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3
c=9;
for i=1:size(Y,1)
    new_Y(i,:) = (U')*(Y(i,:)');
end
% Plot data
figure('Name','Q3: Rotated Data','NumberTitle','off'),plot(new_Y(:,1),new_Y(:,2),'r.');
axis('equal');
hold on;
line([0,-sqrt(c*V(1,1))],[0,0]);
line([0,0],[0,sqrt(c*V(2,2))]);
ellipse(sqrt(c*V(1,1)),sqrt(c*V(2,2)),0,0,0,'g');
hold off;
figure('Name','Q3: Histogram of uncorrelated data','NumberTitle','off'),hist(new_Y(:,2),20);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4
clc;
clear all;
close all;
load 'given_data/xtraj.txt';
load 'given_data/ytraj.txt';
figure('Name','Q4: X,Y Trajectories','NumberTitle','off'),plot(xtraj, ytraj, 'r.-');
axis('equal');
xlabel('X'); ylabel('Y');
figure('Name','Q4: X,Y Trajectories vs T','NumberTitle','off'),plot(xtraj, 'r.-');
hold on;
plot(ytraj, 'b.-');
hold off;
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 5
Ga = [1,4,6,4,1];
G = (Ga./16)';
xtraj_s = imfilter(xtraj,G,'replicate');
ytraj_s = imfilter(ytraj,G,'replicate');
figure('Name','Q5: Smoothed X,Y Trajectories','NumberTitle','off'),plot(xtraj_s, ytraj_s, 'r.-');
axis('equal');
figure('Name','Q5: Smoothed X,Y Trajectories vs T','NumberTitle','off'),plot(xtraj_s, 'r.-');
hold on;
plot(ytraj_s, 'b.-');
hold off;
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 6
[k]=spatioTemporalC(xtraj,ytraj);
figure('Name','Q6: Curvature(Unsmoothed Trajectories)','NumberTitle','off'), plot(k,'b');
[k_s]=spatioTemporalC(xtraj_s,ytraj_s);
figure('Name','Q6: Curvature(Smoothed Trajectories)','NumberTitle','off'), plot(k_s,'g');
hold on;
[f,loc] = findpeaks(k_s);
peaks(:,1) = f';
peaks(:,2) = loc';
m = max(peaks);
thresh = m(1,1);
for i=1:size(peaks,1)
    if(peaks(i,1)>=(0.9*thresh))
       plot(peaks(i,2),peaks(i,1),'r.');
    end
end
hold off;

%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
