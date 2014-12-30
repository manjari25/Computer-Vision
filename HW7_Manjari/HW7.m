% Manjari Akella
% CSE5524 - HW7
% 10/14/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
clear all;
close all;
clc;
load 'kfdata.txt'
%% state equation 
% state update 
G = [1]; 
% Process noise variance 
Q = .001; 
% Q = .001*100;
% Q = .001/100;

%% observation equation 
% transformation 
H = [1]; 
% observation noise variance 
R = .05; 
% R = .05*100;
% R = .05/100;

%% initial guesses 
x = 0; % initial state guess 
P = [1]; % initial state error variance guess 

% Estimate at t=1 without seeing observed value
x_(1,1) = G*x;
P_(1,1) = G*P*G' + Q;

for i=1:size(kfdata,1)
    K(i,1) = (P_(i,1)*H')*((H*P_(i,1)*H'+R)^-1);
    x(i,1) = x_(i,1) + (K(i,1)*(kfdata(i,1)-H*x_(i,1)));
    P(i,1) = (1-K(i,1))*P_(i,1);
    z_filtered(i,1) = H*x(i,1);
    % If i=100, no need to estimate next observation
    if(i~=size(kfdata,1))
        x_(i+1,1) = G*x(i,1);
        P_(i+1,1) = G*P(i,1)*G'+Q;
    end
end
figure;
subplot(311),plot(kfdata,'r');
title('Given data');
xlabel('t');
ylabel('kfdata');
hold on;
plot(1:size(kfdata,1),kfdata,'r*');
hold off;
subplot(312),plot(kfdata,'r');
title('Estimated mean and Error bar');
xlabel('t');
ylabel('kfdata/z_f');
hold on;
plot(1:size(kfdata,1),kfdata,'r*');
plot(1:size(z_filtered,1),z_filtered,'bo');
errorbar(z_filtered,sqrt(P),'b');
hold off;
subplot(313),plot(P,'g');
title('Error variance');
xlabel('t');
ylabel('P');

fprintf('Final (x,P) = %f,%f',z_filtered(100,1),P(100,1));
