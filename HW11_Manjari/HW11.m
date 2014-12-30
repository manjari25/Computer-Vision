% Manjari Akella
% CSE5524 - HW11
% 11/04/2013
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
clear all;
close all;
clc;
% Define k
k=[1,5,11,15];
% Load train and test files
tr = load('given_data/train.txt');
te = load('given_data/test.txt');
% Define number of training and testing instances
Ntr = size(tr,1);
Nte = size(te,1);
fprintf('KNN\n');
% Loop for each k
for itr=1:size(k,2)
    % For classes, to count no. of misclassifications, distance
    count = 0;
    dist = [];
    class = [];
    for i=1:Nte
        for j=1:Ntr
            % Euclidian Distance   
            dist(1,j) = sqrt(sum(((te(i,1:2)-tr(j,1:2)).^2),2));
        end
        % Sort in ascending order
        [sdist,idx] = sort(dist,'ascend');
        % Nearest k neighbours indices
        idx = idx(1:k(itr));
        % Get classification label of nearest neighbours             
        nearest = tr(idx,3);
        % Assign mode of labels of nearest neighbours to the test point
        class(i,1) = mode(nearest);
    end
    figure('Name',strcat('Q1: k= ',num2str(k(itr))),'NumberTitle','off');
    hold on;
    for i=1:Nte
            % Plot the points according to class label
            if(class(i,1)==1)
                plot(te(i,1),te(i,2),'b.');
            else
                plot(te(i,1),te(i,2),'r.');
            end
            
            % (Re)-plot misclassified points
            if(class(i,1)~=te(i,3))
                plot(te(i,1),te(i,2),'ko');
                count = count+1;
            end
    end
    hold off;
    accuracy = (Nte-count)/Nte;
    fprintf('k=%d,accuracy=%f\n',k(itr),accuracy);
end
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2
clear all;
close all;
clc;
% Load train and test files and define test data size
tr = load('given_data/train.txt');
te = load('given_data/test.txt');
Nte = size(te,1);
% Level of pruning
L = [0,2,9,10,11];
fprintf('D-TREE\n');
% Loop for each L
for itr=1:size(L,2)
    tree = ClassificationTree.fit(tr(:,1:2),tr(:,3)); 
    tprune = prune(tree,'level',L(itr)); 
    class = predict(tprune,te(:,1:2));
    view(tprune,'mode','graph');
    figure('Name',strcat('Q2: L= ',num2str(L(itr))),'NumberTitle','off');
    hold on;
    % To count number of misclassified instances
    count=0;
    for i=1:Nte
            % Plot the points according to class label
            if(class(i,1)==1)
                plot(te(i,1),te(i,2),'b.');
            else
                plot(te(i,1),te(i,2),'r.');
            end
            
            % (Re)-plot misclassified points
            if(class(i,1)~=te(i,3))
                plot(te(i,1),te(i,2),'go');
                count = count+1;
            end
    end
    hold off;
    accuracy = (Nte-count)/Nte;
    fprintf('L=%d,accuracy=%f\n',L(itr),accuracy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%