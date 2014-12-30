% Manjari Akella
% CSE5524 - HW3
% 09/16/2013

mkdir('Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1

refresh(); 
R = double(imread('given_pics/bg000.bmp'));
I = double(imread('given_pics/walk.bmp'));
% Computes a threshold for the image
normlevel = graythresh(I);
% Convert normalized value from [0-1] to [0-255] range
level = (normlevel*(max(max(I))-min(min(I))))-min(min(I));
% Check the image for threshold values for -50 to +20 of the returned threshold
% value
figure('Name','Q1: Background Subtraction 1 Thresholds','NumberTitle','off');
for i=-50:10:20
    fprintf('Threshold=%d\n',level+i);
    bs1Im = (abs((R-I)))>(level+i);
    imagesc(bs1Im);
    colormap('gray');
    pause(1);
end
% Pick Threshold = 76.0039
bs1Im = (abs((R-I)))>(76.0039);
imwrite(bs1Im,'Output/Q1_B_Sub1.bmp');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2

refresh();
I = double(imread('given_pics/walk.bmp'));
% Read all background images into a 3D matrix (Each slice is an
% image(240x320)
N=30;
for  i=1:N
filename = sprintf('given_pics/bg%03d.bmp', i-1);
Im(:,:,i) = double(imread(filename));
end

% Model mean and standard deviation statistics for all background images
meanIm = mean(Im,3);
sdIm = std(Im,0,3);

% Threshold the image to subtract the background
figure('Name','Q2: Background Subtraction 2 Thresholds','NumberTitle','off');
for i=3:20
    fprintf('Threshold=(%d*sigma)^2\n',i);
    bs2Im =(((I-meanIm).^2)./(sdIm.^2))>((i.*sdIm).^2);
    imagesc(bs2Im);
    colormap('gray');
    axis('image');
    pause(1);
end
% Pick Threshold = (20*sigma)^2
bs2Im =(((I-meanIm).^2)./(sdIm.^2))>((20.*sdIm).^2);
imwrite(bs2Im,'Output/Q2_B_Sub2.bmp');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3

d_bs2Im = bwmorph(bs2Im, 'dilate');
figure('Name','Q3: Dialated Image','NumberTitle','off'),imagesc(d_bs2Im);
colormap('gray');
axis('image');
imwrite(d_bs2Im,'Output/Q3_Dilated_B_Sub2.bmp');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4

% Label the regions
[L, num] = bwlabel(d_bs2Im, 8);
% Find out more about the connected components
CC = bwconncomp(L);
% Find number of pixels in each label
numPixels = cellfun('size',CC.PixelIdxList,1);
% Find label with maximum number of pixels
[biggest,idx] = max(numPixels);
% Label those indices 1, rest all 0
for i=1:size(L,1)
    for j=1:size(L,2)
        if (L(i,j)==idx)
            L(i,j) = 1;
        else
            L(i,j) = 0;
        end
    end
end
figure('Name','Q4: Largest Component','NumberTitle','off'),imagesc(L);
colormap('gray');
axis('image');
imwrite(L,'Output/Q4_Largest_Connected_Component.bmp');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 5

refresh();
Im = im2bw(double(imread('given_pics/boxIm1.bmp')));
% 'Basic' returns area,centroid and bounding box
s = regionprops(Im,'basic');
c = s.Centroid;
a = s.Area;
bb = s.BoundingBox;
% Coordinates displayed are w.r.t image axis
fprintf('Area = %d (pixels)\n',a);
fprintf('Centroid Coordinates (column average,row average)  = (%d,%d)\n',c(1,1),c(1,2));
fprintf('Bounding Box (x,y,width,height) = (%d,%d,%d,%d)\n',bb(1,1),bb(1,2),bb(1,3),bb(1,4));
figure('Name','Q5: Centroid & Bounding Box','NumberTitle','off'),imagesc(Im);
colormap('gray');
axis('image');
hold on
rectangle('Position',bb,'EdgeColor','g','LineWidth',1);
plot(s.Centroid(1,1), s.Centroid(1,2), 'r*');
pause;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 6

refresh();
% Loop for each box image
for i=1:4
    Im =  im2bw(double(imread(strcat('given_pics/boxIm',num2str(i),'.bmp'))));
    % Find out the seven moment descriptors
    n02 = similitudeMoments(Im,0,2);
    n03 = similitudeMoments(Im,0,3);
    n11 = similitudeMoments(Im,1,1);
    n12 = similitudeMoments(Im,1,2);
    n20 = similitudeMoments(Im,2,0);
    n21 = similitudeMoments(Im,2,1);
    n30 = similitudeMoments(Im,3,0);
    fprintf('boxIm%d\n',i);
    N = [n02, n03, n11, n12, n20, n21, n30];
    disp(N);
end
%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
