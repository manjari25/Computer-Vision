% Manjari Akella
% CSE5524 - HW2
% 09/09/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1

refresh();          
sigma = input('Enter the value of sigma\n');
[Gx,Gy] = gaussDeriv2D(sigma);
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2

close all;
%  For checkerboard image and a group photograph
pics = {'my_pictures/1.bmp','my_pictures/2.jpg'};
for i=1:2
    close all;
%     Group photograph is not grayscale
    if (i==2)
        GIm = rgb2gray(imread(pics{i}));
        Im = double(GIm);
    else
        GIm = imread(pics{i});
        Im = double(GIm);
    end
    figure('Name','Q2:Original Image','NumberTitle','off'),imagesc(uint8(Im));
    axis('image');
    colormap('gray');
    gxIm = imfilter(Im,Gx,'replicate');
    figure('Name','Q2:x-gradient','NumberTitle','off'),imagesc(gxIm);
    axis('image');
    colormap('gray');
    gyIm = imfilter(Im,Gy,'replicate');
    figure('Name','Q2:y-gradient','NumberTitle','off'),imagesc(gyIm);
    axis('image');
    colormap('gray');
    magIm = sqrt(gxIm.^2+gyIm.^2);
    figure('Name','Q2:Gradient Magnitude','NumberTitle','off'),imagesc(magIm);
    axis('image');
    colormap('gray');
    pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3

close all;
max_val = max(max(magIm));
% Greater than percentage of max value
for i=0.1:0.2:0.9
    tIm = magIm>(i*max_val);
    figure('Name',strcat('Q3:Threshold ',num2str(i*100),' percent'),'NumberTitle','off'),imagesc(uint8(tIm));
    colormap('gray');
    axis('image');
end
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4

refresh();
%  For checkerboard image and a group photograph
pics = {'my_pictures/1.bmp','my_pictures/2.jpg'};
for i=1:2
    close all;
%     Group photograph is not grayscale
    if (i==2)
        GIm = rgb2gray(imread(pics{i}));
        Im = double(GIm);
    else
        GIm = imread(pics{i});
        Im = double(GIm);
    end
    fx = -fspecial('sobel')';
    fxIm = imfilter(Im,fx,'replicate');
    figure('Name','Q4:Sobel x-direction','NumberTitle','off'),imagesc(fxIm);
    colormap('gray');
    axis('image');
    fy = fx';
    fyIm = imfilter(Im,fy,'replicate');
    figure('Name','Q4:Sobel y-direction','NumberTitle','off'),imagesc(fyIm);
    colormap('gray');
    axis('image');
    magIm = sqrt(fxIm.^2+fyIm.^2);
    figure('Name','Q4:Sobel Gradient Magnitude','NumberTitle','off'),imagesc(magIm);
    colormap('gray');
    axis('image');
    pause;
    close all;
end
    max_val = max(max(magIm));
    % Greater than percentage of max value
    for i=0.1:0.2:0.9
        tIm = magIm>(i*max_val);
        figure('Name',strcat('Q4:Threshold ',num2str(i*100),' percent'),'NumberTitle','off'),imagesc(tIm);
        colormap('gray');
        axis('image');
    end
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 5

refresh();

% Group photo (more detail in image)
% Im = double(rgb2gray(imread('my_pictures/2.jpg')));

% Less detail in image
Im = double(rgb2gray(imread('my_pictures/3.jpg')));
EGIm = edge(Im,'canny');
figure('Name','Q5_Canny default result','NumberTitle','off'),imagesc(EGIm);
colormap('gray');
axis('image');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 6

refresh();
a=0.4;
w1=[0.25-0.5*a,0.25,a,0.25,0.25-0.5*a];
w2 = w1';

%Reduce operation
% Level 1
Im=double(rgb2gray(imread('my_pictures/4.jpg')));
figure('Name','Q6_Gaussian_Pyramid','NumberTitle','off'),subplot(2,2,1),imagesc(Im);
colormap('gray');
axis('image');

% Level 2
I2 = reduce_G(Im,w1,w2);
subplot(2,2,2),imagesc(I2);
colormap('gray');
axis('image');

% Level 3 
I3 = reduce_G(I2,w1,w2);
subplot(2,2,3),imagesc(I3);
colormap('gray');
axis('image');

% Level 4 
I4 = reduce_G(I3,w1,w2);
subplot(2,2,4),imagesc(I4);
colormap('gray');
axis('image');

% Expand operation
% Level 4
E4 = I4;
figure('Name','Q6_Laplacian_Pyramid','NumberTitle','off'),subplot(2,2,1),imagesc(E4);
colormap('gray');
axis('image');

% Level 3
I = expand(I4);
E3 = I3-I;
subplot(2,2,2),imagesc(E3);
colormap('gray');
axis('image');

% Level 2
I = expand(I3);
E2 = I2-I;
subplot(2,2,3),imagesc(E2);
colormap('gray');
axis('image');

% Level 1
I = expand(I2);
E1 = Im-I;
subplot(2,2,4),imagesc(E1);
colormap('gray');
axis('image');

% Reconstruction
L3 = expand(E4);
temp = L3+E3;
L2 = expand(temp);
temp = L2+E2;
L1 = expand(temp);
result = L1+E1;

figure('Name','Q6_Reconstructed','NumberTitle','off'),imagesc(result);
colormap('gray');
axis('image');
figure('Name','Q6_Original','NumberTitle','off'),imagesc(Im);
colormap('gray');
axis('image');

%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
