% Manjari Akella
% CSE5524 - HW4
% 09/23/2013

mkdir('Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
clear all;
close all;
clc;
N=22;
for  i=1:N
filename = sprintf('given_pics/aerobic2.%03d.bmp', i);
Im(:,:,i) = double(imread(filename));
end
figure('Name','Q1: Image Differencing','NumberTitle','off');
for i=2:N
    I(:,:,(i-1)) = abs(Im(:,:,i)-Im(:,:,i-1))>15;
    diff_I(:,:,(i-1)) = bwareaopen(I(:,:,(i-1)),140);
    imwrite(diff_I(:,:,(i-1)),strcat('Output/diff_I',num2str(i-1),'.bmp'));
    subplot(3,7,i-1),imagesc(diff_I(:,:,(i-1)));
    colormap('gray');
    axis('image');
end
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2
% Compute Motion History Image
MHI = zeros(320,240);
for i=2:N
    MHI(diff_I(:,:,i-1)~=0)=(i-1);
    MHI(MHI<((i-1)-21))=0;
end
figure('Name','Q2: Motion History Image','NumberTitle','off'),imagesc(MHI);
colormap('gray');
axis('image');

% Compute Motion Energy Image
MEI = MHI>0;
figure('Name','Q2: Motion Energy Image','NumberTitle','off'),imagesc(MEI);
colormap('gray');
axis('image');

% Normalize Motion History Image
NMHI = normalize(MHI);

% Find out the seven moment descriptors for NMHI
n02 = similitudeMoments(NMHI,0,2);
n03 = similitudeMoments(NMHI,0,3);
n11 = similitudeMoments(NMHI,1,1);
n12 = similitudeMoments(NMHI,1,2);
n20 = similitudeMoments(NMHI,2,0);
n21 = similitudeMoments(NMHI,2,1);
n30 = similitudeMoments(NMHI,3,0);
fprintf('MHI\n');
N = [n02, n03, n11, n12, n20, n21, n30];
disp(N);
    
% Find out the seven moment descriptors for MEI
n02 = similitudeMoments(MEI,0,2);
n03 = similitudeMoments(MEI,0,3);
n11 = similitudeMoments(MEI,1,1);
n12 = similitudeMoments(MEI,1,2);
n20 = similitudeMoments(MEI,2,0);
n21 = similitudeMoments(MEI,2,1);
n30 = similitudeMoments(MEI,3,0);
fprintf('MEI\n');
N = [n02, n03, n11, n12, n20, n21, n30];
disp(N);

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3
clear all;
close all;
clc;
Im1 = double(imread('my_pictures/box.bmp')); 
Im2 = double(imread('my_pictures/box_shifted.bmp'));
%Set up filters
fx_s = (-fspecial('sobel')')/8;
fy_s = (fx_s)';
%Run filters on Image at time t
fx = imfilter(Im2,fx_s,'replicate');
fy = imfilter(Im2,fy_s,'replicate');
figure('Name','Q3: Filtered Image x-gradient','NumberTitle','off'),imagesc(fx);
colormap('gray');
axis('image');
figure('Name','Q3: Filtered Image y-gradient ','NumberTitle','off'),imagesc(fy);
colormap('gray');
axis('image');
% Compute ft
ft = Im2-Im1;
t = fx.^2+fy.^2;
d = sqrt(t);
% Compute vectors
for i=1:100
    for j=1:100
        if(d(i,j)~=0)
            u(i,j) = (-ft(i,j)/d(i,j))*(fx(i,j)/d(i,j));
            v(i,j) = (-ft(i,j)/d(i,j))*(fy(i,j)/d(i,j));
        else
            u(i,j) = 0;
            v(i,j) = 0;
        end
    end
end
figure('Name','Q3: Motion vector plot','NumberTitle','off'),quiver(u,v);
axis('ij');

%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
