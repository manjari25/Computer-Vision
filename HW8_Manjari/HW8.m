% Manjari Akella
% CSE5524 - HW8
% 10/21/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1
clear all;
close all;
clc; 
%Define standard deviation values and alpha
sigmaD = 0.7;
sigma = 1;
a=0.05;
% Read image
I = double(imread('given_pics/checker.jpg'));
% Calculate gaussian derivatives
[Gx,Gy] = gaussDeriv2D(sigmaD);
% Run filter on image
Ix = imfilter(I,Gx,'replicate');
Iy = imfilter(I,Gy,'replicate');
% Compute Ix^2,Iy^2,Ix*Iy
Ix2 = Ix.^2;
Iy2 = Iy.^2; 
Ixy = Ix.*Iy;
% Weighting gaussian
WG = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
GIx2 = imfilter(Ix2,WG,'replicate');
GIy2 = imfilter(Iy2,WG,'replicate');
GIxy = imfilter(Ixy,WG,'replicate');
% Compute R
R = (GIx2.*GIy2)-((GIxy).^2)-(a.*((GIx2+GIy2).^2));
% Display R
Rimg = R - min(0, min(R(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q1: R','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
% Threshold R
R(R<1e6) = 0;
% Display R
Rimg = R - min(0, min(R(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q1: R(threshold=1e6)','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
% Non-maximal Suppression 
R_new = nlfilter(R, [3 3], @(x) all(x(5) >= x([1:4 6:9])));
RR = R_new.*R;
% Display RR
Rimg = RR - min(0, min(RR(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q1: R(threshold=1e6 + Non-maximal suppression)','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
% Show on original image
[r,c] = find(RR~=0);
figure('Name','Q1: Corner Points','NumberTitle','off');  
imagesc(I);
colormap gray;
axis('image');
hold on;
for i=1:size(r,1)
    plot(c(i),r(i),'r.');
end
hold off;
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2
clear all;
close all;
clc; 
%Define standard deviation values
sigmaD = 0.7;
sigma = 1;
% Read image and pre-allocate R
I = double(imread('given_pics/tower.jpg'));
R = zeros(size(I,1),size(I,2));
% Calculate gaussian derivatives
[Gx,Gy] = gaussDeriv2D(sigmaD);
% Run filter on image
Ix = imfilter(I,Gx,'replicate');
Iy = imfilter(I,Gy,'replicate');
% Compute Ix^2,Iy^2,Ix*Iy
Ix2 = Ix.^2;
Iy2 = Iy.^2; 
Ixy = Ix.*Iy;
% Weighting gaussian
WG = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
GIx2 = imfilter(Ix2,WG,'replicate');
GIy2 = imfilter(Iy2,WG,'replicate');
GIxy = imfilter(Ixy,WG,'replicate');
% Compute R (SHI-TOMASI)
for i=1:size(I,1)
    for j=1:size(I,2)
        M = [GIx2(i,j),GIxy(i,j);GIxy(i,j),GIy2(i,j)];
        d = eig(M);
        R(i,j) = min(d);
    end
end
% Display R
Rimg = R - min(0, min(R(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q2: R','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
% Non-maximal Supression
mask = nlfilter(R, [3 3], @(x) all(x(5) >= x([1:4 6:9])));
R_supp = R.*mask;
% Pick top 150 points
temp=sort(R_supp(:),'descend');
R_new=(R_supp>=temp(150)).*R_supp;
% Display
Rimg = R_new - min(0, min(R_new(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q2: R(Non-maximal suppression+top 150)','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
% Display
[r,c] = find(R_new~=0);
figure('Name','Q2/Q3: Corner Points(Shi-Tomasi+Harris)','NumberTitle','off'); 
subplot(121),imagesc(I);
colormap gray;
axis('image');
title('Shi-Tomasi Function');
hold on;
for i=1:size(r,1)
    plot(c(i),r(i),'b.');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3
% Compute R (HARRIS FUNCTION)
a=0.05;
RH = (GIx2.*GIy2)-((GIxy).^2)-(a.*((GIx2+GIy2).^2));
% Display
Rimg = RH - min(0, min(RH(:))); 
Rimg = Rimg/max(Rimg(:)); 
% Non-maximal Supression
maskH = nlfilter(RH, [3 3], @(x) all(x(5) >= x([1:4 6:9])));
R_suppH = RH.*maskH;
% Pick top 150 points
temp=sort(R_suppH(:),'descend');
R_newH=(R_suppH>=temp(150)).*R_suppH;
% Display
[rH,cH] = find(R_newH~=0);
subplot(122),imagesc(I);
colormap gray;
axis('image');
title('Harris Function');
hold on;
for i=1:size(rH,1)
    plot(cH(i),rH(i),'r.');
end
hold off;
figure('Name','Q3: R','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
Rimg = R_newH - min(0, min(R_newH(:))); 
Rimg = Rimg/max(Rimg(:)); 
figure('Name','Q3: R(Non-maximal Suppression+top 150)','NumberTitle','off'); 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
imshow(Rimg, 'Border', 'Tight'); 
colormap gray; 
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4
clc;
clear all;
close all;
% Load image 
Im = double(imread('given_pics/tower.jpg'));
% Define radius, no. of consecutive points, intensity thresholds
N = 9;
T = 10;
fastX=0;
fastY=0;
count = 1;
nr = size(Im,1)-4;
nc = size(Im,2)-4;
% For each pixel, check 16 pixel neighbourhood(1=above,-1=below,0=pixel
% doesn't exist or neither above nor below
for i=4:nr
        for j=4:nc
% 1 i-3,j
            if(Im(i-3,j)>(Im(i,j)+T))
                A(1,1)=1;
            elseif(Im(i-3,j)<(Im(i,j)-T))
                A(1,1)=-1;
            else
                A(1,1)=0;
            end
        
% 2 i-3,j+1
            if(Im(i-3,j+1)>(Im(i,j)+T))
                A(2,1)=1;
            elseif(Im(i-3,j+1)<(Im(i,j)-T))
                A(2,1)=-1;
            else
                A(2,1)=0;
            end

% 3 i-2,j+2
            if(Im(i-2,j+2)>(Im(i,j)+T))
                A(3,1)=1;
            elseif(Im(i-2,j+2)<(Im(i,j)-T))
                A(3,1)=-1;
            else
                A(3,1)=0;
            end
% 4 i-1,j+3
            if(Im(i-1,j+3)>(Im(i,j)+T))
                A(4,1)=1;
            elseif(Im(i-1,j+3)<(Im(i,j)-T))
                A(4,1)=-1;
            else
                A(4,1)=0;
            end
% 5 i,j+3
            if(Im(i,j+3)>(Im(i,j)+T))
                A(5,1)=1;
            elseif(Im(i,j+3)<(Im(i,j)-T))
                A(5,1)=-1;
            else
                A(5,1)=0;
            end
% 6 i+1,j+3
            if(Im(i+1,j+3)>(Im(i,j)+T))
                A(6,1)=1;
            elseif(Im(i+1,j+3)<(Im(i,j)-T))
                A(6,1)=-1;
            else
                A(6,1)=0;
            end
%  7 i+2,j+2
            if(Im(i+2,j+2)>(Im(i,j)+T))
                A(7,1)=1;
            elseif(Im(i+2,j+2)<(Im(i,j)-T))
                A(7,1)=-1;
            else
                A(7,1)=0;
            end
  
% 8 i+3,j+1
            if(Im(i+3,j+1)>(Im(i,j)+T))
                A(8,1)=1;
            elseif(Im(i+3,j+1)<(Im(i,j)-T))
                A(8,1)=-1;
            else
                A(8,1)=0;
            end
  
% 9 i+3,j
            if(Im(i+3,j)>(Im(i,j)+T))
                A(9,1)=1;
            elseif(Im(i+3,j)<(Im(i,j)-T))
                A(9,1)=-1;
            else
                A(9,1)=0;
            end
  
% 10 i+3,j-1
            if(Im(i+3,j-1)>(Im(i,j)+T))
                A(10,1)=1;
            elseif(Im(i+3,j-1)<(Im(i,j)-T))
                A(10,1)=-1;
            else
                A(10,1)=0;
            end
% 11 i+2,j-2
            if(Im(i+2,j-2)>(Im(i,j)+T))
                A(11,1)=1;
            elseif(Im(i+2,j-2)<(Im(i,j)-T))
                A(11,1)=-1;
            else
                A(11,1)=0;
            end
% 12 i+1,j-3
            if(Im(i+1,j-3)>(Im(i,j)+T))
                A(12,1)=1;
            elseif(Im(i+1,j-3)<(Im(i,j)-T))
                A(12,1)=-1;
            else
                A(12,1)=0;
            end
% 13 i,j-3
            if(Im(i,j-3)>(Im(i,j)+T))
                A(13,1)=1;
            elseif(Im(i,j-3)<(Im(i,j)-T))
                A(13,1)=-1;
            else
                A(13,1)=0;
            end
 % 14 i-1,j-3
            if(Im(i-1,j-3)>(Im(i,j)+T))
                A(14,1)=1;
            elseif(Im(i-1,j-3)<(Im(i,j)-T))
                A(14,1)=-1;
            else
                A(14,1)=0;
            end
% 15 i-2,j-2
            if(Im(i-2,j-2)>(Im(i,j)+T))
                A(15,1)=1;
            elseif(Im(i-2,j-2)<(Im(i,j)-T))
                A(15,1)=-1;
            else
                A(15,1)=0;
            end
% 16 i-3,j-1
            if(Im(i-3,j-1)>(Im(i,j)+T))
                A(16,1)=1;
            elseif(Im(i-3,j-1)<(Im(i,j)-T))
                A(16,1)=-1;
            else
                A(16,1)=0;
            end
       
        % All 0 means no pixel breaches threshold    
        if(nnz(A)~=0)
           flag = check(A',N);
           if (flag~=0)
                fastX(1,count)=j;
                fastY(1,count)=i;
                count=count+1;
           end
        end
       
        end
end
figure('Name','Q4: T=10','NumberTitle','off'); imagesc(Im); 
colormap('gray');
axis('image');
hold on; 
for i=1:size(fastX,2)
    plot(fastX(1,i),fastY(1,i),'r.'); 
end
hold off;
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 5
clc;
clear all;
% Load image 
Im = double(imread('given_pics/tower.jpg'));
% Define radius, no. of consecutive points, intensity thresholds
N = 9;
Threshold = [20,30,50];
nr = size(Im,1)-4;
nc = size(Im,2)-4;
% For each pixel, check 16 pixel neighbourhood(1=above,-1=below,0=pixel
% doesn't exist or neither above nor below
for itr=1:size(Threshold,2)
    T = Threshold(1,itr);
    fastX=0;
    fastY=0;
    count = 1;
    for i=4:nr
        for j=4:nc
% 1 i-3,j
            if(Im(i-3,j)>(Im(i,j)+T))
                A(1,1)=1;
            elseif(Im(i-3,j)<(Im(i,j)-T))
                A(1,1)=-1;
            else
                A(1,1)=0;
            end
        
% 2 i-3,j+1
            if(Im(i-3,j+1)>(Im(i,j)+T))
                A(2,1)=1;
            elseif(Im(i-3,j+1)<(Im(i,j)-T))
                A(2,1)=-1;
            else
                A(2,1)=0;
            end

% 3 i-2,j+2
            if(Im(i-2,j+2)>(Im(i,j)+T))
                A(3,1)=1;
            elseif(Im(i-2,j+2)<(Im(i,j)-T))
                A(3,1)=-1;
            else
                A(3,1)=0;
            end
% 4 i-1,j+3
            if(Im(i-1,j+3)>(Im(i,j)+T))
                A(4,1)=1;
            elseif(Im(i-1,j+3)<(Im(i,j)-T))
                A(4,1)=-1;
            else
                A(4,1)=0;
            end
% 5 i,j+3
            if(Im(i,j+3)>(Im(i,j)+T))
                A(5,1)=1;
            elseif(Im(i,j+3)<(Im(i,j)-T))
                A(5,1)=-1;
            else
                A(5,1)=0;
            end
% 6 i+1,j+3
            if(Im(i+1,j+3)>(Im(i,j)+T))
                A(6,1)=1;
            elseif(Im(i+1,j+3)<(Im(i,j)-T))
                A(6,1)=-1;
            else
                A(6,1)=0;
            end
%  7 i+2,j+2
            if(Im(i+2,j+2)>(Im(i,j)+T))
                A(7,1)=1;
            elseif(Im(i+2,j+2)<(Im(i,j)-T))
                A(7,1)=-1;
            else
                A(7,1)=0;
            end
  
% 8 i+3,j+1
            if(Im(i+3,j+1)>(Im(i,j)+T))
                A(8,1)=1;
            elseif(Im(i+3,j+1)<(Im(i,j)-T))
                A(8,1)=-1;
            else
                A(8,1)=0;
            end
  
% 9 i+3,j
            if(Im(i+3,j)>(Im(i,j)+T))
                A(9,1)=1;
            elseif(Im(i+3,j)<(Im(i,j)-T))
                A(9,1)=-1;
            else
                A(9,1)=0;
            end
  
% 10 i+3,j-1
            if(Im(i+3,j-1)>(Im(i,j)+T))
                A(10,1)=1;
            elseif(Im(i+3,j-1)<(Im(i,j)-T))
                A(10,1)=-1;
            else
                A(10,1)=0;
            end
% 11 i+2,j-2
            if(Im(i+2,j-2)>(Im(i,j)+T))
                A(11,1)=1;
            elseif(Im(i+2,j-2)<(Im(i,j)-T))
                A(11,1)=-1;
            else
                A(11,1)=0;
            end
% 12 i+1,j-3
            if(Im(i+1,j-3)>(Im(i,j)+T))
                A(12,1)=1;
            elseif(Im(i+1,j-3)<(Im(i,j)-T))
                A(12,1)=-1;
            else
                A(12,1)=0;
            end
% 13 i,j-3
            if(Im(i,j-3)>(Im(i,j)+T))
                A(13,1)=1;
            elseif(Im(i,j-3)<(Im(i,j)-T))
                A(13,1)=-1;
            else
                A(13,1)=0;
            end
 % 14 i-1,j-3
            if(Im(i-1,j-3)>(Im(i,j)+T))
                A(14,1)=1;
            elseif(Im(i-1,j-3)<(Im(i,j)-T))
                A(14,1)=-1;
            else
                A(14,1)=0;
            end
% 15 i-2,j-2
            if(Im(i-2,j-2)>(Im(i,j)+T))
                A(15,1)=1;
            elseif(Im(i-2,j-2)<(Im(i,j)-T))
                A(15,1)=-1;
            else
                A(15,1)=0;
            end
% 16 i-3,j-1
            if(Im(i-3,j-1)>(Im(i,j)+T))
                A(16,1)=1;
            elseif(Im(i-3,j-1)<(Im(i,j)-T))
                A(16,1)=-1;
            else
                A(16,1)=0;
            end
       
        % All 0 means no pixel breaches threshold    
        if(nnz(A)~=0)
           flag = check(A',N);
           if (flag~=0)
                fastX(1,count)=j;
                fastY(1,count)=i;
                count=count+1;
           end
        end
       
        end
    end
    figure('Name',strcat('Q4: T=',num2str(T)),'NumberTitle','off'); 
    imagesc(Im); 
    colormap('gray');
    axis('image');
    hold on; 
    for i=1:size(fastX,2)
        plot(fastX(1,i),fastY(1,i),'r.'); 
    end
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%