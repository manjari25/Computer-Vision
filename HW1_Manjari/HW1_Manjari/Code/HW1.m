% Manjari Akella
% CSE5524 - HW1
% 09/02/2013

mkdir('Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 1

refresh();                                                       % clear workspace, command window and close all windows
grayIm = imread('given_pics/buckeyes_gray.bmp');                 
figure('Name','Question 1','NumberTitle','off'), subplot(121), imagesc(grayIm);      
axis('image');                                                   % Set axis to 'image' type(upper right corner is {0,0})
colormap('gray');                                                % Sets color map to gray
imwrite(grayIm, 'Output/Q1_buckeyes_gray.jpg');                  % Write back to disk in jpg fromat
pause;
rgbIm = imread('given_pics/buckeyes_rgb.bmp');
subplot(122),imagesc(rgbIm);
axis('image');
imwrite(rgbIm, 'Output/Q1_buckeyes_rgb.jpg');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2

refresh();
rgbIm = imread('Output/Q1_buckeyes_rgb.jpg');                                          % Read rgb image
grayIm = rgb2gray(rgbIm);                                                              % Apply conversion formula
figure('Name','Question 2','NumberTitle','off'), subplot(121),imshow(rgbIm);           % Show rgb image
subplot(122),imshow(grayIm);                                                           % Show converted grayscale image
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3

refresh();
zBlock = zeros(10,10);                                    % block of black pixel 
oBlock = ones(10,10)*255;                                 % block of white pixels
pattern = [zBlock oBlock; oBlock zBlock];                 % pattern of black and white blocks
checkerIm = repmat(pattern, 5, 5);                        % chess board by repeating the pattern
imwrite(uint8(checkerIm), 'Output/Q3_checkerIm.bmp');        % write file to disk
Im = imread('Output/Q3_checkerIm.bmp');                      
figure('Name','Question 3','NumberTitle','off'), imagesc(Im); % Display in grayscale
colormap('gray')
axis('image');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4

% face.bmp given on the website

refresh();
faceIm=double(imread('given_pics/face.bmp'));
i=1;
figure('Name','Question 4: Sigma values Variation','NumberTitle','off');
for sigma=20:-1.5:0.5
    G = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
    gIm = imfilter(faceIm, G, 'replicate');
    
% Follow each transition of sigma
%   sigma  
%   colormap('gray');
%   imagesc(gIm);
%   pause(0.5);

% Display in single window

    subplot(4,4,i), imagesc(gIm);
    colormap('gray');
    imagesc(gIm);
    i=i+1;
    
end
pause;

% Identify the movie poster !
 
refresh();
pics = dir('my_pictures/*.jpg'); 
figure('Name','Question 4: Sigma values Variation(My_pictures)','NumberTitle','off');
for i = 1:size(pics,1) % loop for each picture in directory
    Im = double(imread(strcat('my_pictures/', pics(i).name)));
        for sigma=20:-1.5:0.5  %loop for each sigma
            G = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
            gIm =(imfilter(Im, G, 'replicate'));
            imshow(uint8(gIm));
            pause(0.5);          
        end
end

%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
