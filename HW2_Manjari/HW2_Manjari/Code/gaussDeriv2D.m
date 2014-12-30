function [ Gx,Gy ] = gaussDeriv2D( sigma )
% Creates a 2D Gaussian Derivative mask in x and y directions
% Returns and displays them
range = ceil(3*sigma);
hsize = 2*ceil(3*sigma)+1;
for i = 1:hsize
    for j = 1:hsize
        x = i-range-1;
        y = j-range-1;
        t_1=x/(2*pi*(sigma^4));
        t_2=exp(-(x^2+y^2)/(2*sigma*sigma));
        Gx(j,i) = t_1*t_2;
    end
end

%generate a 2-D Gaussian kernel along y direction
Gy=Gx';

figure('Name','Q1: Gx (x-direction)','NumberTitle','off'),surf(Gx);
figure('Name','Q1: Gy (y-direction)','NumberTitle','off'),surf(Gy);
end