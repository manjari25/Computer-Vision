function [ I2 ] = reduce_G( I1,w1,w2 )
%Reduction to lower levels (Gaussian Pyramid)
I = imfilter(I1,w1,'replicate');
I2full = imfilter(I,w2,'replicate');
I2 = I2full(1:2:size(I2full,1),1:2:size(I2full,2));
end

