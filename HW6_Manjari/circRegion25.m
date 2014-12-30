function [ info ] = circRegion25( I,center )
% Takes an image and center coordinates of terget patch as input
% and returns a matrix info having (c,r,R,G,B)values of pixels in 
% a 25 radius neighbourhood
n=1;
for i=1:size(I,1)
    for j=1:size(I,1)
        if(sqrt((i-center(1,2))^2+(j-center(1,1))^2)<=25)
            info(1,n) = j;
            info(2,n) = i;
            info(3,n) = I(i,j,1);
            info(4,n) = I(i,j,2);
            info(5,n) = I(i,j,3);
            n=n+1;
        end    
    end
end

end

