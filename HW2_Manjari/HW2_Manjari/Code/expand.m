function [ I ] = expand( E2 )
% Interpolate and expand (Laplaican Pyramid)
I = zeros(2*size(E2,1)-1,2*size(E2,2)-1,'double');
% Fill odd rows and columns in new image using old one
x=0;
for i=1:2:size(I,1)
y=0;
    for j=1:2:size(I,2)
        I(i,j)= E2(i-x,j-y);
        y=y+1;
    end
    x=x+1;
end

% Interpolate even columns
for i=1:size(I,1)
    for j=1:size(I,2)
        if ((I(i,j)==0)&&(mod(i,2)~=0))
            a = I(i,j-1);
            b = I(i,j+1);
            I(i,j)= (a+b)/2;
        end
    end
end

% Interpolate even rows
for i=1:size(I,1)
    for j=1:size(I,2)
        if (I(i,j)==0)
            a = I(i-1,j);
            b = I(i+1,j);
            I(i,j)= (a+b)/2;
        end
    end
end

end

