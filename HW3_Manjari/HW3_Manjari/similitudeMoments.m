function [ Nval ] = similitudeMoments( Im, i, j )
% Computes row and column averages of the given image and then the similitude moment for
% i,j and returns the value
n=0;
m00=sum(sum(Im));
m10=0;
m01=0;
e=((i+j)/2)+1;
% Calculate row and column averages
for r=1:size(Im,1)
    for c=1:size(Im,2)
        m10 = m10 + r*Im(r,c);
        m01 = m01 + c*Im(r,c);
    end
end
% Row average and column average (Centriod)
rc = m10/m00; 
cc = m01/m00;
% Calculate the numerator of the moment
for r=1:size(Im,1)
    for c=1:size(Im,2)
       n = n + ((c-cc)^i)*((r-rc)^j)*Im(r,c);
    end
end
% Calculate the value of the moment by dividing by m00
Nval = n/(m00^e); 
end


