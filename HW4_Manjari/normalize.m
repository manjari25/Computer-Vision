function [ I ] = normalize( Im )
%Normalize Im to [0-1] and return normalized image
for i=1:size(Im,1)
    for j=1:size(Im,2)
        I(i,j)  = (Im(i,j) - min(min(Im)))/(max(max(Im))-min(min(Im)));
    end
end

