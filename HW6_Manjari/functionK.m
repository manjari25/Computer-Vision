function [ k_r ] = functionK( r )
% The Epanechnikov kernel with r(normalized radius) as input 
% Returns k(r)
    if (r<=1)
        k_r = 1-r;
    else
        k_r = 0;
    end
end

