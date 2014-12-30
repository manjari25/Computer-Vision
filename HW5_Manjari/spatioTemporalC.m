function [ k ] = spatioTemporalC(x,y)
% Computes and returns the spatio temporal curvature matrix
difft1 = 1;
difft2 = 0;
    % Compute x'
    for i=1:(size(x,1)-1)
        diffx1(i,1) = x(i+1)-x(i);
    end
    % Compute x''
    for i=1:(size(x,1)-2)
        diffx2(i,1) = diffx1(i+1)-diffx1(i);
    end
    % Compute y'
    for i=1:(size(y,1)-1)
        diffy1(i,1) = y(i+1)-y(i);
    end
    % Compute y''
    for i=1:(size(y,1)-2)
        diffy2(i,1) = diffy1(i+1)-diffy1(i);
    end
    % Compute curvature    
    for i=1:(size(y,1)-2)
        m1=[diffy1(i,1),difft1;diffy2(i,1),difft2];
        m2=[diffx1(i,1),difft1;diffx2(i,1),difft2];
        m3=[diffx1(i,1),diffy1(i,1);diffx2(i,1),diffy2(i,1)];
        den = (diffx1(i,1)^2+diffy1(i,1)^2+difft1^2)^(3/2);
        num = sqrt(det(m1)^2+det(m2)^2+det(m3)^2);
        k(i,1)=num/den;
    end
end
