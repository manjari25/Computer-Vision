function [ flag ] = check( A,N)
% Checks to see if largest contiguous pixel patch under consideration is >= N or
% not. Returns flag = 0 if not else returns flag = 1
% A is a row vector, N is the threshold for number of contiguous pixels
% Values in A are : -1 = pixel is below threshold(T)
%                   +1 = pixel is above threshold(T)
%                   0 = either not a valid pixel or not above/below
%                   threshold
    countA = 0;
    countB = 0;
    % Break matrix so that circularity problem can be avoided
    ind = find(A==0,1,'first');
    % Check if ind is empty (i.e-no element is 0)
    if (isempty(ind))
        ind=0;
        % Find ind for which A(i),A(i+1) are different        
        for i=(1:size(A,2)-1)
             if(A(i)~=A(i+1))
                ind=i+1;
                break;
             else
                continue;
             end
        end 
       % If ind still 0, this indicates, entire vector has same value
       if(ind==0)
            % if entire vector is 1(above)
            if(A(1,1)==1)
                countA=size(A,2);
            % else if entire vector is -1(below)
            elseif(A(1,1)==-1)
                countB=size(A,2);
            end
       else
            temp = A(1:(ind-1));
            newA = [A(ind:end),temp];
            i1 = find(newA==1,1,'first');
            i2 = find(newA==-1,1,'first');
            for i=i1:size(newA,2)
                if(newA(1,i)==0)
                    break;
                elseif(newA(1,i)==-1) 
                    break;
                else
                    countA = countA+1;
                end
            end
            for i =i2:size(newA,2)
                if(newA(1,i)==0)
                    break;
                elseif(newA(1,i)==1) 
                    break;
                else
                    countB = countB+1;
                end
            end
       end            
    else
        % if ind is not empty, break vector 
        temp = A(1:(ind-1));
        newA = [A(ind:end),temp];
        i1 = find(newA==1,1,'first');
        i2 = find(newA==-1,1,'first');
        for i=i1:size(newA,2)
            if(newA(1,i)==0)
                break;
            elseif(newA(1,i)==-1) 
                break;
            else
                countA = countA+1;
            end
        end
        for i =i2:size(newA,2)
            if(newA(1,i)==0)
                break;
            elseif(newA(1,i)==1) 
                break;
            else
                countB = countB+1;
            end
        end
    end
    
    % Compare to N and set flag
    if ((countA>=N)||(countB>=N))
        flag=1;
    else
        flag=0;
    end   
end
