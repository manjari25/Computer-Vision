function [ score ] = compare( modelCovMatrix,candCovMatrix )
% Takes the model and candidate covariance matrices and returns
% a score value representing the similarity between the two
[~,Eval] = eig(modelCovMatrix,candCovMatrix);
for i=1:size(Eval,1)
    for j=1:size(Eval,2)
        if(Eval(i,j)==0)
            lEvalsq(i,j)=0;
        else
            lEvalsq(i,j) = (log(Eval(i,j)))^2;
        end
    end
end
lEvalsum = sum(sum(lEvalsq));
score = sqrt(lEvalsum);
end

