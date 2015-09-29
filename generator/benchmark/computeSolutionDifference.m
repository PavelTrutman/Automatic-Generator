% Prokop Silhavy, silhapro@fel.cvut.cz, August 2015
%
% Utility for validation for benchmark.
% This function computes sinus of solutions computed by automatic generator
% and correct solutions defined by user.
%
% correctSolution - solutions defined by user
% solution - computed solutions by solver
function [err] = computeSolutionDifference( correctSolution, solution )

err = NaN(1,size(solution,2));

bestIndex = 0;
actError = 0;
d = 0;

for i = 1:size(solution,2)
    for j = 1:size(correctSolution,2)
        %actError = sum(abs(solution(:,i) - correctSolution(:,j)))/length(solution(:,i));
        d = dot(solution(:,i), correctSolution(:,j))/(sqrt(sum(solution(:,i).^2))*sqrt(sum(correctSolution(:,j).^2)))
        if abs(d)>1
            d = 1*sign(d);
        end
        actError = sin(acos(d));
        if isnan(err(i)) || actError < err(i)
            err(i) = actError;
        end
    end
end
end

