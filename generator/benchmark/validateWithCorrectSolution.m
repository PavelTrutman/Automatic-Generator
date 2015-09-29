% Prokop Silhavy, silhapro@fel.cvut.cz, August 2015
%
% Validation function used by benchmark utility.
% This function compares solutions computed by automatic generator with 
% correct solutions defined by user.
%
% correctOutput - solutions defined bz user
% solution - computed solutions by solver

function [err] = validateWithCorrectSolution(~, correctOutput, solution, ~, ~, ~)

    err = NaN(size(solution, 1), 0);
    
    reverseStr = '';

    solutionErr = [];
    for i = 1:length(solution)
        
        solutionErr = computeSolutionDifference(correctOutput{i}, solution{i});
        
        % enlarge matrix if needed
        if length(solutionErr) > size(err, 2)
            oldErr = err;
            err = NaN(size(solution, 1), length(solutionErr));
            err(:, 1:size(oldErr, 2)) = oldErr;
        end
        
        err(i,1:length(solutionErr)) = solutionErr;
        
        msg = sprintf('  %2.0f %%%% done', i/size(solution, 1)*100);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg) - 1);
    end

    fprintf(reverseStr);

end