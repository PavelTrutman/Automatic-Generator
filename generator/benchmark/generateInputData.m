% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
%
% Function generates inputData to be used by benchmark. Data are generated
% by Normal distribution with mean 0 and standard deviation 1.
%
% known - list of knowns
% kngroups - groups of knowns
% maxInputs - number of generated inputs

function [inputData] = generateInputData(~, known, ~, kngroups, maxInputs)
  
  % initialize variables
  inputData = cell(maxInputs, 1);
  
  % sort knowns into kngroups
  if isempty(kngroups)
    numKnowns = ones(length(known), 1);
  else
    numKnowns = zeros(0, 1);
    for i = 1:length(known)
      if size(numKnowns, 1) < kngroups(i)
        numKnowns(kngroups(i), 1) = 1;
      else
        numKnowns(kngroups(i)) = numKnowns(kngroups(i)) + 1;
      end
    end
  end
  
  % get dimensions
  rows = length(numKnowns);
  cols = max(numKnowns);
  
  % generate sets of input data
  for i = 1:maxInputs

    % for each known generate random number
   inputData{i, 1} = random('Normal', 0, 1, rows, cols);
   
  end
  
end