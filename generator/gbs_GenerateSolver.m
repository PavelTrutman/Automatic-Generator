% This function generate solver for the given problem specified as
% function. Minimal problems are typicaly defined in folder
% 'minimalProblems'. 
% 
% Minimal problem is every function of the form 
%   '[eq, known, unknown, kngroups, cfg, algB] = nameOfMinimalProblem()'.
% 
% Pavel Trutman, pavel.trutman@fel.cvut.cz, March 2015
% 
% problemName - name of the minimalProblem; function specifies the problem
%   must have the same name

function [res, export] = gbs_GenerateSolver(problemName)
  
  % get minimal problem definition
  problem = str2func(problemName);
  [eq, known, unknown, kngroups, cfg, algB] = problem();
  
  %create code
  [res, export] = gbs_CreateCode(problemName, eq, known, unknown, kngroups, cfg, algB);

end