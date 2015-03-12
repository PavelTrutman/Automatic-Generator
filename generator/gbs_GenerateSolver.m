% This function generate solver for the given problem specified as
% function. Minimal problems are typicaly defined in folder
% 'minimalProblems'. 
% 
% Minimal problem is every function of the form 
%   '[codename, eq, known, unknown, kngroups, cfg, algB] = nameOfMinimalProblem()'.
% 
% Pavel Trutman, pavel.trutman@fel.cvut.cz, March 2015
% 
% problemName - name of the minimalProblem; function specifies the problem
%   must have the same name

function [res, export] = gbs_GenerateSolver(problemName)
  
  problem = str2func(problemName);
  [codename, eq, known, unknown, kngroups, cfg, algB] = problem();
  [res, export] = gbs_CreateCode(codename, eq, known, unknown, kngroups, cfg, algB);

end