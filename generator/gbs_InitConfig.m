% Config file.
% (GBsolver subroutine)
%
% by Martin Bujnak, mar 2008
% last edit by Pavel Trutman, oct 2014

function [cfg] = gbs_InitConfig()

    cfg.ordering = 3;   %grevlex - Do not modify!
    
    % prime field generator
    cfg.prime = 30097;

    % the total degree of the polynomials, we are generating polynomials up to, is inceased by GJstep after each GJ elimination
    % 0 means perform only one GJ elimination at the end
    cfg.GJstep = 0;
    
    %
    % list of equations which should be used in Groebner basis solver (cfg.GBSolver).
    % e.g. set cfg.GBrestrictEq = 1:4; to use equations 1 to 4 only.
    cfg.GBrestrictEq = [];

    %
    % groebner basis solver
    % any function of the form 
    %     "[algB res] = GBSolver(cfg, eq, known, unknown);"
    
    % maple solver
    cfg.GBSolver = @gbs_findAlgB_maple;         
    
    % Macaulay2 gb solver
    cfg.GBSolver = @gbs_findAlgB_macaulay;      
    
    % instance generators
    % any function of the form 
    %     "function [eqi] = gbs_RandomInstance(cfg, eq, known, unknown)"
    %
    % define your own function if you want to study a special scene configurations
    % like planar scenes etc.

    % random instance within Zp field
    cfg.InstanceGenerator = @gbs_RandomInstanceZp;   
    cfg.ZpGeneratorMaxNumber = 999;             

    % random instance (coeffs may be bigger than cfg.prime)
    % cfg.InstanceGenerator = @gbs_RandomInstance;    
    % cfg.ZpGeneratorMaxNumber = 55;             
    
    % polynomial removing step
    cfg.bdoReduce = true;
    cfg.bIncremental = false; 
    
    % crash recovery from a log file. 
    % copy&paste matlab console output to a log file and run solver again
    % while setting cfg.crashlog = 'consoleout.log'. Copy just solver output.
    cfg.crashlog = [];
    
    % export coefficients * monomials vectors (debug)
    cfg.exportEqs = false;
    cfg.bGeneratorSolverResult = false;
    
    % code to export
    cfg.exportCode = {'matlab' 'maple'};
end