% Config file.
% (GBsolver subroutine)
%
% by Martin Bujnak, mar 2008
% last edit by Pavel Trutman, February 2015

function [cfg] = gbs_InitConfig()

    cfg.ordering = 3;   %grevlex - Do not modify!
    
    % prime field generator
    cfg.prime = 30097;
    
    % list of equations which should be used in Groebner basis solver (cfg.GBSolver).
    % e.g. set cfg.GBrestrictEq = 1:4; to use equations 1 to 4 only.
    cfg.GBrestrictEq = [];

    % groebner basis solver
    % any function of the form 
    %     "[algB res] = GBSolver(cfg, eq, known, unknown);"
    
    % Macaulay2 gb solver
    cfg.GBSolver = @gbs_findAlgB_macaulay;         
    
    % Maple gb solver
    %cfg.GBSolver = @gbs_findAlgB_maple;
    
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
    %cfg.bdoReduce = true;
    %cfg.bIncremental = false;
    
    
    % generator of polynomials
    % generate polynomials to selected degree
    %cfg.PolynomialsGenerator = 'Primitive';
    % config of this algorithm:
      % the total degree of the polynomials, we are generating polynomials up to, is inceased by GJstep after each GJ elimination
      % 0 means perform only one GJ elimination at the end
      %cfg.PolynomialsGeneratorCfg.GJstep = 0;
      
    % use strategies from F4 algorithm
    cfg.PolynomialsGenerator = 'F4';
    % config of this algorithm
      % define selection strategy of the F4 algorithm
      cfg.PolynomialsGeneratorCfg.Sel = @F4_SelNormal;
    
    
    % use matrix partitioning (by PaToH)
    % how to set up this external library see the 'installation.txt', this library is not available for Windows
    % possible values
    %   'none' - no matrix partitioning is used
    %   'last' - only the last elimination is done by using partitioning
    %   'all'  - for all eliminations is used partitioning
    cfg.matrixPartitioning = 'none';
    
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