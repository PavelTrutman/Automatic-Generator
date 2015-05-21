% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
% 
% Benchmark templates for benchmarking the matrix partitioning.

function [ret] = bench_matrixPartitioning(cfg)
  
  % generate MATLAB code only
  cfg.exportCode = {'matlab'};
  % use systematic generator
  cfg.PolynomialsGenerator = 'systematic';
  cfg.PolynomialsGeneratorCfg.GJstep = 1;
  
  % no partitioning
  ret{1}.info = 'No matrix partitioning.';
  ret{1}.abbrev = 'none';
  ret{1}.cfg = cfg;
  ret{1}.cfg.matrixPartitioning = 'none';
  
  % partitioning of the last elimination
  ret{2}.info = 'Matrix partitioning of the last elimiantion.';
  ret{2}.abbrev = 'last';
  ret{2}.cfg = cfg;
  ret{2}.cfg.matrixPartitioning = 'last';
  
  % partitioning of the last elimination
  ret{3}.info = 'Matrix partitioning of all elimiantions.';
  ret{3}.abbrev = 'all';
  ret{3}.cfg = cfg;
  ret{3}.cfg.matrixPartitioning = 'all';
  
end
