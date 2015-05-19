% Pavel Trutman, pavel.trutman@fel.cvut.cz, May 2015
% 
% Benchmark templates for benchmarking methods of generating polynomials.

function [ret] = bench_polynomialGenerator(cfg)
  
  % generate MATLAB code only
  cfg.exportCode = {'matlab'};
  
  % solver using systematic method to add polynomials
  ret{1}.info = 'Systematic solver.';
  ret{1}.abbrev = 'systematic';
  ret{1}.cfg = cfg;
  ret{1}.cfg.PolynomialsGenerator = 'systematic';
  ret{1}.cfg.PolynomialsGeneratorCfg.GJstep = 0;
  
  % solver using methods of the F4 Algorithm to add polynomials
  ret{2}.info = 'F4 Algorithm solver without matrix partitioning.';
  ret{2}.abbrev = 'F4_none';
  ret{2}.cfg = cfg;
  ret{2}.cfg.PolynomialsGenerator = 'F4';
  ret{2}.cfg.PolynomialsGeneratorCfg.Sel = @F4_SelNormal;
  
  % solver using methods of the F4 Algorithm to add polynomials
  ret{3}.info = 'F4 Algorithm solver with matrix partitioning.';
  ret{3}.abbrev = 'F4_last';
  ret{3}.cfg = cfg;
  ret{3}.cfg.PolynomialsGenerator = 'F4';
  ret{3}.cfg.PolynomialsGeneratorCfg.Sel = @F4_SelNormal;
  ret{3}.cfg.matrixPartitioning = 'last';

end
