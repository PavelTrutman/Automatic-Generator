function [ret] = bench_polynomialAddingMethods(cfg)
  
  % generate MATLAB code only
  cfg.exportCode = {'matlab'};
  
  % solver using primitive method to add polynomials
  ret{1}.info = 'Primitive solver.';
  ret{1}.abbrev = 'primitive';
  ret{1}.cfg = cfg;
  ret{1}.cfg.PolynomialsGenerator = 'Primitive';
  ret{1}.cfg.PolynomialsGeneratorCfg.GJstep = 0;
  
  % solver using methods of the F4 Algorithm to add polynomials
  ret{2}.info = 'F4 Algorithm solver.';
  ret{2}.abbrev = 'F4';
  ret{2}.cfg = cfg;
  ret{2}.cfg.PolynomialsGenerator = 'F4';
  ret{2}.cfg.PolynomialsGeneratorCfg.Sel = @F4_SelNormal;

end
