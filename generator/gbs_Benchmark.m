% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
% 
% problemName - name of the minimalProblem; function specifies the problem
%   must have the same name
% benchmark - function, which defines benchmark procedures
% inputData - set of data on which we will test the solvers
% correctOutput - correct solutions of the problem for the inputData. It is
%   used to compare the correctness of the generated solvers.
% validationFunction - function called to evaluate the results
% renderFunction - function to render results

function [results, inputData] = gbs_Benchmark(problemName, benchmarkFunction, inputData, correctOutput, validationFunction, renderFunction)
  
  % get minimal problem definition
  problem = str2func(problemName);
  [eq, known, unknown, kngroups, cfg, algB] = problem();
  
  % parse input parameters
  if nargin < 3
    % if no inputData, we have to generate random inputs
    inputData = generateInputData(eq, known, unknown, kngroups, cfg.benchmark.maxInputs);
  end
  if nargin < 4
    % default zero of polynomials
    validationFunction = @validateZeroPolynomials;
    correctOutput = [];
  end
  if nargin < 6
    % default histogram
    renderFunction = @renderHistogram;
  end
  
  % enable benchmark mode
  cfg.benchmark.enable = 1;
  
  % get benchmark definitions
  benchmarkConfig = benchmarkFunction(cfg);
  
  % do benchmarks
  fprintf('\nBENCHMARK STARTED\n');
  results = cell(length(benchmarkConfig), 1);
  for i = 1:length(benchmarkConfig)
    benchmark = benchmarkConfig{i};
    
    % print info
    fprintf('\n################################################################\n');
    fprintf(['\n', benchmark.info, '\n\n\n']);
    
    % create code
    gbs_CreateCode([problemName, '_', benchmark.abbrev], eq, known, unknown, kngroups, benchmark.cfg, algB);
    
    % update MATLAB cache to recognize new generated functions
    rehash;
    
    solver = str2func(['solver_', problemName, '_' ,benchmark.abbrev]);
    
    fprintf('\nSolver generated.\n');
    fprintf('Benchmarking solver.\n');
    
    % call solver for all inputData and save results
    results{i}.solution = cell(size(inputData, 1), 1);
    results{i}.benchData = cell(size(inputData, 1), 1);
    results{i}.time = zeros(size(inputData, 1), 1);
    reverseStr = '';
    for j = 1:min(size(inputData, 1), cfg.benchmark.maxInputs)
      try
        [results{i}.solution{j}, results{i}.time(j)] = solver(inputData{j});
        msg = sprintf('  %2.0f %%%% done', j/min(size(inputData, 1), cfg.benchmark.maxInputs)*100);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg) - 1);
      catch
        fprintf(reverseStr);
        reverseStr = '';
        fprintf('  solver executed on data %d failed!\n', j);
        results{i}.solution{j} = zeros(length(unknown), 0);
        results{i}.time(j) = NaN;
      end
    end
    fprintf(reverseStr);
    
    fprintf('Evaluating results.\n');
    % get errors of solutions
    results{i}.err = validationFunction(inputData, correctOutput, results{i}.solution, eq, unknown, known);
    
    % render results
    renderFunction(results{i}, benchmark, length(benchmarkConfig), i, false);
    
  end

end