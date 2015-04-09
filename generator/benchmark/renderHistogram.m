% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
% 
% Default function for rendering errors for benchmarking.
%
% results - results of benchmark
% benchmark = settings of benchmark

function [] = renderHistogram(results, benchmark)
  
  % reshape into vector and get rid of NaNs
  vector = reshape(results.err, [], 1);
  vector = abs(vector(isnan(vector) == 0));
  
  % logarithm
  vector = log10(vector);
  
  % get rid of -Inf
  vector = vector(vector ~= -Inf);
  
  % prepare data for histogram
  range = (min(vector):(max(vector)-min(vector))/100:max(vector))';
  count = histc(vector, range);
  
  % plot histogram
  figure;
  f = fit(range, count, 'smoothingspline');
  plot(f, '-b');
  legend('off');
  axis([min(range) max(range) min(count) max(count)]);
  title(['Histogram of error. (', benchmark.info, ')']);
  ylabel('Frequency');
  xlabel('Log10 of error');
  drawnow;
  
end