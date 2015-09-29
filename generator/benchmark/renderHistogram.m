% Pavel Trutman, pavel.trutman@fel.cvut.cz, April 2015
% 
% Default function for rendering errors for benchmarking.
%
% results - results of benchmark
% benchmark - settings of benchmark
% numFigures - total count of figures to be displayed on the screen
% i - index of current subfigure
% fitHistogram - boolean if graph should be plot or histogram 

function [] = renderHistogram(results, benchmark, numFigures, i, fitHistogram)
  
  % reshape into vector and get rid of NaNs
  vector = reshape(results.err, [], 1);
  vector = abs(vector(isnan(vector) == 0));
  
  % logarithm
  vector = log10(vector);
  
  % get rid of -Inf
  ix = vector ~= -Inf;
  vector = vector(ix);
  
  % prepare data for histogram
  range = (min(vector):(max(vector)-min(vector))/100:max(vector))';
  count = histc(vector, range);
  
  % divide screen into adequate number of subfigures
  if numFigures <= 4
    rows = 2;
    cols = 2;
  elseif numFigures <= 6
    rows = 2;
    cols = 3;
  else
    rows = 3;
    cols = 4;
  end
  
  % plot histogram
  subfig(rows, cols, i);
  if fitHistogram      
      f = fit(range, count, 'smoothingspline');
      plot(f, '-b');
  else
      histogram(vector,100);
  end
  
  legend('off');
  axis([min(range) max(range) min(count) max(count)]);
  title(['Histogram of error. (', benchmark.info, ')']);
  ylabel('Frequency');
  xlabel(sprintf('Log10 of |error| for |error| > 0, %.2e of error = 0', nnz(~ix)/length(ix)));
  drawnow;
  
end