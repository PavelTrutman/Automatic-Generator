% GJ elimination with PaToH partitioning
% Pavel Trutman, pavel.trutman@fel.cvut.cz, February 2015

function [res, workflow] = myPaToHNew(matrix, amCols, prime)

  %separete algB collumns
  amCols = sort(amCols);
  noAmCols = setdiff(1:size(matrix, 2), amCols);
  algBcnt = size(amCols, 2);
  redMatrix = matrix(:, noAmCols);
  
  %use PaToH to patition the matrix
  [P, PRows, PCols] = PaToHMatrixPart(redMatrix, 2, 'RWS');
  
  %extract rows and cols for each part
  [PRows1, ~] = find(PRows == 1);
  [PRows2, ~] = find(PRows == 2);
  [PCols1, ~] = find(PCols == 1);
  [PCols2, ~] = find(PCols == 2);
  
  tmp = P;
  tmp(PRows1, PCols1) = zeros(size(PRows1, 1), size(PCols1, 1));
  tmp(PRows2, PCols2) = zeros(size(PRows2, 1), size(PCols2, 1));
  [~, cols] = find(tmp);
  BCols = unique(cols);
  
  ACols1 = setdiff(PCols1, BCols);
  ACols2 = setdiff(PCols2, BCols);
  
  %elimination of the first matrix
  mat1 = [redMatrix(PRows1, [ACols1; BCols]), matrix(PRows1, amCols)];
  mat1 = gjzpsp(mat1, prime);
  if size(ACols1, 1) ~= 0
    [lastRow1, ~] = find(mat1(:, size(ACols1, 1)), 1, 'last');
  else
    lastRow1 = 0;
  end
  
  %elimination of the second matrix
  mat2 = [redMatrix(PRows2, [ACols2; BCols]), matrix(PRows2, amCols)];
  mat2 = gjzpsp(mat2, prime);
  if size(ACols2, 1) ~= 0
    [lastRow2, ~] = find(mat2(:, size(ACols2, 1)), 1, 'last');
  else
    lastRow2 = 0;
  end
  
  %assemble the top of the final matrix
  tmp = zeros(lastRow1 + lastRow2, size(matrix, 2) - algBcnt);
  tmp(1:lastRow1, [ACols1; BCols]) = mat1(1:lastRow1, 1:end-algBcnt);
  tmp(lastRow1+1:lastRow1 + lastRow2, [ACols2; BCols]) = mat2(1:lastRow2, 1:end-algBcnt);
  
  topRows = 1:size(tmp, 1);
  bottomRows = size(tmp, 1)+1:size(matrix, 1);
  
  res = zeros(size(matrix));
  res(topRows, noAmCols) = tmp;
  res(topRows, amCols) = [mat1(1:lastRow1, end - algBcnt + 1:end); mat2(1:lastRow2, end - algBcnt + 1:end)];
  
  %assemble the bottom of the final matrix
  tmp = zeros(size(bottomRows, 2), size(matrix, 2) - algBcnt);
  tmp(1:size(mat1, 1) - lastRow1, [ACols1; BCols]) = mat1(lastRow1 + 1:end, 1:end-algBcnt);
  tmp(size(mat1, 1) - lastRow1 + 1:end, [ACols2; BCols]) = mat2(lastRow2 + 1:end, 1:end-algBcnt);
  res(bottomRows, noAmCols) = tmp;
  res(bottomRows, amCols) = [mat1(lastRow1 + 1:end, end - algBcnt + 1:end); mat2(lastRow2 + 1:end, end - algBcnt + 1:end)];
  
  %again eliminate the bottom part
  nonzero = find(sum(res(bottomRows, :), 1) ~= 0);
  res(bottomRows, nonzero) = gjzpsp(res(bottomRows, nonzero), prime);
  
  %save workflow
  workflow.enable = 1;
  workflow.res = res;
  workflow.amCols = amCols;
  workflow.noAmCols = noAmCols;
  workflow.PRows1 = PRows1;
  workflow.PRows2 = PRows2;
  workflow.ACols1 = ACols1;
  workflow.ACols2 = ACols2;
  workflow.BCols = BCols;
  workflow.mat1TopRows = 1:lastRow1;
  workflow.mat2TopRows = 1:lastRow2;
  workflow.mat1BottRows = lastRow1 + 1:size(mat1, 1);
  workflow.mat2BottRows = lastRow2 + 1:size(mat2, 1);
  workflow.resMat1TopRows = 1:lastRow1;
  workflow.resMat2TopRows = (1:lastRow2) + lastRow1;
  lastRow = lastRow1 + lastRow2;
  if size(workflow.mat1BottRows, 2) == 0
    workflow.resMat1BottRows = 1:0;
  else
    workflow.resMat1BottRows = workflow.mat1BottRows - workflow.mat1BottRows(1) + lastRow + 1;
    lastRow = lastRow + size(workflow.resMat1BottRows, 2);
  end
  if size(workflow.mat2BottRows, 2) == 0
    workflow.resMat2BottRows = 1:0;
  else
    workflow.resMat2BottRows = workflow.mat2BottRows - workflow.mat2BottRows(1) + lastRow + 1;
  end
  workflow.topRows = topRows;
  workflow.bottomRows = bottomRows;
  workflow.resBottNonzeroCols = nonzero;
  
end