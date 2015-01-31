%[P, PRows, PCols] = PaToHMatrixPart(matrix, 2, 'RWU');
prime = 30097;

[PCols1, ~] = find(PCols == 1);
[PCols2, ~] = find(PCols == 2);
[PRows1, ~] = find(PRows == 1);
[PRows2, ~] = find(PRows == 2);

tmp = P;
tmp(PRows1, PCols1) = zeros(size(PRows1, 1), size(PCols1, 1));
tmp(PRows2, PCols2) = zeros(size(PRows2, 1), size(PCols2, 1));
[~, cols] = find(tmp);
BCols = unique(cols);

A11Cols = setdiff(PCols1, BCols);
A22Cols = setdiff(PCols2, BCols);

New = zeros(size(matrix));
New(1:size(PRows1, 1), 1:size(A11Cols, 1)) = matrix(PRows1, A11Cols);
New(end-size(PRows2, 1)+1:end, size(A11Cols, 1)+1:size(A11Cols, 1)+size(A22Cols, 1)) = matrix(PRows2, A22Cols);
New(1:size(PRows1, 1), end-size(BCols, 1)+1:end) = matrix(PRows1, BCols);
New(end-size(PRows2, 1)+1:end, end-size(BCols, 1)+1:end) = matrix(PRows2, BCols);

mat1 = New(1:size(PRows1, 1), :);
mat1 = gjzpsp(mat1, prime);
[lastRow1, ~] = find(mat1(:, size(A11Cols, 1)), 1, 'last');

mat2 = New(end-size(PRows2, 1)+1:end, :);
mat2 = gjzpsp(mat2, prime);
[lastRow2, ~] = find(mat2(:, size(A11Cols, 1)+size(A22Cols, 1)), 1, 'last');

Res = [mat1(1:lastRow1, :); mat2; mat1(lastRow1+1:end, :)];

Res(lastRow1+lastRow2+1:end, size(A11Cols, 1)+size(A22Cols, 1)+1:end) = gjzpsp(Res(lastRow1+lastRow2+1:end, size(A11Cols, 1)+size(A22Cols, 1)+1:end), prime);

firstRow = lastRow1+lastRow2+1;
lastRow = size(matrix, 1);
firstCol = size(A11Cols, 1)+size(A22Cols, 1)+1;
lastCol = size(matrix, 2);
previousRow = firstRow - 1;
for i = firstCol:lastCol
  [row, ~] = find(Res(firstRow:lastRow, i) == 1);
  row = row + firstRow - 1;
  if size(row, 1) == 1
    if row > previousRow
      previousRow = row;
      [rows, ~, val] = find(Res(1:firstRow-1, i));
      for j = 1:size(rows)
        Res(rows(j), :) = mod(-val(j)*Res(row, :) + Res(rows(j), :), prime);
      end
    end
  end
end

Fin = zeros(size(matrix));
Fin(:, A11Cols) = Res(:, 1:size(A11Cols, 1));
Fin(:, A22Cols) = Res(:, size(A11Cols, 1)+1:size(A11Cols, 1)+size(A22Cols, 1));
Fin(:, BCols) = Res(:, end-size(BCols, 1)+1:end);

spy(Fin);