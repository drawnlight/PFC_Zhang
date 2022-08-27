function [theta, Phi, Kec, theta_u, H, G, X] = GetCoefMat(Ae, Be, Ce, P, N)
% vector programming version
% row_Ae = size(Ae, 1);
% col_Be = size(Be, 2);
% 
% base_cell = cellfun(@(cell, idx) mpower(cell, idx-1), ...
%     repmat({Ae}, 1, P), num2cell(1 : P), 'uniformoutput', false).';
% base_mat = cell2mat(base_cell);
% 
% 
% % theta matrix
% theta = base_mat * Ae;
% 
% tmp_cell = repmat({base_mat}, 1, P);
% trans_cell = cellfun(@(data, idx, row) ...
%     cell2mat({[repmat(zeros(row), idx - 1, 1);(data(1 : end - (idx - 1) * row, :))]}), ...
%     tmp_cell, num2cell(1 : P), repmat({row_Ae}, 1, P), 'uniformoutput', false);
% 
% % psi matrix
% Phi = cell2mat(cellfun(@(data, factor) data * factor, ...
%     trans_cell, repmat({Be}, 1, P),  'uniformoutput', false));
% 
% % xi matrix
% Kec = cell2mat(cellfun(@(data, factor) data * factor, ...
%     trans_cell, repmat({Ce}, 1, P),  'uniformoutput', false));
% 
% % theta_u matrix
% theta_u = Phi(:, 1 : col_Be);
% 
% % H matrix
% H = cell2mat(arrayfun(@(base, index) baseFunc(base, index), ...
%     0 : P -1, N * ones(1, P),  'uniformoutput', false)');
% 
% % G matrix
% firstRow = baseFunc(0, N);
% otherRow = cell2mat(arrayfun(@(base, index) baseFunc(base, index) - baseFunc(base - 1, index), ...
%     1 : P -1, N * ones(1, P - 1),  'uniformoutput', false)');
% G = [firstRow; otherRow];
% 
% % chi matrix
% X = Phi * G;

% % forloop version
[row_Ae, col_Ae] = size(Ae);
baseMat = eye(row_Ae * P, col_Ae);
for ii = 2 : P
    baseMat((ii - 1) * row_Ae + 1: ii * row_Ae, :) ...
        = baseMat((ii - 2) * row_Ae + 1 : (ii - 1) * row_Ae, :) * Ae;
end
% matrix theta
theta = baseMat * Ae;
% matrix Phi
[row_Be, col_Be] = size(Be);
Phi = zeros(row_Be * P, col_Be * P);
for ii = 1 : P
    Phi((ii - 1) * row_Be + 1 : row_Be * P, (ii - 1) * col_Be + 1 : ii * col_Be) ...
        = baseMat(1 : (P - ii + 1) * row_Be, :) * Be;
end
% matrix iterKec
[row_Ce, col_Ce] = size(Ce);
Kec = zeros(row_Ce * P, col_Ce * P);
for ii = 1 : P
    Kec((ii - 1) * row_Ce + 1 : end, (ii - 1) * col_Ce + 1 : ii * col_Ce) ...
        = baseMat(1 : (P - ii + 1) * row_Ce, :) * Ce;
end
% matrix theta_u
theta_u = Phi(:, 1 : col_Be);
% matrix H
H = zeros(P, N);
for ii = 0 : (P - 1)
    H(ii + 1, :) = baseFunc(ii, N);
end
% matrix G
G = zeros(P, N);
G(1, :) = baseFunc(0, N);
for ii = 1 : (P - 1)
    G(ii + 1, :) = baseFunc(ii, N) - baseFunc(ii - 1, N);
end
% X matrix
X = Phi * G;
end