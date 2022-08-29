function ec = GetEc(Y_P_cell, Yr_P_cell, K)
ec = 0;
for k = 1 : K - 1
ec = ec + (Y_P_cell{:, k} - Yr_P_cell{:, k});
end
row_ec = size(ec, 1);
mat = zeros(row_ec * 2, row_ec); 
mat(1 : 2 : end, :) = zeros(row_ec);
mat(2 : 2 : end, :) = eye(row_ec);
ec = mat * ec;
end