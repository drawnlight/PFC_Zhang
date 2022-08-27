function Diag_A = GetDiagMarix(A, N)
A_cell_array = repmat({A}, 1, N);
Diag_A = blkdiag(A_cell_array{:});
end