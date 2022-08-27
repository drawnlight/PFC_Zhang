% base function
function res = baseFunc(base, N)
res =[1, arrayfun(@(base, index) (base ^ index), base * ones(1, N - 1), 1 : N - 1)];
end