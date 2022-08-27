function Yr_P = GetYr_P(strat_idx, w, Yp_t, c_t, P)
tmpCell = repmat({Yp_t}, P, 1);
for i = 1 : P
    tmpCell{i, :} = GetYr(strat_idx, w, Yp_t, c_t);
    strat_idx = strat_idx + 1;
end
Yr_P = cell2mat(tmpCell);
end