function Yp_P = GetYp_P(Xp_t, Yp_t, predU, P, A, B, C)
YpCell = repmat({Yp_t}, P, 1);
Xp = cell2mat(Xp_t);
U = cell2mat(predU);
for i = 1 : P
    Xp = A * Xp + B * U(i);
    YpCell{i, 1} = C * Xp;
end
Yp_P = cell2mat(YpCell);
end