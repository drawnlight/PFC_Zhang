function Yp_P = GetY_P(Xp_t, Yp_t, predU, P, A, B, C)
YpCell = repmat({Yp_t}, P, 1);
Xp = Xp_t;
U = predU;
for i = 1 : P
    Xp = A * Xp + B * U(i);
    YpCell{i, 1} = C * Xp;
end
Yp_P = cell2mat(YpCell);
end