% get reference tracjectory
function yr = GetYr(idx, w, y, c)
yr = power(w, idx) * y + (1 - power(w, idx)) * c;
end

% function yr = GetYr(idx, t, w, Yp, c)
% yr = c(t + idx) - (w ^ idx) * (c(t) - Yp);
% end