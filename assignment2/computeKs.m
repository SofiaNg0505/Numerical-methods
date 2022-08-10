function k = computeKs(h, v)
    k = [-(1/h + v/2)/h, 2/h/h, (v/2 - 1/h)/h];
end