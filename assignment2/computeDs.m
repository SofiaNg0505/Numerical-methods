function d = computeDs(h, v, alpha0, Tout)
    alpha = sqrt(v*v/4 + alpha0*alpha0) - v/2;
%     fprintf("alpha = %f\n", alpha);
    denom = -3 - 2*h*alpha;
    d = [-2*h*alpha/denom*Tout, -4/denom, 1/denom];
end