function [t, x, y] = getPoint(T, X, Y, h)
%%t = getPoint(X, Y, T, h)
%%returns the temperature at the point (3, 1)

i = 3/h + 1;
j = 1/h + 1;

t = T(i, j); 
x = X(i, j);
y = Y(i, j);
end