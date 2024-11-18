function S = triArea(X1,X2,X3)
x1 = X1(1);
y1 = X1(2);
x2 = X2(1);
y2 = X2(2);
x3 = X3(1);
y3 = X3(2);
S = 0.5*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
end