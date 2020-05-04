function area = amv(d1, d2)

[f1, x1] = ecdf(d1);
[f2, x2] = ecdf(d2);

x1(1) = x1(1) - eps;
x2(1) = x2(1) - eps;

[xMin, xMax] = bounds([x1; x2]);
area = integral(@(a) abs(interp1(x1, f1, a, 'nearest', 'extrap') - interp1(x2, f2, a, 'nearest', 'extrap')), xMin, xMax);

end
