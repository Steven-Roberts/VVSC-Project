x = [398 292 352 575 568 450 550 408 484 350 503 600 600]';
y = [0.15 0.05 0.23 0.43 0.23 0.4 0.44 0.44 0.45 0.09 0.59 0.63 0.6]';

z = linspace(min(x), max(x));

f = fit(x, y, 'poly1');
pInt = predint(f, z);
plot(f, x, y);
hold on;
plot(z, pInt);
