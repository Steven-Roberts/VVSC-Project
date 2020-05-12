% The first value comes from HW5 and the others are made up
alphas = [1, 1.1, 1.35, 1.42, 1.59, 1.7, 1.88, 1.9, 1.92, 2.1]';
avms = [0.0788, 0.0801, 0.0955, 0.0922, 0.100, 0.109, 0.132, 0.130, 0.148, 0.155]';
aPred = 2.5;

z = linspace(0.9, 2.7);

f = fit(alphas, avms, 'poly1');
pInt = predint(f, z);
scatter(alphas, avms);
hold on;
plot(z, feval(f, z));
plot(z, pInt, 'k--');
scatter(aPred, feval(f, aPred), 60, 'filled');
text(aPred, feval(f, aPred), 'Extrapolated AVM');
xlabel('\alpha');
ylabel('AVM');
legend({'AVM', 'Linear Regression', '95% prediction interval'}, 'Location', 'NorthWest');