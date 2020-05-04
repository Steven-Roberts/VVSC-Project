rng(1234);
n = 64;

betaMean = 1.25;
betaStd = 0.5;

sol = solve_allen_cahn(n, n, betaMean - betaStd);
a = sol.srq;

sol = solve_allen_cahn(n, n, betaMean + betaStd);
b = sol.srq;

samples = 100;
betas = sort(lhsnorm(betaMean, betaStd^2, samples));
srqs = zeros(size(betas));
for s = 1:samples
    sol = solve_allen_cahn(n, n, betas(s));
    srqs(s) = sol.srq;
end

quadFit = polyfit(betas, srqs, 2);

t = linspace(betas(1) - 0.1, betas(end) + 0.1);

figure;
plot(betas, srqs, '--o', 'LineWidth', 2);
hold all;
plot(t, polyval(quadFit, t));
legend({'SRQ from samples', sprintf('%.3fx^2 + %.3fx + %.3f', quadFit)}, 'Location', 'NorthWest');
xlabel('Model parameter \beta');
ylabel('SRQ');
hold off;


% chi = [0.55, 0.95, 1.0, 1.1, 1.5];
chi = [0.1, 0.4, 0.6, 0.75, 0.8, 0.9, 0.91, 0.97, 1.3, 1.6];
synthData = a + chi * (b - a);

figure;
ecdf(synthData);
hold all;
ecdf(srqs);
xlabel('SRQ');
ylabel('Cumulative probibility');
legend({'Measurements', 'Model'}, 'Location', 'NorthWest');

amv(synthData, srqs)
