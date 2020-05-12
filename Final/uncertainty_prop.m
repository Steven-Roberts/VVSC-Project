%% Config
rng(1234);
n = 64;

alphaMin = 0.1;
alphaMax = 1;
betaMean = 1.25;
betaStd = 0.5;

Na = 25;
Ne = 50;
alphas = linspace(alphaMin, alphaMax, Ne + 1);

srqs = zeros(Ne + 1, Na);

%% Compute SRQs
for i = 1:length(alphas)
    alpha = alphas(i);
    betas = sort(lhsnorm(betaMean, betaStd^2, Na));
    parfor j = 1:length(betas)
        sol = solve_allen_cahn(n, n, alpha, betas(j));
        srqs(i, j) = sol.srq;
    end
end

%% Plot results
fig = figure;
ax = axes(fig);
hold(ax, 'all');

for i = 1:length(alphas)
    ecdf(ax, srqs(i, :));
    ax.Children(1).Color = 'g';
end

ecdf(ax, min(srqs));
ax.Children(1).LineWidth = 2;
ax.Children(1).Color = 'k';
ax.Children(1).LineStyle = '--';

ecdf(ax, max(srqs));
ax.Children(1).LineWidth = 2;
ax.Children(1).Color = 'k';
ax.Children(1).LineStyle = '--';

ecdf(ax, mean(srqs));
ax.Children(1).LineWidth = 2;
ax.Children(1).Color = 'b';

legend(ax.Children([1,2,4]), {'Probabilistic (uniform)', 'P-box envelope', 'Sample CDF'}, 'Location', 'NorthWest');

xlabel('SRQ');
ylabel('Cumulative probibility');

hold(ax, 'off');

%% Total uncertainty
fig = figure;
ax = axes(fig);
hold(ax, 'all');

ecdf(ax, mean(srqs));
ax.Children(1).Color = 'k';

ecdf(ax, min(srqs));
ax.Children(1).Color = 'r';

ecdf(ax, max(srqs));
ax.Children(1).Color = 'r';

modelFormUncertainty = 0.176;
ecdf(ax, min(srqs) - modelFormUncertainty);
ax.Children(1).Color = 'g';

ecdf(ax, max(srqs) + modelFormUncertainty);
ax.Children(1).Color = 'g';

numericalUncertainty = 1.01e-4;
ecdf(ax, min(srqs) - modelFormUncertainty - numericalUncertainty);
ax.Children(1).Color = 'b';

ecdf(ax, max(srqs) + modelFormUncertainty + numericalUncertainty);
ax.Children(1).Color = 'b';

xlabel('SRQ');
ylabel('Cumulative probibility');
legend(ax.Children(1:2:end), {'Numerical uncertainty', 'Model form uncertainty', 'Input uncertainty', 'Probabilistic (uniform)'}, 'Location', 'NorthWest');

hold(ax, 'off');