clear;

timesteps = [8, 12, 16, 24, 32, 48, 64];
gridPts = timesteps;

for i = 1:length(timesteps)
    tstep = timesteps(i);
    pts = gridPts(i);
    
    sol(i) = solve_allen_cahn(tstep, pts);
end

h = [sol.h];
err1 = [sol.err1];
err2 = [sol.err2];
errInf = [sol.errInf];

figure;
loglog(h, err1, '-*');
hold on;
loglog(h, err2, '-+');
loglog(h, errInf, '-o');
loglog([2e-2, 4e-2], [2e-2, 8e-2], 'LineWidth', 2);

ylabel('Error');
xlabel('h = \Delta t = \Delta x');
legend({'1 norm', '2 norm', '\infty norm', 'Reference slope 2'}, 'Location', 'NorthWest');
hold off;

hSub = h(2:end);
order1 = diff(log(err1)) ./ diff(log(h));
order2 = diff(log(err2)) ./ diff(log(h));
orderInf = diff(log(errInf)) ./ diff(log(h));

figure;
semilogx(hSub, order1, '-*');
hold on;
semilogx(hSub, order2, '-+');
semilogx(hSub, orderInf, '-o');
xlabel('h = \Delta t = \Delta x');
ylabel('Order');
legend({'1 norm order', '2 norm order', '\infty norm order'}, 'Location', 'NorthWest');
