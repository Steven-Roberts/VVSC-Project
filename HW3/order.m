clear;

timesteps = [8, 12, 16, 24, 32, 48, 64];
gridPts = timesteps;

for i = 1:length(timesteps)
    tstep = timesteps(i);
    pts = gridPts(i);
    
    sol(i) = solve_allen_cahn(tstep, pts);
end

h = [sol.h];
loglog(h, [sol.err1], '-*');
hold on;
loglog(h, [sol.err2], '-+');
loglog(h, [sol.errInf], '-o');
% loglog([h(1), h(end)], 

ylabel('Error');
xlabel('h');
legend('1 norm', '2 norm', '\infty norm');
hold off;