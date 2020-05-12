function sol = solve_allen_cahn(timesteps, spatialGridPts, alpha, beta, manufacture)

if nargin < 3
    alpha = 1;
    beta = 1;
end

if nargin < 5
    manufacture = false;
end

exactSol = @(t, x, y) exp(-t) * (cos(2 * pi * x) + cos(4 * pi * y).^2);

tspan = [0, 1];
x = linspace(0, 1, spatialGridPts);
[x, y] = meshgrid(x, x);
x = x(:);
y = y(:);
y0 = exactSol(tspan(1), x, y);

params.alpha = alpha;
params.beta = max(beta, eps);
params.n = spatialGridPts;
if manufacture
    params.forcing = @(t, x, y) exp(-t) * ((4 * pi^2 - 2) * cos(2 * pi * x) + (32 * pi^2 - 1) * cos(8 * pi * y) - 1) + ...
        exp(-3 * t) * (cos(2 * pi * x) + cos(4 * pi * x).^2).^3;
else
    params.forcing = [];
end

p = otp.allencahn.AllenCahnProblem(tspan, y0, params);

h = diff(p.TimeSpan) / timesteps;

% Use fixed timestepping
opts = MATLODE_OPTIONS('Hmin', h, 'Hmax', h, 'Hstart', h, 'AbsTol', inf, 'RelTol', inf, 'Max_no_steps', inf, ...
    'Jacobian', p.Rhs.Jacobian, 'Method', 'Ros-2', 'MatrixFree', true, 'storeCheckpoint', false);

[~, traj, stats] = MATLODE_ROS_FWD_Integrator(p.Rhs.F, p.TimeSpan, y0, opts);

if stats.ISTATUS.Nstp ~= timesteps
    warning('%d timesteps requested but %d taken', timesteps, stats.ISTATUS.Nstp);
end

sol.h = diff(tspan) / timesteps;
sol.dx = 1 / (spatialGridPts - 1);
sol.y = traj(end, :).';
sol.srq = mean(sol.y);
if manufacture
    sol.yExact = exactSol(tspan(end), x, y);
    sol.err = sol.y - sol.yExact;
    sol.err1 = norm(sol.err, 1) / p.NumVars;
    sol.err2 = rms(sol.err);
    sol.errInf = norm(sol.err, inf);
end

end
