function sol = solve_allen_cahn(timesteps, spatialGridPts, beta)

tspan = [0, 1];
x = linspace(0, 1, spatialGridPts);
[x, y] = meshgrid(x, x);
x = x(:);
y = y(:);
y0 = cos(2 * pi * x) + cos(4 * pi * y).^2;

params.alpha = 1;
params.beta = beta;
params.n = spatialGridPts;    
params.forcing = [];

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

end
