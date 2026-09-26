%%
% Solves Statistical Moment Steering
%%

clear;
close all;

initialization;
params_simulation;

%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Function ("quadratic", or "cvar")
opt.objective = "quadratic";
%opt.objective = "cvar";

% Additional Constraints
opt.constraints.skew = 0; % 0 or 1
opt.constraints.skew_eps = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute Initial Reference
if ~opt.constraints.skew 
    refTraj = generateReference();

    if opt.objective == "cvar"
        refTraj.CVaR_tau = 0;
    end
else
    filename_ref = "refTraj_" + opt.objective + ".mat";
    fprintf("\nLoading "+filename_ref+" as initial reference...\n")
    load(filename_ref, "refTraj")
end

%% SCVX* Params
opt_params.tol_opt     = 1E-4;
opt_params.tol_feas    = 1E-6;

% 1 > eta0 > eta1 > eta2 > 0
opt_params.eta0        = 1; % [0,1]
opt_params.eta1        = 0.75; % [0,1]
opt_params.eta2        = 0.25; % [0,1]

opt_params.alpha1      = 2; % Decrease trust: alpha1 > 1 
opt_params.alpha2      = 2; % Increase trust: alpha2 > 1
opt_params.beta        = 1.5; % beta > 1
opt_params.gamma       = 0.9; % [0,1]

opt_params.r_init      = 0.1; % r > 0
opt_params.r_min       = 1E-10;
opt_params.r_max       = 1;

opt_params.w_init      = 1; % w > 0
opt_params.w_max       = 1E10;

%% Solve Statistical Moment Steering
flag_plot = 0;
[refTraj, opt_hist] = solve_SMS(refTraj, opt, opt_params, flag_plot);
    
%% Propagate in Nonlinear and Update Moments
refTraj = propagateSigmaPoints_nonlinear(refTraj);

%% Save Results
refTraj.objective = opt.objective;

result_filename = "refTraj_" + opt.objective;
if opt.constraints.skew
    result_filename = result_filename + "_skew";
end

save(result_filename + ".mat", "refTraj","opt_hist")
fprintf("Saved results to %s.mat\n", result_filename)

%% Plotting
figure;
plot(refTraj.mu_hist(1,:),refTraj.mu_hist(2,:),'b-')
hold on
for k = 1:N_nodes
    xik = reshape(refTraj.xk_ref(:,k),nx,refTraj.ns);
    plot(xik(1,:),xik(2,:),'kx','LineStyle','none')
    hold on
end
grid on
axis equal

theta_plot = linspace(0,2*pi,1000);
for k = 1:N_nodes
    plot(refTraj.muk_ref(1,k)+constraint.r_3sigma*cos(theta_plot),refTraj.muk_ref(2,k)+constraint.r_3sigma*sin(theta_plot),'r:')
end
plot(constraint.muf(1)+constraint.r_3sigma*cos(theta_plot),constraint.muf(2)+constraint.r_3sigma*sin(theta_plot),'r:')
