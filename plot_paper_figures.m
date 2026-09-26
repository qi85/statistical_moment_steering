%%
% Generates misc plots from Monte Carlo results
%%

clear;
close all;

initialization;
params_simulation;

%% Plot Curve
[~,X_curve] = ode45(@(t,X) nonlinearEOM(t,X,params_ODE), [0,curve.Tf],curve.X0,options.ode);

theta_plot = linspace(0,2*pi,1000);

fig = figure;
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;

h_curve = plot(X_curve(:,1),X_curve(:,2),'b-');
hold on
h_start = plot(curve.X0(1),curve.X0(2),'Color',colors.green, 'Marker','pentagram','MarkerSize',20,'LineStyle','none');
h_target = plot(constraint.muf(1),constraint.muf(2),'rx','MarkerSize',12);
h_const = plot( ...
    constraint.r_3sigma*cos(theta_plot) + constraint.muf(1), ...
    constraint.r_3sigma*sin(theta_plot) + constraint.muf(2),'r:');

grid minor
axis equal

% Direction of motion
quiver(2.4,-0.4,0,0.8,0, ...
    'Color','k','LineWidth',2,'MaxHeadSize',0.9, ...
    'HandleVisibility','off');

xlabel('$x$')
ylabel('$y$')

ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1.5;
set(gcf,'PaperPositionMode','auto')

legend([h_curve,h_start,h_target,h_const], ...
    {'Reference','Start','Mean Constraint','$3\sigma$ Constraint'}, ...
    'Orientation','vertical', ...
    'NumColumns',1, ...
    'Location','eastoutside');

set(fig,'Position',[100,200,700,400])


%% Quadratic and CVaR Control-Cost Histograms
result_filenames = ["results_MC_quadratic.mat", "results_MC_cvar.mat"];
result_titles = ["Quadratic", "CVaR"];
control_costs = cell(size(result_filenames));

for i = 1:numel(result_filenames)
    results = load(result_filenames(i), "deltaV_total_MC");
    control_costs{i} = results.deltaV_total_MC(:);
end

all_control_costs = vertcat(control_costs{:});
%x_limits = [0, 1.02*max(all_control_costs)];
x_limits = [0, 1];

fig = figure;
set(fig, "Position", [100,100,900,600]);
t = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
ax = gobjects(numel(control_costs),1);
legend_handles = gobjects(2,1);

for i = 1:numel(control_costs)
    ax(i) = nexttile;
    histogram(control_costs{i}, ...
        "FaceColor","k", ...
        "Normalization","percentage");
    hold on
    h_mean = xline(mean(control_costs{i}), "-", ...
        "LineWidth",2, "Color",colors.green);
    h_99 = xline(prctile(control_costs{i},99), "--", ...
        "LineWidth",2, "Color",colors.green);
    if i == 1
        legend_handles = [h_mean;h_99];
    end

    title(result_titles(i));
    grid on

    ax(i).FontSize = 24;
    ax(i).LineWidth = 1.5;
end

linkaxes(ax,"x");
xlim(ax(1),x_limits);
xlabel(t,"$\sum_k ||\mathbf{u}_k||_2$", ...
    "FontSize",24, "FontName","Times", "Interpreter","latex");
ylabel(t,"\% in Bin", ...
    "FontSize",24, "FontName","Times", "Interpreter","latex");

lgd = legend(ax(1), legend_handles, ...
    {"Mean (MC)","$99$-th Percentile (MC)"}, ...
    "Orientation","horizontal");
lgd.Layout.Tile = "north";
lgd.FontSize = 20;
lgd.FontWeight = "bold";
lgd.LineWidth = 1.5;
lgd.Box = "on";

set(fig,"PaperPositionMode","auto");

%% SCvx* Convergence
refTraj_filenames = ["refTraj_quadratic.mat", ...
                     "refTraj_cvar.mat", ...
                     "refTraj_quadratic_skew.mat"];
solution_names = {"Quadratic", "CVaR", "Quadratic + Skewness"};
solution_colors = {"b", colors.green, "k"};
solution_markers = {"o", "s", "^"};

opt_hist_list = cell(size(refTraj_filenames));
max_iter = 1;
for i = 1:numel(refTraj_filenames)
    convergence_data = load(refTraj_filenames(i),"opt_hist");
    opt_hist_list{i} = convergence_data.opt_hist;
    max_iter = max(max_iter,max(opt_hist_list{i}.iter));
end
infeas_values = cellfun(@(hist) hist.infeas(:), ...
    opt_hist_list,"UniformOutput",false);
infeas_values = vertcat(infeas_values{:});
deltaJ_values = cellfun(@(hist) abs(hist.deltaJ(:)), ...
    opt_hist_list,"UniformOutput",false);
deltaJ_values = vertcat(deltaJ_values{:});

fig = figure;
t = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
convergence_ax = gobjects(2,1);
convergence_lines = gobjects(1,numel(opt_hist_list));

convergence_ax(1) = nexttile;
convergence_ax(1).YScale = "log";
hold on
yline(opt_hist_list{1}.opt_params.tol_opt,"r-", ...
    {"$\epsilon_{\mathrm{opt}} = 10^{-4}$"}, ...
    "LabelHorizontalAlignment","left", ...
    "LabelVerticalAlignment","bottom",...
    "FontSize",18, ...
    "LineWidth",1, ...
    "Interpreter","latex");
for i = 1:numel(opt_hist_list)
    opt_hist = opt_hist_list{i};
    marker_indices = unique([1:2:numel(opt_hist.iter),numel(opt_hist.iter)]);
    convergence_lines(i) = semilogy(opt_hist.iter,abs(opt_hist.deltaJ),"--", ...
        "Color",solution_colors{i}, ...
        "LineWidth",3, ...
        "Marker",solution_markers{i}, ...
        "MarkerIndices",marker_indices, ...
        "MarkerSize",12);
end
grid on
ylabel("$|\Delta J|$");
xlim([1,max_iter]);
setLogDecadeTicks(convergence_ax(1), ...
    [deltaJ_values;opt_hist_list{1}.opt_params.tol_opt]);

convergence_ax(2) = nexttile;
convergence_ax(2).YScale = "log";
hold on
yline(opt_hist_list{1}.opt_params.tol_feas,"r-", ...
    {"$\epsilon_{\mathrm{feas}} = 10^{-6}$"}, ...
    "LabelHorizontalAlignment","left", ...
    "LabelVerticalAlignment","bottom",...
    "FontSize",18, ...
    "LineWidth",1, ...
    "Interpreter","latex");
for i = 1:numel(opt_hist_list)
    opt_hist = opt_hist_list{i};
    marker_indices = unique([1:2:numel(opt_hist.iter),numel(opt_hist.iter)]);
    semilogy(opt_hist.iter,opt_hist.infeas,"--", ...
        "Color",solution_colors{i}, ...
        "LineWidth",3, ...
        "Marker",solution_markers{i}, ...
        "MarkerIndices",marker_indices, ...
        "MarkerSize",12);
end
grid on
ylabel("$\chi$");
xlabel("Iteration");
xlim([1,max_iter]);
setLogDecadeTicks(convergence_ax(2), ...
    [infeas_values;opt_hist_list{1}.opt_params.tol_feas]);

linkaxes(convergence_ax,"x");
for i = 1:numel(convergence_ax)
    convergence_ax(i).FontSize = 24;   % tick labels
    convergence_ax(i).LineWidth = 1.5; % axis thickness
end

leg = legend(convergence_ax(1),convergence_lines,solution_names, ...
    "Orientation","horizontal");
leg.Layout.Tile = "north";

fig.Position = [100,100,900,600];
set(fig,"PaperPositionMode","auto");

function setLogDecadeTicks(ax,values)
    values = values(isfinite(values) & values > 0);
    lower_exponent = floor(log10(min(values)));
    upper_exponent = ceil(log10(max(values)));
    if upper_exponent - lower_exponent < 2
        upper_exponent = lower_exponent + 2;
    end

    tick_exponents = round(linspace(lower_exponent,upper_exponent,3));

    ylim(ax,10.^[lower_exponent,upper_exponent]);
    yticks(ax,10.^tick_exponents);
    yticklabels(ax,compose("$10^{%d}$",tick_exponents));
    set(ax, ...
        "XMinorGrid","off", ...
        "YMinorGrid","off", ...
        "XMinorTick","off", ...
        "YMinorTick","off");
end
