%%
% Performs Monte Carlo Analysis of Solution 
%%

clear
close all

initialization;
params_simulation;

%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define MC
%MC_indices = 0:100;
MC_indices = 0:10000;
%MC_indices = 0:1000;

%% Select SMS Solution
solution_name = "CS";
%solution_name = "quadratic";
%solution_name = "quadratic_skew";
%solution_name = "cvar";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solution_filename = "refTraj_" + solution_name + ".mat";
load(solution_filename, "refTraj");

% Save MC Results
results_X_sol = cell(length(MC_indices),1);
results_Xk = nan(nx,N_nodes,length(MC_indices));
results_uk_hist = nan(nu,N_edges,length(MC_indices));

%% Run Sim
parfor i = 1:length(MC_indices)
%for i = 1:length(MC_indices)
    mc_index = MC_indices(i);
    fprintf("MC Run: %d\n", mc_index);

    [X_sol, Xk_hist, uk_hist] = main_SIM(refTraj);
    
    % [X_sol, Xk_hist, u_k_hist] = main_SIM(refTraj)
    % results_X_sol{i} = X_sol; % A lot of data
    
    results_Xk(:,:,i) = Xk_hist;
    results_uk_hist(:,:,i)= uk_hist;
end

%% Results Processing
muk_MC = zeros(nx,N_nodes);
Pk_MC = zeros(nx,nx,N_nodes);
gammak_MC = zeros(nx,N_nodes);
kurtk_MC = zeros(nx,N_nodes);
for k = 1:N_nodes
    Xk_MC = squeeze(results_Xk(:,k,:));

    muk_MC(:,k) = mean(Xk_MC,2);
    Pk_MC(:,:,k) = cov(Xk_MC.');
    gammak_MC(:,k) = skewness(Xk_MC,0,2);
    kurtk_MC(:,k) = kurtosis(Xk_MC,0,2);
end

deltaV_mag_hist = squeeze(vecnorm(results_uk_hist,2,1));
deltaV_total_MC = squeeze(sum(vecnorm(results_uk_hist,2,1),2));
ub_k_mag = vecnorm(refTraj.ub_k,2,1);
deltaV_99_hist_MC = zeros(N_edges,1);
for k = 1:N_edges
    deltaV_99_hist_MC(k) = prctile(deltaV_mag_hist(k,:),99);
end

%% Save MC Control Results
results_filename = "results_MC_" + solution_name + ".mat";
save(results_filename, ...
    "solution_name", ...
    "solution_filename", ...
    "MC_indices", ...
    "results_uk_hist", ...
    "deltaV_mag_hist", ...
    "deltaV_total_MC", ...
    "deltaV_99_hist_MC");
fprintf("Saved MC control results to %s\n", results_filename);

%% Plotting
theta_plot = linspace(0,2*pi,1000);

%% Trajectory
fig = figure;
set(fig, 'Position',  [100,200,1300,400])
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
for k = 1:N_nodes
    Xk_MC = squeeze(results_Xk(:,k,:));
    plot(Xk_MC(1,:),Xk_MC(2,:),'k.')
    hold on
end
plot(refTraj.mu_hist(1,:),refTraj.mu_hist(2,:),'b-')
for k = 1:N_nodes
    xik_ref = reshape(refTraj.xk_ref(:,k),nx,refTraj.ns);
    plot(xik_ref(1,:),xik_ref(2,:),'c.')
end
plot(constraint.muf(1),constraint.muf(2),'rx','MarkerSize',12)
plot(constraint.r_3sigma*cos(theta_plot) + constraint.muf(1), constraint.r_3sigma*sin(theta_plot) + constraint.muf(2), 'r:')

grid minor
axis equal

% Direction of motion
quiver(2.4,-0.4,0,0.8,0, ...
    'Color','k','LineWidth',2,'MaxHeadSize',0.9, ...
    'HandleVisibility','off');

xlabel('$x$')
ylabel('$y$')

ax = gca;
ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')

h_mu = plot(refTraj.mu_hist(1,:),refTraj.mu_hist(2,:),'b-');
h_sigma = plot(nan,nan,'c.','MarkerSize',24); 
h_mc = plot(nan,nan,'k.','MarkerSize',24); 
h_target = plot(nan, nan, 'rx', 'MarkerSize',12);
h_const  = plot(nan, nan, 'r:', 'LineWidth',3);
legend([h_mu, h_sigma, h_mc, h_target, h_const], ...
    {'Mean (Predicted)', ...
    'Cubature Points', ...
    'MC Sample' ...
    'Mean Constraint', ...
    '$3\sigma$ Constraint'}, ...
    'Orientation','horizontal', ...
    'Location','northoutside');

%% Joint & Marginal Plot
% Select Results
%k_test = 11;
k_test = N_nodes;

axis_ind = [1,2];
axis_labels = {'$x$','$y$'};
axis_center = refTraj.muk_ref(axis_ind,k_test); % orgin at this point

target = constraint.muf(axis_ind);

Xk_mc = squeeze(results_Xk(axis_ind,k_test,:)); % x & y
data_XY = Xk_mc - axis_center;

fig = figure;
set(fig, 'Position',  [50,60,700,700])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax_mat = gobjects(2,2);   % store axes handles

% Joint Plot
ax = nexttile(3);
% MC Points
plot(data_XY(1,:),data_XY(2,:),'k.')
hold on
% Reference 3-sigma
[X_E_ref,Y_E_ref] = plotCovEllipse2D(refTraj.Pk_ref(axis_ind,axis_ind,k_test),3);
plot((X_E_ref +  refTraj.muk_ref(axis_ind(1),k_test)) - axis_center(1), (Y_E_ref + refTraj.muk_ref(axis_ind(2),k_test)) - axis_center(2), 'b-','LineWidth',3)
hold on
% MC 3-sigma
[X_E_MC,Y_E_MC] = plotCovEllipse2D(Pk_MC(axis_ind,axis_ind,k_test),3);
plot((X_E_MC + muk_MC(axis_ind(1),k_test)) - axis_center(1), (Y_E_MC + muk_MC(axis_ind(2),k_test)) - axis_center(2), '--','Color',colors.green,'LineWidth',3)
% Mean
plot(muk_MC(axis_ind(1),k_test) - axis_center(1),muk_MC(axis_ind(2),k_test) - axis_center(2),'Marker','s','Color',colors.green,'MarkerSize',12)
hold on
% Constraint
if k_test == N_nodes
    plot(target(1) - axis_center(1),target(2) - axis_center(2),'rx','MarkerSize',12)
    plot((constraint.r_3sigma*cos(theta_plot) + target(1)) - axis_center(1), (constraint.r_3sigma*sin(theta_plot) + target(2)) - axis_center(2), 'r','LineStyle',':','LineWidth',3)
end
grid minor
axis equal
xlabel(axis_labels{1})
ylabel(axis_labels{2})

ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')

Xlim_joint = ax.XLim;
Ylim_joint = ax.YLim;

ax_mat(2,1) = ax;

% Marginal (axes 1)
ax = nexttile(1);
pd = fitdist(data_XY(1,:).','Normal');

pdf_values = linspace(Xlim_joint(1),Xlim_joint(2),1000);

y = pdf(pd,pdf_values);
histogram(data_XY(1,:),'FaceColor','k','Normalization','pdf', 'Orientation','vertical');
hold on
plot(pdf_values,y,'r-')
set(gca,'YTickLabel',[]); % We don't need ylabel
skew_title = sprintf("%.3f",gammak_MC(axis_ind(1),k_test));
title("$\gamma="+ skew_title+"$")
grid minor
ylabel("pdf("+axis_labels{1}+")")

ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')

ax_mat(1,1) = ax;

% Marginal (axes 2)
ax = nexttile(4);
pd = fitdist(data_XY(2,:).','Normal');

pdf_values = linspace(Ylim_joint(1),Ylim_joint(2),1000);

y = pdf(pd,pdf_values);
histogram(data_XY(2,:),'FaceColor','k','Normalization','pdf', 'Orientation','horizontal');
hold on
plot(y,pdf_values,'r-')
set(gca,'XTickLabel',[]); % We don't need xlabel
skew_title = sprintf("%.3f",gammak_MC(axis_ind(2),k_test));
title("$\gamma="+ skew_title+"$")
grid minor
xlabel("pdf("+axis_labels{2}+")")

ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')

ax_mat(2,2) = ax;

% ---- Link axes ----
linkaxes([ax_mat(1,1), ax_mat(2,1)], 'x')
linkaxes([ax_mat(2,1), ax_mat(2,2)], 'y')

xlim(ax_mat(2,1), Xlim_joint)
ylim(ax_mat(2,1), Ylim_joint)
axis(ax_mat(2,1),'equal')

% Legend
ax_leg = nexttile(2);   % empty tile (top-right)

hold(ax_leg,'on')

% Create dummy handles matching your plot styles
h_data   = plot(ax_leg, nan, nan, 'k.','MarkerSize',24);
h_ref    = plot(ax_leg, nan, nan, 'b-', 'LineWidth',3);
h_mc     = plot(ax_leg, nan, nan, '--', 'Color',colors.green, 'LineWidth',3);
h_mean   = plot(ax_leg, nan, nan, 's', 'Color',colors.green, 'MarkerSize',12);
h_target = plot(ax_leg, nan, nan, 'rx', 'MarkerSize',12);
h_const  = plot(ax_leg, nan, nan, 'r:', 'LineWidth',3);
h_Gaussian  = plot(ax_leg, nan, nan, 'r-', 'LineWidth',3);

lgd = legend(ax_leg, ...
    [h_data, h_ref, h_mc, h_mean, h_target, h_const, h_Gaussian], ...
    {'MC Sample','$3\sigma$ (Predicted)','$3\sigma$ (MC)','Mean (MC)','Mean Constraint','$3\sigma$ Constraint','Fitted Gaussian'}, ...
    'Location','northwest');

lgd.FontSize = 20;          % bigger text
lgd.FontWeight = 'bold';    % bold text
lgd.LineWidth = 1.5;        % thicker legend box
lgd.Box = 'on';             % keep box (IEEE style)

axis(ax_leg,'off')

%% Control 
% Total
fig = figure;
set(fig, 'Position',  [100,100,800,500]);

tiledlayout(1,1);
ax = nexttile;
histogram(deltaV_total_MC,'FaceColor','k','Normalization','percentage');
hold on
%xline(refTraj.deltaV_avg,'b-','LineWidth',2);
hold on
xline(mean(deltaV_total_MC),'-','LineWidth',2,'Color',colors.green);
hold on
%xline(refTraj.deltaV_99),'b--','LineWidth',2);
hold on
xline(prctile(deltaV_total_MC,99),'--','LineWidth',2,'Color',colors.green);
grid on
xlabel('$\sum_k ||\mathbf{u}_k||_2$', 'Interpreter', 'latex')
xlim([0,inf])
ylabel('\% in Bin')
%lgd = legend('','Mean (CUT)','Mean (MC)','$\Delta V_{99,\mathrm{ub}}$ (CUT)','$\Delta V_{99}$ (MC)','Orientation','horizontal');
lgd = legend('','Mean (MC)','$99$-th Percentile (MC)','Orientation','horizontal');

lgd.Layout.Tile = 'north';
lgd.FontSize = 20;          % bigger text
lgd.FontWeight = 'bold';    % bold text
lgd.LineWidth = 1.5;        % thicker legend box
lgd.Box = 'on';             % keep box (IEEE style)

ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')


%% (Takes too long to plot if 10k)
if 0
%% Control History 
axis_labels = {'$u_1$','$u_2$','$||\cdot||_2$'};

fig = figure;
set(fig, 'Position',  [100,100,700,600]);
t = tiledlayout(3,1);

for axes = 1:2
    ax = nexttile(axes);

    for i = length(MC_indices):-1:1 
        % Make sure zoh lasts until Tf
        stairs([1:N_edges,N_edges*10],[results_uk_hist(axes,:,i),results_uk_hist(axes,end,i)],'Color',colors.gray)
        hold on
    end

    % Nominal
    stairs([1:N_edges,N_edges*10],[refTraj.ub_k(axes,:),refTraj.ub_k(axes,end)],'r--','LineWidth',1,'Marker','pentagram','MarkerFaceColor','r','MarkerSize',12)
    
    xlim([1,N_nodes])
    ylabel(axis_labels{axes});
    grid minor

    ax.FontSize = 24;        % tick labels
    ax.LineWidth = 1.5;      % axis thickness
    set(gcf,'PaperPositionMode','auto')
end
% Magnitude
ax = nexttile(3);
for i = length(MC_indices):-1:1
    stairs([1:N_edges,N_edges*10],[deltaV_mag_hist(:,i);deltaV_mag_hist(end,i)] ,'Color',colors.gray)
    hold on
end
% stairs([1:N_edges,N_edges*10],[refTraj.deltaV_99_hist;refTraj.deltaV_99_hist(end)],'b--','LineWidth',1, ...
%                                                                                                 'Marker','pentagram','MarkerFaceColor','b','MarkerSize',12)
stairs([1:N_edges,N_edges*10],[deltaV_99_hist_MC;deltaV_99_hist_MC(end)],'--','LineWidth',1, 'Color',colors.green,'Marker','pentagram','MarkerFaceColor',colors.green,'MarkerSize',12)
hold on

xlim([1,N_nodes])
ylabel(axis_labels{end});
grid minor

ax.FontSize = 24;        % tick labels
ax.LineWidth = 1.5;      % axis thickness
set(gcf,'PaperPositionMode','auto')

xlabel(t, '$k$','FontSize',24,'FontName', 'Times', 'Interpreter', 'latex')

h_data   = plot(nan, nan,'Color',colors.gray);
h_ref    = plot(nan, nan,'r--','LineWidth',1,'Marker','pentagram','MarkerFaceColor','r','MarkerSize',12);
h_mc     = plot(nan, nan,'--','LineWidth',1, 'Color',colors.green,'Marker','pentagram','MarkerFaceColor',colors.green,'MarkerSize',12);
lgd = legend([h_data, h_ref, h_mc], ...
    {'$\mathbf{u}_k$ (MC)','$\bar{\mathbf{u}}_k$','$99$-th Percentile (MC)'}, ...
    'Orientation','horizontal');

lgd.Layout.Tile = 'north';
lgd.FontSize = 20;          % bigger text
lgd.FontWeight = 'bold';    % bold text
lgd.LineWidth = 1.5;        % thicker legend box
lgd.Box = 'on';             % keep box (IEEE style)
end
