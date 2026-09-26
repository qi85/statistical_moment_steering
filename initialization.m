% Add all subfolders to path 
addpath(genpath('Supporting_Scripts'));
addpath(pwd);

warning('off','MATLAB:deval:NonuniqueSolution'); % deltaV introduces non-continuous deval

%% Colors
colors.gray = [.7 .7 .7];
colors.green = "#00FF00";

%% Figure Formating
% Text Interpreters
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

% Fonts
set(0, 'defaultAxesFontName', 'Times');
set(0, 'defaultTextFontName', 'Times');

% Fontsizes
set(0, 'defaultTextFontSize', 16);
set(0, 'defaultAxesFontSize', 18);
set(0, 'defaultLegendFontSize', 14);
set(0, 'defaultColorbarFontSize', 14);
set(0, 'defaultTextarrowshapeFontSize', 14);

set(0, 'DefaultLineLineWidth', 2);