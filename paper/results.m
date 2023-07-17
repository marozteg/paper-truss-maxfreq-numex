clc
clear
close all

fprintf('shortcant\n')
load("data_shortcant","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

fprintf('shortcant2\n')
load("data_shortcant2","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

fprintf('mbbbeam\n')
load("data_mbbbeam","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

fprintf('halfwheel\n')
load("data_halfwheel","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

fprintf('hempcant\n')
load("data_hempcant","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

% fprintf('benchmark01\n')
% load("data_benchmark01","m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff")
% table(m,n,time_CPAS,time_PEN,lambda_CPAS,lambda_PEN,vol_CPAS,vol_PEN,compl_CPAS,compl_PEN,x_diff)

% % drawing all trusses
% problem = "data_benchmark01";
% load(problem,"x_CPAS_cell","x_PEN_cell")
% prob_sizes = length(x_PEN_cell);
% for i = 1:prob_sizes  
%     data_problem = str2func(problem);
%     [coord, elem, dof, adof, E, rho, f] = data_problem(i);
%     draw_truss(coord,elem,x_PEN_cell{i},[1,1],1e-6);
% end

