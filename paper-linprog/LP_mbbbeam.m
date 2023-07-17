clear
close all
clc
addpath('/home/marozteg/MATLAB-Drive/paper/examples')

% Example from 2022-EO-Weldeyesus-Gilbert-Gondzio
fprintf('given v and lambda\n')
fprintf('min gamma\n')
fprintf('s.t. gamma*K(x)-ff''>=0,\n')
fprintf('     sum(x)<=v and\n') 
fprintf('     x>=0.\n')
[coord, elem, dof, adof, E, rho, f] = data_mbbbeam(1);
m = size(elem,1);
n = size(f,1);

B = mat_B_vol(E*ones(m,1),elem,coord,dof,adof);
%B = mat_B_area(E*ones(nel,1),elem,coord,dof,adof);

A = [B';-B'];
b = ones(2*m,1);
opt = optimoptions('linprog','Display','iter');
[d,~,~,~,multiplier] = linprog(-f,A,b,[],[],[],[],opt);

v = 1;%7030/2; % m^3
sigma = multiplier.ineqlin;
sumsigma = sigma(1:m) + sigma(m+1:end);
mu = (1/v)*ones(m,1)'*sumsigma;
u = mu*d;
xLP = 1/mu*sumsigma;

draw_truss(coord,elem,xLP,[1 1],1e-6);

Mdx = mat_Mdx_vol(rho*ones(m,1),elem,coord,dof,adof);
Kdx = mat_Kdx_vol(E*ones(m,1),elem,coord,dof,adof);
M0 = sparse(n,n);
K0 = sparse(n,n);
Mfun = @(x) sum_mat(1,m,M0,Mdx,x);
Kfun = @(x) sum_mat(1,m,K0,Kdx,x);

fprintf('Compliance at solution: %1.3e\n', mincompl(Kfun(xLP),f))
fprintf('smallest freq at solut: %1.3e\n', maxpsdpair(Kfun(xLP),Mfun(xLP)))