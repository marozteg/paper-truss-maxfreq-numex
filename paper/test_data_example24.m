%example 2.4 2008-SJO-Achtziger-Kocvara.pdf
clc
clear
close all

[coord, elem, dof, adof, E, rho, f] = data_example24();
m = size(elem,1);

%draw_elem(coord,elem,ones(m,1));
cathetus = coord(elem(:,2),:)-coord(elem(:,1),:); % nel x 3
elem_length = vecnorm(cathetus'); % 1 x nel

e = eye(m);
n = size(f,1);
Kdx = cell(m,1);
Mdx = cell(m,1);
for i = 1:m
    area = 1/elem_length(i)*e(:,i); % area for unitary volume of elem i
    Kdx{i} = mat_K(E*ones(m,1),area,elem,coord,dof,adof);
    mass = rho*1*e(:,i); % mass for unitary volume of elem i 
    Mdx{i} = mat_M(mass,elem,coord,dof,adof);
end
K0 = zeros(n,n);
M0 = zeros(n,n);
Kfun = @(x) sum_mat(1,m,K0,Kdx,x);
Mfun = @(x) sum_mat(1,m,M0,Mdx,x);

K1 = Kfun(ones(m,1));
M1 = Mfun(ones(m,1));

% in example 2.4 2008-SJO-Achtziger-Kocvara.pdf:
K2 = diag([2 2 1.28 0.32]);
M2 = diag([2.83 2.83 4.47 4.47]);
% 4.47 4.47 is impossible since cos and sin of element are differente in
% absolute value
[K1 K2]
[M1 M2]


