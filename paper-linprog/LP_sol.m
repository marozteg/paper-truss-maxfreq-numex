clear
close all

% lembrar de ajustar o volume v !!!!!!!!!!!!!!!
%[coord, elem, dof, adof, E, rho, f] = data_cantilever_a1(); % Example 7.1 a1 1992-ICSE-Achziger
%[coord, elem, dof, adof, E, rho, f] = data_cantilever_b1(); % Example 7.1 b1 1992-ICSE-Achziger
%[coord, elem, dof, adof, E, rho, f] = data_michell_a();     % Example 7.4 a 1992-ICSE-Achziger 
%[coord, elem, dof, adof, E, rho, f] = data_michell_c();     % Example 7.4 c 1992-ICSE-Achziger
[coord, elem, dof, adof, E, rho, f] = data_shortcant(2); % EVEN NUMBER!! % Example from 2022-EO-Weldeyesus-Gilbert-Gondzio
m = size(elem,1);
n = size(f,1);
%draw_truss(coord,elem,ones(m,1),[2 0]);

B = mat_B_vol(E*ones(m,1),elem,coord,dof,adof);
%B = mat_B_area(E*ones(nel,1),elem,coord,dof,adof);

A = [B';-B'];
b = ones(2*m,1);
opt = optimoptions('linprog','Display','iter');
[d,~,~,~,multiplier] = linprog(-f,A,b,[],[],[],[],opt);

v = 1;%7030;%5;
sigma = multiplier.ineqlin;
sumsigma = sigma(1:m) + sigma(m+1:end);
mu = (1/v)*ones(m,1)'*sumsigma;
u = mu*d;
xlp = 1/mu*sumsigma;

draw_truss(coord,elem,xlp,[2 1],1e-6);

Mdx = mat_Mdx_vol(rho*ones(m,1),elem,coord,dof,adof);
Kdx = mat_Kdx_vol(E*ones(m,1),elem,coord,dof,adof);
M0 = zeros(n,n);
Mfun = @(x) sum_mat(1,m,M0,Mdx,x);
K0 = zeros(n,n);
Kfun = @(x) sum_mat(1,m,K0,Kdx,x);

fprintf('Compliance at solution: %1.3e\n', f'*u)
fprintf('smallest freq at solut: %1.3e\n', maxpsdpair(Kfun(xlp),Mfun(xlp)))

% checking the idea of initial relaxation !!!!
[V,D] = eig([f'*u -f';-f Kfun(xlp)]);
% eigenval = diag(D);
% j_pos_eigenval = find(eigenval>1e-7);
% V = V(:,j_pos_eigenval);
r = size(V,2); % the same as rank(D)=rank([gammaLP -f';-f KK])
A = zeros(r,m+1); b = zeros(r,1);
for i = 1:r
    h = V(1,i);
    w = V(2:end,i);
    for j = 1:m
        A(i,j) = - w'*Kdx{i}*w;
    end
    A(i,m+1) = - h*h;
    b(i) = w'*K0*w - 2*h*f'*w;
end
c = [zeros(m,1);1];
A = [A;[ones(1,m) 0]];
b = [b;v];
lb = [zeros(m,1);0];
ub = [];
opt = optimoptions('linprog','Display','iter');
[xsol2,~,~,~,~] = linprog(c,A,b,[],[],lb,ub,opt);
[elem xlp xsol2(1:m)]

