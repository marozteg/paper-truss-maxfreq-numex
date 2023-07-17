clc
clear
close all

nsize = 5;
areas = zeros(nsize,1);
volume = zeros(nsize,1);
gamma1 = zeros(nsize,1);
gamma2 = zeros(nsize,1);
for i = 1:nsize
    [coord, elem, dof, adof, E, rho, f] = data_benchmark01(i);

    n = length(f);
    m = size(elem,1);
    
    % draw_truss(coord,elem,ones(m,1),[1 2],0.0);
    % continue
    
    K0 = zeros(n,n);
    Kdx = mat_Kdx_vol(E*ones(m,1),elem,coord,dof,adof);
    M0 = zeros(n,n);
    Mdx = mat_Mdx_vol(rho*ones(m,1),elem,coord,dof,adof);
    Kfun = @(x) sum_mat(1,m,K0,Kdx,x);
    Mfun = @(x) sum_mat(1,m,M0,Mdx,x);
    
    cath = coord(elem(:,2),:)-coord(elem(:,1),:); % nel x 3
    L = vecnorm(cath'); % 1 x nel

    for area = 0.001:0.01:100
        xA = area*L;
        u = Kfun(xA)\f;
        if norm(u) <= 0.1*max(max(coord))
            break
        end
    end
    areas(i) = area;
    gamma1(i) = f'*u;
    volume(i) = sum(xA);

    B = mat_B_vol(E*ones(m,1),elem,coord,dof,adof);
    %B = mat_B_area(E*ones(nel,1),elem,coord,dof,adof);
    
    A = [B';-B'];
    b = ones(2*m,1);
    opt = optimoptions('linprog','Display','iter');
    [d,~,~,~,multiplier] = linprog(-f,A,b,[],[],[],[],opt);
    
    sigma = multiplier.ineqlin;
    sumsigma = sigma(1:m) + sigma(m+1:end);
    mu = (1/volume(i))*ones(m,1)'*sumsigma;
    u = mu*d;
    xlp = 1/mu*sumsigma;

    gamma2(i) = f'*u;
    
    %draw_truss(coord,elem,xlp,[1 1],1e-6);

end
[volume areas gamma1 gamma2]
