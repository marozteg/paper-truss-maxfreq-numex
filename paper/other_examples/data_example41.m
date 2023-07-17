function [coord, elem, dof, adof, E, rho, f] = data_example41()
    % example 4.1 2008-SJO-Achtziger-Kocvara.pdf
    
    nx = 3;
    ny = nx;
    nz = 1;
    index = 1;
    coord = zeros(nx*ny*nz,3);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                coord(index,:) = [i-1 j-1 k-1];
                index = index + 1;
            end
        end
    end
    nn = size(coord,1);
    
    % all [a b] with a and b in 1:nn 
    [n1,n2] = meshgrid(1:nn,1:nn);
    elem = [n1(:) n2(:)];
    
    % Remove [a a]
    elem(elem(:,1)==elem(:,2),:) = [];
    
    % Remove [a b] if a>b
    elem = sort(elem,2);
    elem = unique(elem,'rows');
    
    % Remove bars longer than distmax
    n1 = elem(:,1);
    n2 = elem(:,2);
    c1x = coord(n1,1);
    c1y = coord(n1,2);
    c1z = coord(n1,3);
    c2x = coord(n2,1);
    c2y = coord(n2,2);
    c2z = coord(n2,3);
    dx = c2x-c1x;
    dy = c2y-c1y;
    dz = c2z-c1z;
    lengths = sqrt(dx.^2 + dy.^2 + dz.^2);
    distmax = max(lengths); % this way nothing is removed
    elem(lengths>distmax,:) = [];
    
    % Remove elements longer and in the same line than smaller ones
    nel = size(elem,1);
    elem2rem = false(nel,1);
    for i = 1:nel
        element = elem(i,:);
        n1 = element(1);
        n2 = element(2);
        v = [coord(n2,1)-coord(n1,1)
             coord(n2,2)-coord(n1,2)
             coord(n2,3)-coord(n1,3)];
        for n3 = 1:nn
            if n3 ~= n1 && n3 ~= n2
                w = [coord(n3,1)-coord(n1,1)
                     coord(n3,2)-coord(n1,2)
                     coord(n3,3)-coord(n1,3)];
                if norm(cross(v,w)) < 1e-6 % v and w are parallel
                    [~,j] = max(abs(v));
                    if 0 < w(j)/v(j) && w(j)/v(j) < 1
                        elem2rem(i) = true;
                    end
                end
            end
        end
    end
    elem(elem2rem,:) = [];

    % Remove elements with nodes x=0 (***)
    n1 = elem(:,1);
    n2 = elem(:,2);
    nel = size(elem,1);
    elem2rem = false(nel,1);
    for i = 1:nel
        if coord(n1(i),1) == 0 && coord(n2(i),1) == 0
            elem2rem(i) = true;
        end
    end
    elem(elem2rem,:) = [];

    % degree of freedom
    nel = size(elem,1);
    dof = zeros(nel,3*2);
    for i = 1:nel
        no1 = elem(i,1);
        no2 = elem(i,2);
        o1 = (no1-1)*3;
        o2 = (no2-1)*3;
        dof(i,:) = [o1+1 o1+2 o1+3 o2+1 o2+2 o2+3];
    end
    ndof = max(max(dof));
    if nn*3 ~= ndof
        keyboard
    end
        
    % nodes with x=0 are fixed (***)
    fdof = [];
    for i = 1:nn
        offs = (i-1)*3;
        if coord(i,1) == 0
            fdof = [fdof offs+1 offs+2 offs+3];
        else
            fdof = [fdof offs+3]; % all z are fixed
        end
    end
    fdof = unique(fdof);
    adof = 1:ndof;
    adof(fdof) = [];
    
    E = 1.0; % Pa
    
    rho = 1.0; % kg/m^3
    
    f = zeros(ndof,1);
    f((6-1)*3+1,1) = -0.461426;%-1.0; % N
    f = f(adof,:);
end