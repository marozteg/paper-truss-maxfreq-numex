function [coord, elem, dof, adof, E, rho, f] = data_michell()
% 7.4. Discretized Plate (Michell Truss)
% The ground structures are discretized plates in which each node is 
% connected to each other node. 
% The structures are fixed on the left side and a single vertical force is 
% applied at the mid righthand node. 

    coord = gen_coord(3,9,1);
    nn = size(coord,1);

    elem = gen_elem(coord,1);

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
        
    fdof = [];
    for i = 1:nn
        offs = (i-1)*3;
        if coord(i,1) == 0
            fdof = [fdof offs+1 offs+2];
        end
        if coord(i,3) == 0
            fdof = [fdof offs+3]; % all z are fixed
        end
    end
    fdof = unique(fdof);
    adof = 1:ndof;
    adof(fdof) = [];
    
    E = 1.0; % Pa
    
    rho = 1.0; % kg/m^3
    
    f = zeros(ndof,1);
    f((15-1)*3+2,1) = -0.1; % N
    f = f(adof,:);
end