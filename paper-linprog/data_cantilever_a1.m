function [coord, elem, dof, adof, E, rho, f] = data_cantilever_a1()

Lx = 4;%Lx = 40;%6; % must be 2 4 6 8 ...
Ly = 1;
nx = Lx + 1; % must be 3 5 7 9...
deltax = 0.25*Ly;
nn = nx*2;
coord = zeros(nn,3);
k = 1;
for i = 1:nx
    coordx = (i-1)*deltax;
    coord(k,1) = coordx;
    coord(k+1,1) = coordx;
    coord(k+1,2) = Ly;
    k = k + 2;
end

nel = 4*(nx-1) + nx - 1;
elem = zeros(nel,2);
r1 = 1:5;
for i = 1:nx-1
    elem(r1,:) = [1 3;2 4;1 4;2 3; 3 4] + 2*(i-1);
    r1 = r1 + 5;
end


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
    if coord(i,1) == 0 % left nodes are fixed
        fdof = [fdof offs+1 offs+2];
    end
    fdof = [fdof offs+3]; % all z are fixed
end
fdof = unique(fdof);
adof = 1:ndof;
adof(fdof) = [];

E = 1.0; % Pa

rho = 1.0; % kg/m^3

f = zeros(ndof,1);
f(((nn-1)-1)*3+2) = -1.0; % N
f = f(adof);