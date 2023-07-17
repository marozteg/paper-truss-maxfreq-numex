function [coord, elem, dof, adof, E, rho, f] = data_bridge_21()
coord = [
    0   0   0
    1   0   0
    2   0   0
    3   0   0
    4   0   0
    0   1   0
    1   1   0
    2   1   0
    3   1   0
    4   1   0
];
nn = size(coord,1);

elem = [
    1   2 % horizontais
    2   3
    3   4
    4   5
    6   7
    7   8
    8   9
    9   10
    1   6 % verticais
    2   7
    3   8
    4   9
    5   10
    1   7 % diagonais 45
    2   8
    3   9
    4   10
    2   6 % diagonais 90+45
    3   7
    4   8
    5   9
    % 1   8 % outras
    % 1   9
    % 1   10
    % 2   9
    % 2   10
    % 3   6
    % 3   10
    % 4   6
    % 4   7
    % 5   6
    % 5   7
    % 5   8
];
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

fdof = [1 2 3 13 14 15];
for i = 1:nn
    offs = (i-1)*3;
    fdof = [fdof offs+3]; % all z are fixed
end
fdof = unique(fdof);
adof = 1:ndof;
adof(fdof) = [];

E = 1.0; % Pa
rho = 1.0; % kg/m^3
f = zeros(ndof,1);
f(5) = -1.0;
f(8) = -1.0; % N
f(11) =-1.0;
f = f(adof);