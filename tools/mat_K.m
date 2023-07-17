function Matrix = mat_K(young_module,area,elem,coord,dof,adof)
    ndof = max(max(dof));
    Matrix = zeros(ndof,ndof);
    cathetus = coord(elem(:,2),:)-coord(elem(:,1),:); % nel x 3
    elem_length = vecnorm(cathetus'); % 1 x nel
    p = cathetus'./elem_length; % (3 x nel)./(1 x nel) = 3 x nel
    coef = young_module.*area./elem_length';
    for i = 1:size(elem,1)
        Q = [p(:,i) zeros(3,1);zeros(3,1) p(:,i)];
        MatElem = coef(i)*Q*[1 -1;-1 1]*Q';
        gl12 = dof(i,:);
        Matrix(gl12,gl12) = Matrix(gl12,gl12) + MatElem;
    end
    Matrix = Matrix(adof,adof); % reduzed global matrix
end