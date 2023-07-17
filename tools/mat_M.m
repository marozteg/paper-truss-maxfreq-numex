function Matrix = mat_M(mass,elem,coord,dof,adof)
    ndof = max(max(dof));
    Matrix = zeros(ndof,ndof);
    r = coord(elem(:,2),:) - coord(elem(:,1),:);
    elem_length = vecnorm(r')';
    coef = mass.*elem_length/6;
    for i = 1:size(elem,1)
        MatElem = coef(i)*( diag(2*ones(6,1)) + diag(ones(3,1),3) + diag(ones(3,1),-3) );
        gl12 = dof(i,:);
        Matrix(gl12,gl12) = Matrix(gl12,gl12) + MatElem;
    end
    Matrix = Matrix(adof,adof); % reduzed global matrix 
end
