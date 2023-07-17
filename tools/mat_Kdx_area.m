function Kdx = mat_Kdx_area(young,elem,coord,dof,adof)
    B = mat_B_area(young, elem, coord, dof, adof);
    nel = size(elem,1);
    Kdx = cell(nel,1);
    for i = 1:size(elem,1)
        Kdx{i} = B(:,i)*B(:,i)';
    end
end