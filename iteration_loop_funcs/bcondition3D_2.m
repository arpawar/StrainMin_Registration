function BBM = bcondition3D_2(BBmatrix,BBV)

% In this function the dirichlet boundary condition is applied at the
% control points located at the boundary of the image domain
% The end condition is fixed, i.e. the displacement of the control points
% at the boundary is zero.

% INPUT:
% BBvector: the RHS matrix containing the the value of displacement of the
% control points, x,y,z directions. The fourth index corresponds to whether
% the control point lies on the boundary of the image domain (1) or not
% (0).

% OUTPUT:
% BBV: the output RHS matrix after applying the Dirichlet BC

%%
BBmatrix1 = BBmatrix;
for i =1:size(BBV,1)
    
    if(BBV(i,4)==1)
        
        BBmatrix1(i,:) = 0.0;
        BBmatrix1(:,i) = 0.0;
        BBmatrix1(i,i) = 1.0;
    end
    
    BBM = BBmatrix1;
    
end