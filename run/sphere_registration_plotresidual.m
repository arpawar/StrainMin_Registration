clc;
clear all;
close all;

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();

%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_sphere_local();

%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-04;
%maximum number of iterations for each refinement level
itermax = 500;

%% 1: Read image data
%Enter file name here
disp('Reading image data...');
sphere_img = zeros(50,50,50);
sphere_2x_img = zeros(50,50,50);

radius1 = 15;
radius2 = 18;

%Image size
sizeImage = size(sphere_img);
Nx = sizeImage(1,1);
Ny = sizeImage(1,2);
Nz = sizeImage(1,3);

for i =1:Nx
    for j=1:Ny
        for k=1:Nz
            if((i-Nx/2)^2+(j-Ny/2)^2+(k-Nz/2)^2<=radius1^2)
                sphere_img(i,j,k) = 1;
            end
            if((i-Nx/2)^2+(j-Ny/2)^2+(k-Nz/2)^2<=radius2^2)
                sphere_2x_img(i,j,k) = 1;
            end
        end
    end
end

img2 = sphere_2x_img;
img1 = sphere_img;

%Image size
sizeImage = size(img1);
Nx = sizeImage(1,1);
Ny = sizeImage(1,2);
Nz = sizeImage(1,3);

Img1 = img1;
Img2 = img2;

Img_source = Img1;
Img_target = Img2;

Img_source_initial = Img_source; %store initial moving image
Img_target_initial = Img_target; %store initial target image

%Image size
sizeImage = size(Img_source);
Nx = sizeImage(1,1);
Ny = sizeImage(1,2);
Nz = sizeImage(1,3);

PlotImage(1,Img_source,Img_target,Img_source_initial+Img_target_initial);

%compute initial RS value
disp('Compute initial RS...');
Idiff = Img_source-Img_target;
Idiff2 = (Idiff).^2;
residual_initial = sqrt(sum(sum(sum(Idiff2,3),2),1));
RS_initial = 1;
disp(sprintf('initial residual %f\n',residual_initial));
disp(sprintf('initial RS value %f\n',RS_initial));

%Store the pixel coordinates
disp('Storing pixel coordinates...');
[pixX, pixY, pixZ]= ndgrid(linspace(1,Nx,Nx),linspace(1,Ny,Ny),linspace(1,Nz,Nz));
X = pixX(:);
Y = pixY(:);
Z = pixZ(:);
pix = [X,Y,Z];

[Xvtk,Yvtk,Zvtk] =  meshgrid(1:Nx,1:Ny,1:Nz);

%Initialize the parametric space
disp('Store the knotvectors...');
[CP, knotvectorU, knotvectorV, knotvectorW, uknotvectorU, uknotvectorV, uknotvectorW] = storeKnotArray(param,sizeImage,param.maxlevel);

%% Compute Elem, Basis, Control point data structure for THBS splines
disp('Set B-spline grid...');

tic
[Dm,Em] = setBsplineGrid_func(knotvectorU, knotvectorV, knotvectorW, uknotvectorU,uknotvectorV,uknotvectorW,param,param.maxlevel);
toc


%% Store level 1 control points
disp('Compute control points at level 1 using Greville coordinates...');
knotu = knotvectorU{1,1};
knotv = knotvectorV{1,1};
knotw = knotvectorW{1,1};

pp = zeros(param.nobU(1,1)*param.nobV(1,1)*param.nobW(1,1),3);
for k = 1:param.nobW(1,1)
    for j =1:param.nobV(1,1)
        for i = 1:param.nobU(1,1)
            index =(k-1)*param.nobV(1,1)*param.nobU(1,1) + (j-1)*param.nobU(1,1) + i;
            %Greville Abscissae
            coordx = sum(knotu(i+1:i+param.pU))./param.pU;
            coordy = sum(knotv(j+1:j+param.pV))./param.pV;
            coordz = sum(knotw(k+1:k+param.pW))./param.pW;
            pp(index,:) = [coordx,coordy,coordz];
        end
    end
end

%Control points stored for each refienement level
CP(1).pts = pp;
Pm = CP;
Pm_undeformed = Pm;

ActiveNodes = [];
Node = [];

for levels = 1:param.maxlevel,
    [nx1,ny1,nz1] = meshgrid(uknotvectorV{levels,1},uknotvectorU{levels,1},uknotvectorW{levels,1});
    Node = [Node;ny1(:),nx1(:),nz1(:)];
end

figure
p = patch(isosurface(Xvtk,Yvtk,Zvtk,Img_source,0.5));
isonormals(Xvtk,Yvtk,Zvtk,Img_source,p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3);
axis tight
camlight
lighting gouraud
figure
p = patch(isosurface(Xvtk,Yvtk,Zvtk,Img_target,0.5));
isonormals(Xvtk,Yvtk,Zvtk,Img_target,p);
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3);
axis tight
camlight
lighting gouraud

vtkwrite('post_processing/source_sphere.vtk', 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_source);
vtkwrite('post_processing/target_sphere.vtk', 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_target);

iterct = 0;
iteration = [];
rresidual_img = [];
rresidual_ssd = [];
rresidual_error = [];
rresidual_reg = [];
rresidual_regl = [];
rresidual_total = [];

% Loop over each refinement level
disp('Loop over each refinement level...');
% AdrianQ: does registration begin from scratch for every refinement level?
% Aishwarya: This was what I was doing initially, restarting the grid, but
% it was wrong, removed the line for reinitializing control grid at the
% start of each refinement level
pxx = Xvtk;
pyy = Yvtk;
pzz = Zvtk;
figure
for multilev = 0:1:param.maxlevel-1
    
    % this is not correct
    % only gradient of target is needed
    [DIITX,DIITY,DIITZ] = gradient(Img_target);
    
    % confused about this line and the one below where Pm_old = Pm
    % is Pm used in the THB_refinement function? if not, if it is
    % just an output from the THB_refinement, do you really need this line?
    % Aishwarya: removed the reinitialization part
    %Pm = Pm_old; %Reinitialization of the control points done here
    
    %tic
    fprintf('Refinement at level %i...\n',multilev+1);
    if(multilev>0)
        [Em,Dm,Pm,ActiveNodes] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multilev,ActiveNodes);
    end
    %toc
    
    disp('Collecting active elements, control points and basis functions...');
    [ac, bf, ACP, RHS,Em,Dm,ActiveNodes] = storeActiveElem(Em,Dm,Pm,multilev,ActiveNodes);
    
    ac_ct = size(ac,1);
    bf_ct = size(bf,1);
    disp(sprintf('active elements = %d, active DOF =  %d at level = %d \n',ac_ct,bf_ct,multilev+1));
    
    ActiveNodes = unique(ActiveNodes);
    %store undeformed control points
    Pm_undeformed = Pm;
    %store the initial RHS value in RHS_init
    RHS_init = RHS;
    
    disp('Computing the non-zeros spline over each active element and storing coefficient matrices...');
    tic
    [Jm, Coeff] = computeNonZeroSplines(ac, param, Em, Dm);
    toc
    
    disp('Computing the basis functions at pixel coordinates...');
    tic
    numPixels = int64(prod(sizeImage));
    [Pixel, Pix2] = storePixelPhi_mex(numPixels, multilev,pix, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param);
    for p_ind = 1:numPixels
        Pixel(p_ind).phi = Pix2{p_ind,1};
        Pixel(p_ind).phiu = Pix2{p_ind,2};
        Pixel(p_ind).phiv = Pix2{p_ind,3};
        Pixel(p_ind).phiw = Pix2{p_ind,4};
    end
    clear Pix2
    toc;
    
    tic;
    disp('Computing the basis functions at gaussian points...');
    %compute the gaussian points and weights of the given gauss order
    [Gu,Wu] = ggquad(param.orderGauss);
    [Gv,Wv] = ggquad(param.orderGauss);
    [Gw,Ww] = ggquad(param.orderGauss);
    
    [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi_mex(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param);
    
    Node(:,4) = 0;
    for i=1:size(ActiveNodes,1)
        Node(ActiveNodes(i,1),4) = i;
    end

    % interpolate the intensity values of the source image at the gauss
    % points stored in BIGX, BIGY, BIGZ
    % AdrianQ: the statement is incorrect, we want to evaluate the source
    % at the undeformed gauss points. I'm changing the name to be a little
    % more clear
    
    cII_SX = interp3(pixY, pixX, pixZ, Img_source_initial, BIGY, BIGX, BIGZ,'*linear',min(Img_target(:)));
    
    %cII2 = interpolateBsplineFunc_gausspts(Nx,Ny,Nz,PHI,BIGY, BIGX, BIGZ,Jm,Pm,Dm,phi_cp);
    % AdrianQ: I dont know why you are clearing these, I think they might
    % be used later on?
    % Aishwarya: It is not used later, only deformed gauss points are used. Thus,
    % clearing it to not use too much RAM
    % clear('BIGX','BIGY','BIGZ');
    %toc;
    
    %% Start the image registration for each iteration
    disp('Starting the iteration loop for dynamic update of control points...');
    % start the iteration loop
    %% Update the iteration loop here
    
    % Gauss order
    orderGauss = param.orderGauss;
    
    % gamma term in g(x)
    smallNumber = param.smallNumber;
    
    % timestep for each refinement level
    timestep = param.timestep(multilev+1);
    
    % convert the cell array to struct array
    PHI1 = cell2struct(PHI,'mat',2);
    PHIU1 = cell2struct(PHIU,'mat',2);
    PHIV1 = cell2struct(PHIV,'mat',2);
    PHIW1 = cell2struct(PHIW,'mat',2);
    
    iterct_level = 0;
    RS_final = 2;
    level_i=150;
    
    PlotGrid_new
    % while the stopping criterion is satisfied
    while (iterct_level < itermax)
        % counting total iterations
        iterct = iterct +1;
        
        % iterations for each refinement level
        iterct_level = iterct_level+1;
        
        RS_initial = RS_final;
        RHS = RHS_init;
        
        % compute the spatial transformation function f(x), f_u(x), f_v(x), f_w(x)
        [BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ] = computenewPoints_mex(Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,orderGauss);
        
        % interpolate the intensity and grdient at f(x) at the deformed positions of the gauss points
        
        % AdrianQ: interpolating the target image function at the deformed
        % points T(Y)
        cII_TY = interp3(pixY, pixX, pixZ, Img_target, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
        
        % AdrianQ: interpolating the gradient of the target at Y
        cDII_TY_X = interp3(pixY, pixX, pixZ, DIITX,BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
        cDII_TY_Y = interp3(pixY, pixX, pixZ, DIITY, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
        cDII_TY_Z = interp3(pixY, pixX, pixZ, DIITZ, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
        
        % denominator of the fidelity term (g(x))
        % AdrianQ: this I don't know what it is, what is this used for?
        % Is it on the paper? This was never discussed
        %Aishwarya: speed term to speed up near higher gradient regions.
        %can remove it
        
        %denominate = sqrt((cDII_TY_X.^2) + (cDII_TY_Y.^2) + (cDII_TY_Z.^2)+ smallNumber); %g(x)
        
        % fidelity term in x, y, z directions
        Bterm1 = (cII_TY - cII_SX).*2.*cDII_TY_Y; %./denominate;
        Bterm2 = (cII_TY - cII_SX).*2.*cDII_TY_X; %./denominate;
        Bterm3 = (cII_TY - cII_SX).*2.*cDII_TY_Z; %./denominate;
        
        %Bterm = [Bterm1;Bterm2;Bterm3];
        img_error = (cII_TY-cII_SX).^2;
        
        % Now compute the integrals for each basis function in the support of the active cell.
        [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error] = compute_Integ_Domain_hyper_elastic_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,Wu,Wv,Ww,H,cII_SX);
        %RHS = compute_Integ_Domain(Jm,Bterm1,Bterm2,Bterm3,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,lambda_1,lambda_1,Wu,Wv,Ww,H);
        
        % apply dirichlet boundary condition on the control points at the
        % boundary
        RHS = bcondition3D(RHS);
        RHS_reg = bcondition3D(RHS_reg);
        RHS = bcondition3D(RHS);
        RHS_image = bcondition3D(RHS_image);
        LRHS_reg = bcondition3D(LRHS_reg);
        
        vec_reg = RHS_reg(:,1:3);
        vec_reg = vec_reg(:);
        norm_reg = norm(vec_reg);
        
        vec_total = RHS(:,1:3);
        vec_total = vec_total(:);
        norm_total = norm(vec_total);
        
        vec_lreg = LRHS_reg(:,1:3);
        vec_lreg = vec_lreg(:);
        norm_lreg = norm(vec_lreg);
        
        vec_img = RHS_image(:,1:3);
        vec_img = vec_img(:);
        norm_img = norm(vec_img);
        
        timestep
        % update the control points P_new = P_old - epsilon*delta(E)
        ACP(:,1:3) = ACP(:,1:3) - timestep.*RHS(:,1:3);
        for ilevel = 1:bf_ct
            level_b = bf(ilevel,2);
            Pmlevel = Pm(level_b).pts;
            index = bf(ilevel,1);
            Pmlevel(index,:) = ACP(ilevel,1:3);
            Pm(level_b).pts = Pmlevel;
        end
        %compute the new positions of the pixel coordiantes using the spatial
        %transformation function
        [pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy] = tripleIterLoop_mex(sizeImage, Pixel, Jm, ACP, Img_source,param.mu,param.lambda);
        
        %Using the new position of the pixels, compute the new moving image
        % AdrianQ: this is not correct, why are you interpolating both?
        % do we even need to interpolate? I don't think so, the above code
        % already computes the source image at the undeformed Gauss points
        % CII_SX; it also computes the image at the deformed Gauss points
        % CII_TY, why are you recomputing the image? I don't think this is
        % needed at all
        %Img_source = interp3(pixY, pixX, pixZ, Img_source_initial,pyy,pxx,pzz,'*linear');
        %Img_source(isnan(Img_source)) = min(Img_target(:));
        
        %Img_target = interp3(pixY, pixX, pixZ, Img_target_initial,pyy,pxx,pzz,'*linear');
        %Img_target(isnan(Img_target)) = min(Img_target(:));
        
        %compute the image gradient
        % AdrianQ: I also don't think this gradient is needed, you already
        % have the gradient function above, I dont think that changes, only
        % thing that changes is at which points you evaluate
        %[DI1X,DI1Y,DI1Z] = gradient(Img_target);
        
        %compute the residual and RS value
        % AdrianQ: the error can be integrated in the function that
        % calculates the residual, it can also be approximated from CII_SX
        % and CII_TY, which I know are not at the Gauss points but at the
        % pixel locations, but I am confident that the error is the
        % integral of (CII_TY-CII_SX).^2, if you want the error based on
        % the pixels that is also fine, but I don't think you should be
        % overwriting the original functions S() and T(), otherwise if you
        % change the functions S() and T() then the algorithm is doing
        % something different, it should be finding Y(X) such that
        % T(Y)=S(X), that's what it should be doing, given fixed functions
        % T() and S(), in the next few lines (and the previous ones) you
        % are actually changing the functions themselves, no?
        %Idiff2 = Img_source-Img_target;
        %Idiff22 = (Idiff2).^2;
        %Sum1 = sum(sum(sum(Idiff22,3),2),1);
        %final_residual = sqrt(Sum1);
        %RS_final = final_residual/residual_initial;
        % Same as the way the residual is computed, but interpolating at
        % the deformed pixels instead of deformed Gauss points
        
        Img_target_reparameterized_at_pixels = interp3(pixY, pixX, pixZ, Img_target, pyy, pxx, pzz,'*linear',min(Img_target(:)));
        %Idiff2 = Img_source-Img_target_reparameterized_at_pixels;
        %Idiff22 = (Idiff2).^2;
        %Sum1 = sum(sum(sum(Idiff22,3),2),1);
        %final_residual = sqrt(Sum1);
        %RS_final = final_residual/residual_initial;
        
        % AdrianQ
        % IF what you want is the deformed image you can just use the
        % original image values with the new pixel coordinates
        %Img_source_deformed = array[pxx;pyy;pzz;Img_source]
        % no interpolation needed, this is your new image, the same pixel
        % intensity, just different pixel position
        
        % Note that this is how we compute the residual, we take a pixel,
        % we see where it lands, we know that this pixel takes with it its
        % intensity, we want to see if the target image has the same
        % intensity at that point, that's why we do the T(Y), but the
        % deformed source is simple the same intensities per pixels, just
        % moved pixels
        
        disp(sprintf('Iteration = %d, residual = %f',iterct_level, final_error));
        rs(iterct,1) = final_error;
        iter_reg(iterct,1) = norm_reg;
        iter_total(iterct,1) = norm_total;
        iter_lreg(iterct,1) = norm_lreg;
        iter_img(iterct,1) = norm_img;
        iterations(iterct,1) = iterct;
        
        save('mu_1000.mat','iter_reg','iter_total','iter_lreg','iter_img','iterations');
        subplot(3,2,1)
        plot(iterations, rs);
        subplot(3,2,2)
        plot(iterations, iter_img);
        subplot(3,2,3)
        plot(iterations, iter_reg);
        subplot(3,2,4)
        plot(iterations, iter_lreg);
        subplot(3,2,5)
        plot(iterations, iter_total);
        drawnow
        
        if(mod(iterct_level,500)==0)
            filename_vtk = sprintf('post_processing/evolve_sphere_mu_1000_%d.vtk',iterct);
            vtkwrite(filename_vtk, 'structured_grid',pyy,pxx,pzz,'scalars','StrainEnergy',strain_energy,'scalars','J',Jdet,'scalars','F11',F11,'scalars','F22',F22,'scalars','F33',F33);          
            PlotGrid_new
        end
        vtkwrite('post_processing/sphere_deformed_mu_1000.vtk', 'structured_grid',pyy,pxx,pzz,'scalars','Intensity',Img_source);
    end
    
    
    
    %the centroids of the B-spline grid
    cell_co = zeros(ac_ct,3);
    for i = 1:ac_ct
        cell_id = ac(i,1);
        cell_le = ac(i,2);
        cell_co(i,:) = Em(cell_le).cell_centre(cell_id,:);
    end
    
    %compute the absolute difference between the pair of images at centroid
    %of elements
    % AdrianQ: I think this is not correct now that I changed your code
    % above for the Img_source and Img_target
    [CellGrad, meanGrad] = computeDiffGradImage(cell_co,Img_source,Img_target_reparameterized_at_pixels, pixX, pixY, pixZ);
    fprintf('Mean Gradient = %f\n',meanGrad);
    
    
end
