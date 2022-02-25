clear all;
close all;
clc;

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();

%set the parameters for running registration process
disp('Setting parameters...');

%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-04;
%maximum number of iterations for each refinement level
itermax = 4000;
updateFlag = 0;
alpha = [1,2,2,3,2];
load theta_embryo.mat;

theta_g_embryo = theta_g;
theta_g_embryo1 = [theta_g_embryo(1);theta_g_embryo(5);theta_g_embryo(9);theta_g_embryo(13);theta_g_embryo(17)];
disp('Reading image data...');
numberIndex = [240;299;359;419;479];
myVars = {'Embryo240','Embryo299','Embryo359','Embryo419','Embryo479'};

%% 1: Read image data
%Enter file name here
iterct = 0;
iteration = [];
rresidual_img = [];
rresidual_ssd = [];
rresidual_error = [];
rresidual_reg = [];
rresidual_regl = [];
rresidual_total = [];

%Store the pixel coordinates
disp('Storing pixel coordinates...');
Nx = 100;
Ny = 100;
Nz = 100;

[pixX, pixY, pixZ]= ndgrid(linspace(1,Nx,Nx),linspace(1,Ny,Ny),linspace(1,Nz,Nz));
X = pixX(:);
Y = pixY(:);
Z = pixZ(:);
pix = [X,Y,Z];

[Xvtk,Yvtk,Zvtk] =  meshgrid(1:Nx,1:Ny,1:Nz);

pxx = Yvtk;
pyy = Xvtk;
pzz = Zvtk;

beta_1 = 10;
sizeImage = [Nx,Ny,Nz];
load 'thetag_result_4.mat'
vars = {'theta_plot_avg', 'theta_plot_max','theta_plot_center'};
clear(vars{:})

for nsteps = 1:1
    
    fprintf('Nsteps = %i\n',nsteps);
    
    disp('Reading image data...');
    str1 = 'Embryo_img_240.mat';
    str2 = sprintf('Embryo_img_%d.mat',numberIndex(nsteps+1));
    
    img1 = load(str1);
    img2 = load(str2);
    
    Img11 = img1.Img_source;
    Img22 = img2.Img_target;
    
    Img_source = Img11;
    Img_target = Img22;
    
    bw = 30;
    
    Img11 = zeros(size(Img_source,1)+bw,size(Img_source,2)+bw,size(Img_source,3)+bw);
    Img22 = zeros(size(Img_source,1)+bw,size(Img_source,2)+bw,size(Img_source,3)+bw);
    
    st = floor(bw/2);
    Img11(st+1:size(Img11,1)-st,st+1:size(Img11,2)-st,st+1:size(Img11,3)-st) = Img_source;
    Img22(st+1:size(Img22,1)-st,st+1:size(Img22,2)-st,st+1:size(Img22,3)-st) = Img_target;
    
    Img1 = resize3Dmatrix(100,100,100,Img11);
    Img2 = resize3Dmatrix(100,100,100,Img22);
    
    Img_source = Img1;
    Img_target = Img2;
    
    Img_source_initial = Img_source; %store initial moving image
    Img_target_initial = Img_target; %store initial target image
    
    filename1 = sprintf('post_processing/source_%d.vtk',nsteps);
    filename2 = sprintf('post_processing/target_%d.vtk',nsteps);
    vtkwrite(filename1, 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_source);
    vtkwrite(filename2, 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_target);
    
    %compute initial RS value
    disp('Compute initial RS...');
    Idiff = Img_source-Img_target;
    Idiff2 = (Idiff).^2;
    residual_initial = sqrt(sum(sum(sum(Idiff2,3),2),1));
    RS_initial = 1;
    disp(sprintf('initial residual %f\n',residual_initial));
    disp(sprintf('initial RS value %f\n',RS_initial));
    
    if(nsteps<3)
        param = setparameters_embryo_growth(2,alpha(nsteps));
    else
        param = setparameters_embryo_growth(3,alpha(nsteps));
    end
    
    kappa = theta_g_embryo1(nsteps+1)-theta_g_embryo1(nsteps);
    dtime = 1;
    
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
    ActiveNodes = [];
    Node = [];
    
    for levels = 1:param.maxlevel
        [nx1,ny1,nz1] = meshgrid(uknotvectorV{levels,1},uknotvectorU{levels,1},uknotvectorW{levels,1});
        Node = [Node;ny1(:),nx1(:),nz1(:)];
    end
    
    % Loop over each refinement level
    disp('Loop over each refinement level...');
    for multilev = 0:1:param.maxlevel-1
        
        [DIITX,DIITY,DIITZ] = gradient(Img_target);
        
        fprintf('Refinement at level %i...\n',multilev+1);
        if(multilev>0)
            [Em,Dm,Pm,ActiveNodes] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multilev,ActiveNodes);
        end
        
        disp('Collecting active elements, control points and basis functions...');
        [ac, bf, ACP, RHS,Em,Dm,ActiveNodes] = storeActiveElem(Em,Dm,Pm,multilev,ActiveNodes);
        
        ac_ct = size(ac,1);
        bf_ct = size(bf,1);
        disp(sprintf('active elements = %d, active DOF =  %d at level = %d \n',ac_ct,bf_ct,multilev+1));
        
        ActiveNodes = unique(ActiveNodes);
        
        Node(:,4) = 0;
        for i=1:size(ActiveNodes,1)
            Node(ActiveNodes(i,1),4) = i;
        end

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
        
        cII_SX = interp3(pixY, pixX, pixZ, Img_source_initial, BIGY, BIGX, BIGZ,'*linear',min(Img_target(:)));
        
        if(nsteps==1)
            theta_g_new = ones(ac_ct,param.orderGauss,param.orderGauss,param.orderGauss);
            theta_g_old = theta_g_new;
            theta_plot_max(1,1) = 1;
            theta_g = ones(Nx,Ny,Nz);
            theta_g_new_cp = ones(bf_ct,1);
        else
            theta_g_old1 = interp3(pixY, pixX, pixZ, theta_g, BIGY, BIGX, BIGZ,'*linear',1);
            theta_g_old = reshape(theta_g_old1,[ac_ct,param.orderGauss,param.orderGauss,param.orderGauss]);
            theta_g_new_cp = interp3(pixY, pixX, pixZ, theta_g, ACP(:,1), ACP(:,2), ACP(:,3),'*linear',1);
        end
        
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
        filename = sprintf('post_processing/Embryo_mesh_growth_undeformed_%d.vtk',nsteps);
        PlotGrid(filename,ActiveNodes,Node,pxx,pyy,pzz,Em,ac_ct,ac);
        
        if((multilev+1) == param.maxlevel)
            itermax = 500;
        end
        
        % while the stopping criterion is satisfied
        while (iterct_level < itermax)
            % counting total iterations
            iterct = iterct +1;
            
            % iterations for each refinement level
            iterct_level = iterct_level+1;
            RHS = RHS_init;
            
            % compute the spatial transformation function f(x), f_u(x), f_v(x), f_w(x)
            [BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ] = computenewPoints_mex(Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,orderGauss);
            
            % interpolate the intensity and grdient at f(x) at the deformed positions of the gauss points
            cII_TY = interp3(pixY, pixX, pixZ, Img_target, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));

            cDII_TY_X = interp3(pixY, pixX, pixZ, DIITX,BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            cDII_TY_Y = interp3(pixY, pixX, pixZ, DIITY, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            cDII_TY_Z = interp3(pixY, pixX, pixZ, DIITZ, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));

            % fidelity term in x, y, z directions
            Bterm1 = (cII_TY - cII_SX).*2.*cDII_TY_Y; %./denominate;
            Bterm2 = (cII_TY - cII_SX).*2.*cDII_TY_X; %./denominate;
            Bterm3 = (cII_TY - cII_SX).*2.*cDII_TY_Z; %./denominate;
            
            %Bterm = [Bterm1;Bterm2;Bterm3];
            img_error = (cII_TY-cII_SX).^2;
            
            % Now compute the integrals for each basis function in the support of the active cell.
            if(mod(iterct_level,itermax)==0 && multilev == (param.maxlevel-1))
                [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error, theta_g_new] = compute_Integ_Domain_hyper_elastic_growth_projection_1_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,2.5, Wu,Wv,Ww,H,cII_SX,kappa,dtime,theta_g_old);
                clear Mass_matrix;
            else
                [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error] = compute_Integ_Domain_hyper_elastic_theta_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,beta_1,Wu,Wv,Ww,H,cII_SX,theta_g_old);
            end
            
            if(mod(iterct_level,itermax)==0 && multilev == (param.maxlevel-1))
                [pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy,~] = tripleIterLoop_growth2_mex(sizeImage, Pixel, Jm, ACP, Img_source,param.mu,param.lambda, max(theta_g_new(:)));
            end
            
            % apply dirichlet boundary condition on the control points at the
            % boundary
            RHS = bcondition3D(RHS);
            RHS_reg = bcondition3D(RHS_reg);
            RHS_image = bcondition3D(RHS_image);
            
            vec_reg = RHS_reg(:,1:3);
            vec_reg = vec_reg(:);
            norm_reg = norm(vec_reg);
            
            vec_img = RHS_image(:,1:3);
            vec_img = vec_img(:);
            norm_img = norm(vec_img);
            
            disp(sprintf('Iteration = %d, residual = %f',iterct_level, final_error));
            rs(iterct,1) = final_error;
            iter_reg(iterct,1) = norm_reg;
            iter_img(iterct,1) = norm_img;
            iterations(iterct,1) = iterct;
            theta_plot_avg(iterct,1) = mean(theta_g_new(:));
            theta_plot_max(iterct,1) = max(theta_g_new(:));
            temp_theta = theta_g_new(786,:,:,:);
            theta_plot_center(iterct,1) = mean(temp_theta(:));
            
            subplot(2,2,1)
            plot(iterations, rs);
            subplot(2,2,2)
            plot(iterations, iter_reg);
            subplot(2,2,3)
            plot(iterations,theta_plot_avg);
            subplot(2,2,4)
            plot(iterations,theta_plot_max);
            drawnow
            
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
            
            save('post_processing/embryo_growth_residuals.mat','iter_reg','iter_img','iterations');
            if(mod(iterct_level,2000)==0)
                param.alpha = param.alpha*2;
            end
            
            if(mod(iterct_level,500)==0)
                %compute the new positions of the pixel coordiantes using the spatial
                [pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy,theta_g] = tripleIterLoop_growth2_mex(sizeImage, Pixel, Jm, ACP, Img_source,param.mu,param.lambda,max(theta_g_new(:)));
                
                Img_target_reparameterized_at_pixels = interp3(pixY, pixX, pixZ, Img_target, pyy, pxx, pzz,'*linear',min(Img_target(:)));
                
                filename_vtk = sprintf('post_processing/evolve_embryo_growth_%d.vtk',iterct);
                vtkwrite(filename_vtk, 'structured_grid',pyy,pxx,pzz,'scalars','StrainEnergy',strain_energy,'scalars','J',Jdet,'scalars','F11',F11,'scalars','F22',F22,'scalars','F33',F33,'scalars','thetag',theta_g);
            end
        end
        
        save('post_processing/thetag_result_5.mat','theta_g_new','theta_plot_avg','theta_plot_max','theta_plot_center','theta_g');
        
        %the centroids of the B-spline grid
        cell_co = zeros(ac_ct,3);
        for i = 1:ac_ct
            cell_id = ac(i,1);
            cell_le = ac(i,2);
            cell_co(i,:) = Em(cell_le).cell_centre(cell_id,:);
        end
        
        %compute the absolute difference between the pair of images at centroid
        %of elements
        [CellGrad, meanGrad] = computeDiffGradImage(cell_co,Img_source,Img_target, pixX, pixY, pixZ);
        fprintf('Mean Gradient = %f\n',meanGrad);
    end
    
    itermax = itermax+250;

    filename_vtk1 = sprintf('post_processing/embryo_growth_deformed_%d.vtk',nsteps);
    vtkwrite(filename_vtk1, 'structured_grid',pyy,pxx,pzz,'scalars','Intensity',Img_source);
    
    filename1 = sprintf('post_processing/Embryo_mesh_growth_%d.vtk',nsteps);
    PlotGrid(filename1,ActiveNodes,Node,pxx,pyy,pzz,Em,ac_ct,ac);
end
