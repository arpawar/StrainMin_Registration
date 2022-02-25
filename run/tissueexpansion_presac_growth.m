clc;
clear all;
close all;

diary prefill_presac_growth.text

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();

%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-04;
%maximum number of iterations for each refinement level
itermax = 750;
updateFlag = 1;
alpha = [8,0.25,10];
kappa = 1;
theta_crit = 1;
iter_step = 10; 
beta_1 = 10;
dtime = 1;

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
load 1A_L_30cc_prefill_volume.mat;
load 1A_L_30cc_postfill_volume.mat;
load 1A_L_30cc_presac_volume.mat;
load 1A_L_30cc_postsac_volume.mat;

disp('Storing pixel coordinates...');
Nx = 100;
Ny = 100;
Nz = 100;
sizeImage = [Nx,Ny,Nz];

[pixX, pixY, pixZ]= ndgrid(linspace(1,Nx,Nx),linspace(1,Ny,Ny),linspace(1,Nz,Nz));
X = pixX(:);
Y = pixY(:);
Z = pixZ(:);
pix = [X,Y,Z];

[Xvtk,Yvtk,Zvtk] =  meshgrid(1:Nx,1:Ny,1:Nz);

pxx = Yvtk;
pyy = Xvtk;
pzz = Zvtk;

load 'deformed_pts_postfill.mat'
pxx_1 = pxx;
pyy_1 = pyy;
pzz_1 = pzz;

for nsteps = 2:2
    
    fprintf('Nsteps = %i\n',nsteps);
    
    disp('Reading image data...');
    
    Img_source = g1A_postfill_img;
    if(nsteps==1)
        Img_target = g1A_postfill_img;
    elseif(nsteps==2)
        Img_target = g1A_presac_img;
    else
        Img_target = g1A_postsac_img;
    end
    
    Img_source1 = resize3Dmatrix(100,100,100,Img_source);
    Img_target1 = resize3Dmatrix(100,100,100,Img_target);
    
    Img_source1 = (Img_source1-min(Img_source1(:)))/(max(Img_source1(:))-min(Img_source1(:)));
    Img_target1 = (Img_target1-min(Img_target1(:)))/(max(Img_target1(:))-min(Img_target1(:)));
    
    Img_source = Img_source1;
    Img_target = Img_target1;
    
    Img_source_initial = Img_source; %store initial moving image
    Img_target_initial = Img_target; %store initial target image
    
    [Gx_img,Gy_img,Gz_img] = extract_normal(Img_source);
    
    %compute initial RS value
    disp('Compute initial RS...');
    Idiff = Img_source-Img_target;
    Idiff2 = (Idiff).^2;
    residual_initial = sqrt(sum(sum(sum(Idiff2,3),2),1));
    RS_initial = 1;
    disp(sprintf('initial residual %f\n',residual_initial));
    disp(sprintf('initial RS value %f\n',RS_initial));

    param = setparameters_skin(3,alpha(nsteps));
    
    filename1 = sprintf('post_processing/source_%d.vtk',nsteps);
    filename2 = sprintf('post_processing/target_%d.vtk',nsteps);
    vtkwrite(filename1, 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_source);
    vtkwrite(filename2, 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_target);
    
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
        flag_mass =  1;
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
        
        if(updateFlag==1)
            [ACP, Pm, Node] = updateControlPoints(ACP, Pm, pxx_1, pyy_1, pzz_1,bf, bf_ct, Node, ActiveNodes);
        end
        
        %store undeformed control points
        
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
        [Gx_gp,Gy_gp,Gz_gp] = extract_normal_gp(cII_SX,ac_ct,param.orderGauss);
        
        if(updateFlag==1)
            theta_g_new = ones(ac_ct,param.orderGauss,param.orderGauss,param.orderGauss);
            theta_g_old = theta_g_new;
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
        RS_final = 2;
        level_i=150;
        
        filename = sprintf('post_processing/skin_mesh_undeformed_%d.vtk',nsteps);
        PlotGrid(filename,ActiveNodes,Node,pxx,pyy,pzz,Em,ac_ct,ac);
        
        if((multilev+1) == param.maxlevel)
            itermax=180;
        end
        
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
            cII_TY = interp3(pixY, pixX, pixZ, Img_target, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            
            cDII_TY_X = interp3(pixY, pixX, pixZ, DIITX,BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            cDII_TY_Y = interp3(pixY, pixX, pixZ, DIITY, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            cDII_TY_Z = interp3(pixY, pixX, pixZ, DIITZ, BIGYY, BIGXX, BIGZZ,'*linear',min(Img_target(:)));
            
            % denominator of the fidelity term (g(x))
            %denominate = sqrt((cDII_TY_X.^2) + (cDII_TY_Y.^2) + (cDII_TY_Z.^2)+ smallNumber); %g(x)
            
            % fidelity term in x, y, z directions
            Bterm1 = (cII_TY - cII_SX).*2.*cDII_TY_Y; %./denominate;
            Bterm2 = (cII_TY - cII_SX).*2.*cDII_TY_X; %./denominate;
            Bterm3 = (cII_TY - cII_SX).*2.*cDII_TY_Z; %./denominate;
            
            %Bterm = [Bterm1;Bterm2;Bterm3];
            img_error = (cII_TY-cII_SX).^2;
            
            % Now compute the integrals for each basis function in the support of the active cell.
            if(mod(iterct_level,iter_step)==0)
                if(flag_mass==1)
                    [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error, theta_g_new, Mass_matrix, RHS_theta1] = compute_Integ_Domain_skin_growth_projection_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,10, Wu,Wv,Ww,H,cII_SX,kappa,dtime,theta_g_old,theta_crit,flag_mass,Gx_gp,Gy_gp,Gz_gp);
                    flag_mass = 0;
                else
                    [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error, theta_g_new, ~, RHS_theta1] = compute_Integ_Domain_skin_growth_projection_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,10, Wu,Wv,Ww,H,cII_SX,kappa,dtime,theta_g_old,theta_crit,flag_mass,Gx_gp,Gy_gp,Gz_gp);
                end
                RHS_theta = RHS_theta1(:,1);
                theta_g_new_cp = Mass_matrix\RHS_theta;
                theta_g_old = theta_g_new;
            else
                [RHS, RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error] = compute_Integ_Domain_hyper_elastic_skin_theta_mex(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,param.mu, param.lambda, param.alpha, param.beta,beta_1,Wu,Wv,Ww,H,cII_SX,theta_g_old,Gx_gp,Gy_gp,Gz_gp);
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
            
            if(mod(iterct_level,10)==0)
                param.alpha = param.alpha+0.25;
            end
            
            if(mod(iterct_level,iter_step)==0)
                %compute the new positions of the pixel coordiantes using the spatial
                %transformation function
                [pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy,area_change,theta_g,theta_e] = tripleIterLoop_skin_growth_mex(sizeImage, Pixel, Jm, ACP, Img_source,param.mu,param.lambda,Gx_img,Gy_img,Gz_img,theta_g_new_cp);
                %[pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy, area_change] = tripleIterLoop_skin_mex(sizeImage, Pixel, Jm, ACP, Img_source,param.mu,param.lambda,Gx,Gy,Gz);
                
                Img_target_reparameterized_at_pixels = interp3(pixY, pixX, pixZ, Img_target, pyy, pxx, pzz,'*linear',min(Img_target(:)));
                
                filename_vtk = sprintf('post_processing/evolve_skin_%d.vtk',iterct);
                vtkwrite(filename_vtk, 'structured_grid',pyy,pxx,pzz,'scalars','StrainEnergy',strain_energy,'scalars','J',Jdet,'scalars','F11',F11,'scalars','F22',F22,'scalars','F33',F33,'scalars','Area_change',area_change,'scalars','Growth',theta_g,'scalars','Theta_e',theta_e,'scalars','Intensity',Img_source);
            end
            
            disp(sprintf('Iteration = %d, residual = %f, hyperelastic = %f',iterct_level, final_error,norm_reg));
            rs(iterct,1) = final_error;
            iter_reg(iterct,1) = norm_reg;
            iter_img(iterct,1) = norm_img;
            iterations(iterct,1) = iterct;
            theta_plot_avg(iterct,1) = mean(theta_g_new(:));
            theta_plot_max(iterct,1) = max(theta_g_new(:));
            theta_plot_avg_p(iterct,1) = mean(theta_g(:));
            theta_plot_max_p(iterct,1) = max(theta_g(:));
            theta_plot_center(iterct,1) = theta_g(round(Nx/2),round(Ny/2),round(Nz/2));
            
            subplot(3,2,1)
            plot(iterations, rs);
            subplot(3,2,2)
            plot(iterations, iter_reg);
            subplot(3,2,3)
            plot(iterations,theta_plot_avg);
            subplot(3,2,4)
            plot(iterations,theta_plot_avg_p);
            subplot(3,2,5)
            plot(iterations,theta_plot_max);
            subplot(3,2,6)
            plot(iterations,theta_plot_max_p);
            drawnow
            
            save('post_processing/skin_residuals.mat','iter_reg','iter_img','iterations');
            save('post_processing/thetag_result_presac.mat','theta_g_new','theta_plot_avg','theta_plot_max','theta_plot_avg_p','theta_plot_max_p','theta_plot_center','theta_g');
            
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
        [CellGrad, meanGrad] = computeDiffGradImage(cell_co,Img_source,Img_target, pixX, pixY, pixZ);
        fprintf('Mean Gradient = %f\n',meanGrad);
        
        updateFlag=0;
    end
    
    %itermax = itermax+500;
    updateFlag=1;
    
    filename_vtk1 = sprintf('post_processing/skin_deformed_%d.vtk',nsteps);
    vtkwrite(filename_vtk1, 'structured_grid',pyy,pxx,pzz,'scalars','Intensity',Img_source);
    
    filename1 = sprintf('post_processing/Skin_mesh_%d.vtk',nsteps);
    PlotGrid(filename1,ActiveNodes,Node,pxx,pyy,pzz,Em,ac_ct,ac);
end

diary off
