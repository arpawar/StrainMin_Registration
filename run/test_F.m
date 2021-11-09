clc;
clear all;
close all;

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();

%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_biaxial();

%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-04;
%maximum number of iterations for each refinement level
itermax = 200;

%% 1: Read image data
%Enter file name here
disp('Reading image data...');
smallplate_img = zeros(100,100,100);
bigplate_img = zeros(100,100,100);

%Image size
sizeImage = size(smallplate_img);
Nx = sizeImage(1,1);
Ny = sizeImage(1,2);
Nz = sizeImage(1,3);

for i =1:Nx
    for j=1:Ny
        for k=1:Nz
            if(i>=35 && i<=65 && j>=35 && j<=65 && k>=35 && k<=45)
                smallplate_img(i,j,k) = 1;
            end
            if(i>=33 && i<=67 && j>=33 && j<=67 && k>=35 && k<=45)
                bigplate_img(i,j,k) = 1;
            end
        end
    end
end

img2 = bigplate_img;
img1 = smallplate_img;

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
[CP, knotvectorU, knotvectorV, knotvectorW] = storeKnotArray(param,sizeImage);

%% Compute Elem, Basis, Control point data structure for THBS splines
disp('Set B-spline grid...');

tic
[Dm,Em] = setBsplineGrid_func_withoutNode(knotvectorU,knotvectorV,knotvectorW,param);
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
Pm_old = CP;

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

vtkwrite('source_plate.vtk', 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_source);
vtkwrite('target_plate.vtk', 'structured_grid',Xvtk,Yvtk,Zvtk,'scalars','Intensity',Img_target);

%% Compute basis function at pixel coordinates and Gauss points 


% No refinement
multilev = 0;

disp('Collecting active elements, control points and basis functions...');
[ac, bf, ACP, RHS,Em,Dm] = storeActiveElem(Em,Dm,Pm,multilev);
ac_ct = size(ac,1);
bf_ct = size(bf,1);
disp(sprintf('active elements = %d, active DOF =  %d at level = %d \n',ac_ct,bf_ct,multilev+1));
    
disp('Computing the non-zeros spline over each active element and storing coefficient matrices...');
tic
[Jm, Coeff] = computeNonZeroSplines(ac, param, Em, Dm);
toc

disp('Computing the basis functions at pixel coordinates...');
tic
numPixels = int64(prod(sizeImage));
[Pixel, Pix2] = storePixelPhi(numPixels, multilev,pix, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param);
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

[PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param);
% interpolate the intensity values of the target image at the gauss
% points stored in BIGX, BIGY, BIGZ
ImgSource_GPX = interp3(pixY, pixX, pixZ, Img_source_initial, BIGY, BIGX, BIGZ,'*linear',min(Img_target(:)));
toc;
    
%% Apply some affine deformation

nCP = size(ACP,1);
F1 = [1.2,0,0;0,1.2,0;0,0,1];
F2x = [0.75,0,0;0,1,0;0,0,1];
F2y = [1,0,0;0,0.75,0;0,0,1];
for ci=1:nCP
    ACP(ci,1:3) = F1*ACP(ci,1:3)' - [50*0.2;50*0.2;0];
end
for ci=1:nCP
    if ACP(ci,1)>69
        %disp('shrinking this control point in X')
        %disp(ACP(ci,1:3));
        ACP(ci,1:3) = F2x*ACP(ci,1:3)' + [50*0.35;0;0];
    end
    if ACP(ci,2)>69
        %disp('shrinking this control point in Y')
        %disp(ACP(ci,1:3));
        ACP(ci,1:3) = F2y*ACP(ci,1:3)' + [0;50*0.35;0];
    end
    if ACP(ci,1)<31
        %disp('shrinking this control point in X')
        %disp(ACP(ci,1:3));
        ACP(ci,1:3) = F2x*ACP(ci,1:3)' + [50*0.15;0;0];
    end
    if ACP(ci,2)<31
        %disp('shrinking this control point in Y')
        %disp(ACP(ci,1:3));
        ACP(ci,1:3) = F2y*ACP(ci,1:3)' + [0;50*0.15;0];
    end
end
% convert the cell array to struct array
PHI1 = cell2struct(PHI,'mat',2);
PHIU1 = cell2struct(PHIU,'mat',2);
PHIV1 = cell2struct(PHIV,'mat',2);
PHIW1 = cell2struct(PHIW,'mat',2);

% get the deformed positions and deformation gradient at integration points 
orderGauss = param.orderGauss;
[BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ] = computenewPoints(Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,orderGauss);
% initially this should be identity 
ac_ct = size(Jm,1);
bf_ct = size(RHS,1);
xlen = orderGauss;
for i=1:ac_ct

    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;

    gx = BIGXX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gy = BIGYY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gz = BIGZZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    g11 = BIGMUX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g12 = BIGMUY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g13 = BIGMUZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    g21 = BIGMVX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g22 = BIGMVY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g23 = BIGMVZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    g31 = BIGMWX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g32 = BIGMWY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    g33 = BIGMWZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen

                gxx = gx(gg1,gg2,gg3);
                gyy = gy(gg1,gg2,gg3);
                gzz = gz(gg1,gg2,gg3);
                
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);

                %%deformation gradient
                F = [g11(gg1,gg2,gg3), g21(gg1,gg2,gg3), g31(gg1,gg2,gg3); g12(gg1,gg2,gg3), g22(gg1,gg2,gg3), g32(gg1,gg2,gg3); g13(gg1,gg2,gg3), g23(gg1,gg2,gg3), g33(gg1,gg2,gg3)];
                %disp(F)
            end
        end
    end
    
end

%% compute F at the pixels instead of the gauss points 
%compute the new positions of the pixel coordiantes using the spatial
%transformation function
[pxx, pyy, pzz, F11, F22, Jdet, strain_energy] = tripleIterLoop(sizeImage, Pixel, Jm, ACP, Img_target, param);

% save the geometry
filename_vtk = 'Ftest_deformed.vtk';
vtkwrite(filename_vtk, 'structured_grid',pyy,pxx,pzz,'scalars','StrainEnergy',strain_energy,'scalars','J',Jdet,'scalars','F11',F11,'scalars','F22',F22);

