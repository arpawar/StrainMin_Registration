function [pxx, pyy, pzz,F11, F22, F33, Jdet, strain_energy,stress_x, stress_y, stress_z, tstress_x, tstress_y, tstress_z, Fe_x, Fe_y, Fe_z, Fp_x, Fp_y, Fp_z] = tripleIterLoop_stress(sizeImage, Pixel, Jm, ACP, Img_source,mu, lambda, F_prestrain)
%In this function we compute the new positions of the pixel coordinates
%using the spatial transformation function

% INPUT:
% sizeImage: size of the image
% Pixel: the array containing the active element indices and the phi values
% for each pixel
% Jm: the array containing the non-zero splines in each active element
% ACP: the array containing the control points that are set to be active

% OUTPUT:
% pxx, pyy, pzz: the position of the pixel coordinates after spatial transformation

%%
pxx = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
pyy = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
pzz = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

F11 = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
F22 = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
F33 = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

Fe_x = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
Fe_y = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
Fe_z = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

Fp_x = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
Fp_y = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
Fp_z = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

Jdet = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
strain_energy = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

stress_x = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
stress_y = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
stress_z = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

tstress_x = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
tstress_y = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
tstress_z = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

parfor i = 1:sizeImage(1,3)
    
    [temp_pxx, temp_pyy, temp_pzz ,temp_F11, temp_F22, temp_F33, temp_Jdet, temp_strain_energy, temp_stress_x, temp_stress_y, temp_stress_z,temp_tstress_x, temp_tstress_y, temp_tstress_z, temp_Fe_x, temp_Fe_y, temp_Fe_z,temp_Fp_x, temp_Fp_y, temp_Fp_z] = tripleIterLoopBody_stress(i, sizeImage, Pixel, Jm, ACP, Img_source(:,:,i),mu, lambda,F_prestrain);
    
    pxx(:,:,i) = temp_pxx;
    pyy(:,:,i) = temp_pyy;
    pzz(:,:,i) = temp_pzz;
    
    F11(:,:,i) = temp_F11;
    F22(:,:,i) = temp_F22;
    F33(:,:,i) = temp_F33;
    
    Fe_x(:,:,i) = temp_Fe_x;
    Fe_y(:,:,i) = temp_Fe_y;
    Fe_z(:,:,i) = temp_Fe_z;
    
    Fp_x(:,:,i) = temp_Fp_x;
    Fp_y(:,:,i) = temp_Fp_y;
    Fp_z(:,:,i) = temp_Fp_z;
    
    Jdet(:,:,i) = temp_Jdet;
    strain_energy(:,:,i) = temp_strain_energy;
    stress_x(:,:,i) = temp_stress_x;
    stress_y(:,:,i) = temp_stress_y;
    stress_z(:,:,i) = temp_stress_z;
    
    tstress_x(:,:,i) = temp_tstress_x;
    tstress_y(:,:,i) = temp_tstress_y;
    tstress_z(:,:,i) = temp_tstress_z;
end

end

function [ pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy, stress_x, stress_y, stress_z,tstress_x, tstress_y, tstress_z, Fe_x, Fe_y, Fe_z, Fp_x, Fp_y, Fp_z] = tripleIterLoopBody_stress(i, sizeImage, Pixel, Jm, ACP, Img_source,mu,lambda,F_prestrain)

pxx = zeros(sizeImage(1,1),sizeImage(1,2));
pyy = zeros(sizeImage(1,1),sizeImage(1,2));
pzz = zeros(sizeImage(1,1),sizeImage(1,2));

pxxu = zeros(sizeImage(1,1),sizeImage(1,2));
pyyu = zeros(sizeImage(1,1),sizeImage(1,2));
pzzu = zeros(sizeImage(1,1),sizeImage(1,2));

pxxv = zeros(sizeImage(1,1),sizeImage(1,2));
pyyv = zeros(sizeImage(1,1),sizeImage(1,2));
pzzv = zeros(sizeImage(1,1),sizeImage(1,2));

pxxw = zeros(sizeImage(1,1),sizeImage(1,2));
pyyw = zeros(sizeImage(1,1),sizeImage(1,2));
pzzw = zeros(sizeImage(1,1),sizeImage(1,2));

F11 = zeros(sizeImage(1,1),sizeImage(1,2));
F22 = zeros(sizeImage(1,1),sizeImage(1,2));
F33 = zeros(sizeImage(1,1),sizeImage(1,2));

Fe_x = zeros(sizeImage(1,1),sizeImage(1,2));
Fe_y = zeros(sizeImage(1,1),sizeImage(1,2));
Fe_z = zeros(sizeImage(1,1),sizeImage(1,2));

Fp_x = zeros(sizeImage(1,1),sizeImage(1,2));
Fp_y = zeros(sizeImage(1,1),sizeImage(1,2));
Fp_z = zeros(sizeImage(1,1),sizeImage(1,2));

Jdet = zeros(sizeImage(1,1),sizeImage(1,2));
strain_energy = zeros(sizeImage(1,1),sizeImage(1,2));
stress_x = zeros(sizeImage(1,1),sizeImage(1,2));
stress_y = zeros(sizeImage(1,1),sizeImage(1,2));
stress_z = zeros(sizeImage(1,1),sizeImage(1,2));

tstress_x = zeros(sizeImage(1,1),sizeImage(1,2));
tstress_y = zeros(sizeImage(1,1),sizeImage(1,2));
tstress_z = zeros(sizeImage(1,1),sizeImage(1,2));

for j = 1:sizeImage(1,2)
    for k = 1:sizeImage(1,1)
        
        % global index of the pixel
        px = sizeImage(1,1)*sizeImage(1,2)*(i-1)+sizeImage(1,1)*(j-1)+k;
        
        % the active element conatining the pixel
        ac_ind = Pixel(px,1).active_cell;
        
        %phi value associated with the pixel coordinate
        supp = Pixel(px,1).phi;
        supp_u = Pixel(px,1).phiu;
        supp_v = Pixel(px,1).phiv;
        supp_w = Pixel(px,1).phiw;
        
        %the non-zero splines over the activ element containing the pixel
        SB = Jm(ac_ind,1).nzsplines;
        
        %control points over the active element
        pts = ACP(SB,1:3);
        
        %compute the new position of the pixel coordinates
        FXX = pts'*supp;
        FXX_u = pts(:,1)'*supp_u;
        FXX_v = pts(:,1)'*supp_v;
        FXX_w = pts(:,1)'*supp_w;
        
        FYY_u = pts(:,2)'*supp_u;
        FYY_v = pts(:,2)'*supp_v;
        FYY_w = pts(:,2)'*supp_w;
        
        FZZ_u = pts(:,3)'*supp_u;
        FZZ_v = pts(:,3)'*supp_v;
        FZZ_w = pts(:,3)'*supp_w;
        
        pxx(k,j) = FXX(1,1);
        pyy(k,j) = FXX(2,1);
        pzz(k,j) = FXX(3,1);
        
        pxxu(k,j) = FXX_u(1,1);
        pyyu(k,j) = FYY_u(1,1);
        pzzu(k,j) = FZZ_u(1,1);
        
        pxxv(k,j) = FXX_v(1,1);
        pyyv(k,j) = FYY_v(1,1);
        pzzv(k,j) = FZZ_v(1,1);
        
        pxxw(k,j) = FXX_w(1,1);
        pyyw(k,j) = FYY_w(1,1);
        pzzw(k,j) = FZZ_w(1,1);
        
        if(Img_source(k,j,1) > 0.9)
            Fp = F_prestrain;
            Fg = eye(3,3);
            Fg_star = Fg*Fp;
            Fg_star_inv = inv(Fg_star);
            F = [pxxu(k,j), pxxv(k,j), pxxw(k,j); pyyu(k,j), pyyv(k,j), pyyw(k,j); pzzu(k,j), pzzv(k,j), pzzw(k,j)];
            Fe = F*Fg_star_inv;
            F11(k,j) = pxxu(k,j);
            F22(k,j) = pyyv(k,j);
            F33(k,j) = pzzw(k,j);
            
            Fe_x(k,j) = Fe(1,1);
            Fe_y(k,j) = Fe(2,2);
            Fe_z(k,j) = Fe(3,3);
            
            Fp_x(k,j) = Fp(1,1);
            Fp_y(k,j) = Fp(2,2);
            Fp_z(k,j) = Fp(3,3);
            
            Ce = transpose(Fe)*Fe;
            Ce_inv = inv(Ce);
            B = F*transpose(F);
            Jdet(k,j) = det(F); 
            Je = det(Fe);
            
            %%second Piola-Kirchoff stress tensor for Neo-Hookean
            S_e = mu*(eye(3,3)-Ce_inv) + lambda*log(Je)*Ce_inv;
            S = Fg_star_inv*S_e*transpose(Fg_star_inv);
               
            strain_energy(k,j) = (mu/2)*(trace(B)-3)-mu*log(Jdet(k,j))+lambda/2*log(Jdet(k,j))^2;

            tstress_x(k,j) = S(1,1);
            tstress_y(k,j) = S(2,2);
            tstress_z(k,j) = S(3,3);
            
            stress_x(k,j) = S_e(1,1);
            stress_y(k,j) = S_e(2,2);
            stress_z(k,j) = S_e(3,3);
        end
    end
end

end

