function [pxx, pyy, pzz,F11, F22, F33, Jdet, strain_energy, area_change,theta_g_final,theta_e_final] = tripleIterLoop_skin_growth(sizeImage, Pixel, Jm, ACP, Img_source,mu, lambda,Gx,Gy,Gz,theta_g_cp)
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

Jdet = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
strain_energy = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

area_change = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
theta_g_final = ones(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
theta_e_final = ones(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

parfor i = 1:sizeImage(1,3)
    
    [temp_pxx, temp_pyy, temp_pzz ,temp_F11, temp_F22, temp_F33, temp_Jdet, temp_strain_energy, temp_area_change, temp_theta_g, temp_theta_e] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP, Img_source(:,:,i),mu, lambda,Gx(:,:,i),Gy(:,:,i),Gz(:,:,i), theta_g_cp);
    
    pxx(:,:,i) = temp_pxx;
    pyy(:,:,i) = temp_pyy;
    pzz(:,:,i) = temp_pzz;
    
    F11(:,:,i) = temp_F11;
    F22(:,:,i) = temp_F22;
    F33(:,:,i) = temp_F33;
    
    Jdet(:,:,i) = temp_Jdet;
    strain_energy(:,:,i) = temp_strain_energy;
    area_change(:,:,i) = temp_area_change;
    theta_g_final(:,:,i) = temp_theta_g;
    theta_e_final(:,:,i) = temp_theta_e; 
end

end

function [ pxx, pyy, pzz, F11, F22, F33, Jdet, strain_energy, area_change_1, theta_g_1, theta_e_1] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP, Img_source,mu,lambda,Gx,Gy,Gz, theta_g_cp)

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

Jdet = zeros(sizeImage(1,1),sizeImage(1,2));
strain_energy = zeros(sizeImage(1,1),sizeImage(1,2));
area_change_1 = zeros(sizeImage(1,1),sizeImage(1,2));
theta_g_1 = ones(sizeImage(1,1),sizeImage(1,2));
theta_e_1 = ones(sizeImage(1,1),sizeImage(1,2));

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
        theta_cp = theta_g_cp(SB,1);
        
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
        
        if(Img_source(k,j,1) > 0)
            F11(k,j) = pxxu(k,j);
            F22(k,j) = pyyv(k,j);
            F33(k,j) = pzzw(k,j);
            
            F = [pxxu(k,j), pxxv(k,j), pxxw(k,j); pyyu(k,j), pyyv(k,j), pyyw(k,j); pzzu(k,j), pzzv(k,j), pzzw(k,j)];
            B = F*transpose(F);
            Jdet(k,j) = det(F);
            strain_energy(k,j) = (mu/2)*(trace(B)-3)-mu*log(Jdet(k,j))+lambda/2*log(Jdet(k,j))^2;
            C = transpose(F)*F;
            C_inv = inv(C);
            
            norml = [Gx(k,j,1);Gy(k,j,1);Gz(k,j,1)];
            val = Jdet(k,j)*sqrt(norml'*C_inv*norml);
            area_change_1(k,j) = val;
            theta_g_1(k,j) = theta_cp'*supp;
            theta_e_1(k,j) = area_change_1(k,j)/theta_g_1(k,j);
        end
    end
end

end