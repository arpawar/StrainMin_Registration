function [RHS_final,RHS_image, RHS_reg, LRHS_reg, RHS_total, final_error] = compute_Integ_Domain_hyper_elastic_skin_theta(Jm,img_error,Bterm1,Bterm2,Bterm3,BIGXX, BIGYY, BIGZZ, BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,mu, lambda, alpha, beta,beta_1, w1,w2,w3,H,Img_source,theta_g_old,Gx_ip,Gy_ip,Gz_ip)
%#codegen
% This function computes the integral of the energy functional to update
% the position of the control points

% INPUT:
% Jm: the non-zero splines over the active elements
% Bterm1, Bterm2, Bterm3: the fidelity term in x, y, z direction
% respectively
% BIGMUX, BIGMUY, BIGMUZ: f_u(x) computed at the gauss points in the control
% grid
% BIGMVX, BIGMVY, BIGMVZ: f_v(x) computed at the gauss points in the control
% grid
% BIGMWX, BIGMWY, BIGMWZ: f_w(x) computed at the gauss points in the control
% grid
% RHS: the right hand side of the computation of the update of the control
% points
% PHI1,PHIU1,PHIV1,PHIW1:  the phi and first derivative of the phi in x, y,
% z at the Gauss points in the active elements
% lambda1, lambda2: the regularization parameters
% w1,w2,w3: the weights of the gaussian quadrature
% H: array containing the size of the active elements divided by 2

% OUTPUT:
% RHS_final: update of the right hand side of the equation after the
% computation of delta(E).

%%
ac_ct = size(Jm,1);
bf_ct = size(RHS,1);
xlen = size(Bterm1,2);

alpha

RHS_image1 = zeros(bf_ct,4);
RHS_image1(:,4) = RHS(:,4);
RHS_reg1 = RHS_image1;
LRHS_reg1 = RHS_image1;
RHS_total1 = RHS_image1;

error = 0;
sigma = 1.5;

parfor i=1:ac_ct
    RHS1 = zeros(bf_ct,4);
    RHS2 = zeros(bf_ct,4);
    RHS3 = zeros(bf_ct,4);
    RHS4 = zeros(bf_ct,4);
    
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    
    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);
    
    term_img = img_error(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term1 = Bterm1(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term2 = Bterm2(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term3 = Bterm3(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
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
    
    img_source1 = Img_source(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    Gx = Gx_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    Gy = Gy_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    Gz = Gz_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    supp_size = size(supp_phi,1);
    val1 = zeros(supp_size,1);
    val2 = zeros(supp_size,1);
    val3 = zeros(supp_size,1);
    
    val11 = zeros(supp_size,1);
    val22 = zeros(supp_size,1);
    val33 = zeros(supp_size,1);
    
    val11l = zeros(supp_size,1);
    val22l = zeros(supp_size,1);
    val33l = zeros(supp_size,1);
    
    val11r = zeros(supp_size,1);
    val22r = zeros(supp_size,1);
    val33r = zeros(supp_size,1);
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen
                theta_g_o = theta_g_old(i,gg1,gg2,gg3);
                
                gxx = gx(gg1,gg2,gg3);
                gyy = gy(gg1,gg2,gg3);
                gzz = gz(gg1,gg2,gg3);
                
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);
                
                n11 = Gx(gg1,gg2,gg3)*Gx(gg1,gg2,gg3);
                n12 = Gx(gg1,gg2,gg3)*Gy(gg1,gg2,gg3);
                n13 = Gx(gg1,gg2,gg3)*Gz(gg1,gg2,gg3);
                n23 = Gy(gg1,gg2,gg3)*Gz(gg1,gg2,gg3);
                n22 = Gy(gg1,gg2,gg3)*Gy(gg1,gg2,gg3);
                n33 = Gz(gg1,gg2,gg3)*Gz(gg1,gg2,gg3);
                norm_surf = [n11,n12,n13;n12,n22,n23;n13,n23,n33];
                
                %%deformation gradient
                if(img_source1(gg1,gg2,gg3)>=0.5)
                    Fg_inv = (1/sqrt(theta_g_o))*eye(3,3) + (1- (1/sqrt(theta_g_o)))*norm_surf;
                    F = double([g11(gg1,gg2,gg3), g21(gg1,gg2,gg3), g31(gg1,gg2,gg3); g12(gg1,gg2,gg3), g22(gg1,gg2,gg3), g32(gg1,gg2,gg3); g13(gg1,gg2,gg3), g23(gg1,gg2,gg3), g33(gg1,gg2,gg3)]);
                    Fe = F*Fg_inv;
                else
                    Fg_inv = eye(3,3);
                    Fe = eye(3,3);
                    F = eye(3,3);
                end
                
                
                %% Cauchy-Green tensor
                Ce = transpose(Fe)*Fe;
                Ce_inv = inv(Ce);
                B = F*transpose(F);
                J = det(F);
                
                Je = det(Fe);
                
                %%second Piola-Kirchoff stress tensor for Neo-Hookean
                S_e = mu*(eye(3,3)-Ce_inv) + lambda*log(Je)*Ce_inv;
                S = Fg_inv*S_e*transpose(Fg_inv);
                error = error + double(w1(gg1,1).*w2(gg2,1).*w3(gg3,1)*term_img(gg1,gg2,gg3).*hu.*hv.*hw);
                
                for j = 1:size(phi_ui,1)
                    
                    ecoord_x = [1,0,0];
                    ecoord_y = [0,1,0];
                    ecoord_z = [0,0,1];
                    
                    % F = gi dyad Gci
                    % gi = SUM phi_u ctrl_pts
                    % g1 = partial (x,y,z)/partial X
                    % g2 = partial (x,y,z)/partial Y
                    % g3 = partial (x,y,z)/partial Z
                    % dg1 = partial (dx,0,0)/partial X = [dphidU,0,0]
                    % dg2 = partial (0,dy,0)/partial Y = [0,dphidV,0]
                    % dg3 = partial (0,0,dz)/partial X =[0,0,dphidW]
                    % G1 = partial (X,Y,Z)/partial X = [1,0,0]
                    % G2 = partial (X,Y,Z)/partial Y = [0,1,0]
                    % G3 = partial (X,Y,Z)/partial Z = [0,0,1]
                    % delta_gi[incr u] =  SUM phi_u*incr u
                    deltaFx = [phi_ui(j,1)*ecoord_x' , phi_vi(j,1)*ecoord_x', phi_wi(j,1)*ecoord_x'];
                    deltaFy = [phi_ui(j,1)*ecoord_y' , phi_vi(j,1)*ecoord_y', phi_wi(j,1)*ecoord_y'];
                    deltaFz = [phi_ui(j,1)*ecoord_z' , phi_vi(j,1)*ecoord_z', phi_wi(j,1)*ecoord_z'];
                    
                    deltaEx = 0.5*(transpose(F)*deltaFx+transpose(deltaFx)*F);
                    deltaEy = 0.5*(transpose(F)*deltaFy+transpose(deltaFy)*F);
                    deltaEz = 0.5*(transpose(F)*deltaFz+transpose(deltaFz)*F);
                    
                    %delta_Fe = [phi_ui(j,1); phi_vi(j,1); phi_wi(j,1)]; %correct this
                    
                    %delta_E = (Fe*delta_Fe);
                    
                    finalx = S.*deltaEx;
                    sum_finalx = sum(finalx(:));
                    
                    finaly = S.*deltaEy;
                    sum_finaly = sum(finaly(:));
                    
                    finalz = S.*deltaEz;
                    sum_finalz = sum(finalz(:));
                    
                    %.*phi_source(floor(gxx)+1,floor(gyy)+1,floor(gzz)+1)
                    %val1(j,1) = val1(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term1(gg1,gg2,gg3)).*hu.*hv.*hw+2.*lambda_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalx.*hu.*hv.*hw;
                    %val2(j,1) = val2(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term2(gg1,gg2,gg3)).*hu.*hv.*hw+2.*lambda_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finaly.*hu.*hv.*hw;
                    %val3(j,1) = val3(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term3(gg1,gg2,gg3)).*hu.*hv.*hw+2.*lambda_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalz.*hu.*hv.*hw;
                    % AdrianQ: I think the 'terms' in this piece of code
                    % have to be post multiplied by the deformation
                    % gradient, then by ecoord, I still dont know what the
                    % scaling by the magnitude of the gradient is doing
                    % vector with the gradient of target
                    %dTY = [term1(gg1,gg2,gg3);term2(gg1,gg2,gg3);term3(gg1,gg2,gg3)];
                    %dTYFdx = dTY'*F*ecoord_x;
                    %dTYFdy = dTY'*F*ecoord_y;
                    %dTYFdz = dTY'*F*ecoord_z;
                    
                    val1(j,1) = val1(j,1) + alpha*(phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term1(gg1,gg2,gg3)).*hu.*hv.*hw+2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalx.*hu.*hv.*hw + 2.*beta_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(1-img_source1(gg1,gg2,gg3)).*sum_finalx.*hu.*hv.*hw;
                    val2(j,1) = val2(j,1) + alpha*(phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term2(gg1,gg2,gg3)).*hu.*hv.*hw+2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finaly.*hu.*hv.*hw + 2.*beta_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(1-img_source1(gg1,gg2,gg3)).*sum_finaly.*hu.*hv.*hw;
                    val3(j,1) = val3(j,1) + alpha*(phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term3(gg1,gg2,gg3)).*hu.*hv.*hw+2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalz.*hu.*hv.*hw + 2.*beta_1.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(1-img_source1(gg1,gg2,gg3)).*sum_finalz.*hu.*hv.*hw;
                    
                    val11(j,1) = val11(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term1(gg1,gg2,gg3)).*hu.*hv.*hw;
                    val22(j,1) = val22(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term2(gg1,gg2,gg3)).*hu.*hv.*hw;
                    val33(j,1) = val33(j,1) + (phi_i(j,1)).*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*(term3(gg1,gg2,gg3)).*hu.*hv.*hw;
                    
                    val11l(j,1) = val11l(j,1) + 2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalx.*hu.*hv.*hw;
                    val22l(j,1) = val22l(j,1) + 2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finaly.*hu.*hv.*hw;
                    val33l(j,1) = val33l(j,1) + 2.*beta.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalz.*hu.*hv.*hw;
                    
                    val11r(j,1) = val11r(j,1) + 2.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalx.*hu.*hv.*hw;
                    val22r(j,1) = val22r(j,1) + 2.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finaly.*hu.*hv.*hw;
                    val33r(j,1) = val33r(j,1) + 2.*w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*img_source1(gg1,gg2,gg3).*sum_finalz.*hu.*hv.*hw;
                end
            end
        end
    end
    
    RHS1(SB,1) =  val1;
    RHS1(SB,2) =  val2;
    RHS1(SB,3) =  val3;
    
    RHS2(SB,1) =  val11;
    RHS2(SB,2) =  val22;
    RHS2(SB,3) =  val33;
    
    RHS3(SB,1) =  val11r;
    RHS3(SB,2) =  val22r;
    RHS3(SB,3) =  val33r;
    
    RHS4(SB,1) =  val11l;
    RHS4(SB,2) =  val22l;
    RHS4(SB,3) =  val33l;
    
    RHS = RHS + RHS1;
    RHS_image1 = RHS_image1 + RHS2;
    RHS_reg1 = RHS_reg1 + RHS3;
    LRHS_reg1 = LRHS_reg1 + RHS4;
end

RHS_final = RHS;
RHS_image = RHS_image1;
RHS_reg = RHS_reg1;
LRHS_reg = LRHS_reg1;
RHS_total = RHS;
final_error = error;
end