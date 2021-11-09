function phi_interpolated = interpolateBsplineFunc_gausspts(Nx,Ny,Nz,PHI,BIGY, BIGX, BIGZ,Jm,Pm,Dm,phi_cp)

px = 0;
phi_interpolated = zeros(size(BIGX));

for i = 1:ac_ct
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    gx = BIGXX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gy = BIGYY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    gz = BIGZZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    ss = size(SB,1);
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen
                
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                f = 0;
                
                for k = 1:ss
                    CEb = Pm(1,1).pts;
                    BB = Dm(1,1).actB;
                    BB_active = BB(SB(k,1),1);
                    pio = CEb(SB(k,1),1);
                    pj = CEb(SB(k,1),2);
                    f = f + phi_cp(BB_active,1)*phi_i(k,1);
                end
                phi_interpolated(i,gg1,gg2,gg3) = f;
            end
        end
    end
    
end