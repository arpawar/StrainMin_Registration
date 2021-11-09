function phi_interpolated = interpolateBsplineFunc(Nx,Ny,Nz,Pixel,Jm,Pm,Dm,phi_cp)

px = 0;
phi_interpolated = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        for kk= 1:Nz
            px = px +1;
            ac_ind = Pixel(px).active_cell;
            supp = Pixel(px).phi;
            
            SB = Jm(ac_ind,1).nzsplines;
            ss = size(SB,1);
            f = 0;
            
            for k = 1:ss
                CEb = Pm(1,1).pts;
                BB = Dm(1,1).actB;
                BB_active = BB(SB(k,1),1);
                pio = CEb(SB(k,1),1);
                pj = CEb(SB(k,1),2);
                f = f + phi_cp(BB_active,1)*supp(k,1);
            end
            phi_interpolated(j,i,kk) = f;
        end
    end
end

end