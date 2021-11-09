function [ pxx, pyy, pzz, pxxu, pyyu, pzzu, pxxv, pyyv, pzzv, pxxw, pyyw, pzzw ] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP  )
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
for j = 1:sizeImage(1,2)
    for k = 1:sizeImage(1,1)
        px = sizeImage(1,1)*sizeImage(1,2)*(i-1)+sizeImage(1,1)*(j-1)+k;
        ac_ind = Pixel(px,1).active_cell;
        supp = Pixel(px,1).phi;
        suppu = Pixel(px,1).phiu;
        suppv = Pixel(px,1).phiv;
        suppw = Pixel(px,1).phiw;
        SB = Jm(ac_ind,1).nzsplines;
        
        pts = ACP(SB,1:3);
        FXX = pts'*supp;
        FXX_u = pts(:,1)'*suppu;
        FXX_v = pts(:,1)'*suppv;
        FXX_w = pts(:,1)'*suppw;
        
        FYY_u = pts(:,2)'*suppu;
        FYY_v = pts(:,2)'*suppv;
        FYY_w = pts(:,2)'*suppw;
        
        FZZ_u = pts(:,3)'*suppu;
        FZZ_v = pts(:,3)'*suppv;
        FZZ_w = pts(:,3)'*suppw;
        
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
    end
end
end

