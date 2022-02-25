function [Gx_final,Gy_final,Gz_final] = extract_normal(Img)

[Gx,Gy,Gz] = imgradientxyz(Img);

Gx = -Gx./22;
Gy = -Gy./22;
Gz = -Gz./22;

for i=1:size(Gx,1)
    for j=1:size(Gx,2)
        for k=1:size(Gx,3)
            if(Gz(i,j,k)<0)
                Gz(i,j,k) = -Gz(i,j,k);
                Gx(i,j,k) = -Gx(i,j,k);
                Gy(i,j,k) = -Gy(i,j,k);
            end
        end
    end
end

sigma = 1.5;
Gx = imgaussfilt3(Gx,sigma,'FilterSize',5);
Gy = imgaussfilt3(Gy,sigma,'FilterSize',5);
Gz = imgaussfilt3(Gz,sigma,'FilterSize',5);

for i=1:size(Gx,1)
    for j=1:size(Gx,2)
        for k=1:size(Gx,3)
            n_len = sqrt(Gx(i,j,k)^2+ Gy(i,j,k)^2 +  Gz(i,j,k)^2);
            Gx(i,j,k) = Gx(i,j,k)/n_len;
            Gy(i,j,k) = Gy(i,j,k)/n_len;
            Gz(i,j,k) = Gz(i,j,k)/n_len;
        end
    end
end

Gxx = zeros(size(Img));
Gyy = zeros(size(Img));
Gzz = zeros(size(Img));

for ig = 1:size(Img,1)
    for  jg =1:size(Img,2)
        for kg = 1:size(Img,3)
            if(Img(ig,jg,kg)>0)
                Gxx(ig,jg,kg) = Gx(ig,jg,kg);
                Gyy(ig,jg,kg) = Gy(ig,jg,kg);
                Gzz(ig,jg,kg) = Gz(ig,jg,kg);
            end
        end
    end
end

Gx_final = Gxx;
Gy_final = Gyy;
Gz_final = Gzz;

end