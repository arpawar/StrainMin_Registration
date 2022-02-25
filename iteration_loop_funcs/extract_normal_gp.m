function [Gx_final_ip,Gy_final_ip,Gz_final_ip] = extract_normal_gp(Img,ac_ct,xlen)

Gx_final_ip = zeros(size(Img));
Gy_final_ip = zeros(size(Img));
Gz_final_ip = zeros(size(Img));
sigma = 1.5;

for i=1:ac_ct
    
    img_source1 = Img(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    [Gx,Gy,Gz] = imgradientxyz(img_source1);
    
    Gx = -Gx./22;
    Gy = -Gy./22;
    Gz = -Gz./22;
    
    for ig = 1:size(img_source1,1)
        for  jg =1:size(img_source1,2)
            for kg = 1:size(img_source1,3)
                if(Gz(ig,jg,kg)<0)
                    Gx(ig,jg,kg) = -Gx(ig,jg,kg);
                    Gy(ig,jg,kg) = -Gy(ig,jg,kg);
                    Gz(ig,jg,kg) = -Gz(ig,jg,kg);
                end
            end
        end
    end
    
    Gx = imgaussfilt3(Gx,sigma,'FilterSize',5);
    Gy = imgaussfilt3(Gy,sigma,'FilterSize',5);
    Gz = imgaussfilt3(Gz,sigma,'FilterSize',5);
    
    for ig = 1:size(img_source1,1)
        for  jg =1:size(img_source1,2)
            for kg = 1:size(img_source1,3)
                nlen = sqrt(Gx(ig,jg,kg)^2+Gy(ig,jg,kg)^2+Gz(ig,jg,kg)^2);
                Gx(ig,jg,kg) = Gx(ig,jg,kg)/nlen;
                Gy(ig,jg,kg) = Gy(ig,jg,kg)/nlen;
                Gz(ig,jg,kg) = Gz(ig,jg,kg)/nlen;
            end
        end
    end
    
    Gxx = zeros(size(img_source1));
    Gyy = zeros(size(img_source1));
    Gzz = zeros(size(img_source1));
    
    for ig = 1:size(img_source1,1)
        for  jg =1:size(img_source1,2)
            for kg = 1:size(img_source1,3)
                if(img_source1(ig,jg,kg)>0)
                    Gxx(ig,jg,kg) = Gx(ig,jg,kg);
                    Gyy(ig,jg,kg) = Gy(ig,jg,kg);
                    Gzz(ig,jg,kg) = Gz(ig,jg,kg);
                end
            end
        end
    end
    
    Gx = Gxx;
    Gy = Gyy;
    Gz = Gzz;
    
    Gx_final_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen) = Gx;
    Gy_final_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen) = Gy;
    Gz_final_ip(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen) = Gz;
    
end