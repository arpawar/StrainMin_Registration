%% Nearest neighbour interpolation
function phi_cp = compute_phi_cp(Pm,Img)
for i =1:size(Pm(1,1).pts,1)
    
    phi_cp(i,1) = Img(round(Pm(1,1).pts(i,1)), round(Pm(1,1).pts(i,2)), round(Pm(1,1).pts(i,3)));
    
end
end