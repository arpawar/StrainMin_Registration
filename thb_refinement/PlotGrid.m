function [] = PlotGrid(filename,ActiveNodes,Node,pxx,pyy,pzz,Em,ac_ct,ac)


fileID = fopen(filename,'w');
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, 'mesh\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n\n');
fprintf(fileID, 'POINTS %d float\n',size(ActiveNodes,1));

disp('reached here')
for i = 1:size(ActiveNodes,1)
    index11 = max(round(Node(ActiveNodes(i,1),1)),1);
    index22 = max(round(Node(ActiveNodes(i,1),2)),1);
    index33 = max(round(Node(ActiveNodes(i,1),3)),1);
    
    index1 = min(index11,size(pxx ,1));
    index2 = min(index22,size(pxx ,2));
    index3 = min(index33,size(pxx ,3));
    
    xvalue = pxx(index1,index2,index3);
    yvalue = pyy(index1,index2,index3);
    zvalue = pzz(index1,index2,index3);
    
    fprintf(fileID, '%f %f %f\n', xvalue, yvalue, zvalue);
end

fprintf(fileID, 'CELLS  %d %d\n',ac_ct,9*ac_ct);

for i = 1:ac_ct,
    nodes = Em(ac(i,2)).nodes(ac(i,1),:);
    fprintf(fileID, '8 %d %d %d %d %d %d %d %d\n',Node(nodes(1,1),4)-1,Node(nodes(1,2),4)-1,Node(nodes(1,3),4)-1,Node(nodes(1,4),4)-1,Node(nodes(1,5),4)-1,Node(nodes(1,6),4)-1,Node(nodes(1,7),4)-1,Node(nodes(1,8),4)-1);
end

fprintf(fileID, 'CELL_TYPES %d\n',ac_ct);
for i = 1:ac_ct,
    fprintf(fileID, '12 \n');
end

fclose(fileID);

end