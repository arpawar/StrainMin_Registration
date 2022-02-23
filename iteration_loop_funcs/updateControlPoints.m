function [ACP_new, Pm_new, Node_new] = updateControlPoints(ACP,Pm,pxx,pyy,pzz,bf, bf_ct,Node, ActiveNodes)

ACP_new1 = ACP;
Node_new1 = Node;

for i = 1:size(ACP,1)
    if(ACP(i,4)==0)
        index11 = max(round(ACP(i,1)),1);
        index22 = max(round(ACP(i,2)),1);
        index33 = max(round(ACP(i,3)),1);
        
        index1 = min(index11,size(pxx ,1));
        index2 = min(index22,size(pxx ,2));
        index3 = min(index33,size(pxx ,3));
        
        xvalue = pxx(index1,index2,index3);
        yvalue = pyy(index1,index2,index3);
        zvalue = pzz(index1,index2,index3);

        ACP_new1(i,1:3) = [xvalue,yvalue,zvalue];
    end
end
ACP_new = ACP_new1;

for i = 1:size(ActiveNodes,1)
    index11 = max(round(Node(ActiveNodes(i,1),1)),1);
    index22 = max(round(Node(ActiveNodes(i,1),2)),1);
    index33 = max(round(Node(ActiveNodes(i,1),3)),1);
    
    index1 = min(index11,size(pxx ,1));
    index2 = min(index22,size(pxx ,2));
    index3 = min(index33,size(pxx ,3));
    
    if(Node(ActiveNodes(i,1),1)~=0)
        Node_new1(ActiveNodes(i,1),1) = pxx(index1,index2,index3);
    end
    
    if(Node(ActiveNodes(i,1),2)~=0)
        Node_new1(ActiveNodes(i,1),2) = pyy(index1,index2,index3);
    end
    
    if(Node(ActiveNodes(i,1),3)~=0)
        Node_new1(ActiveNodes(i,1),3) = pzz(index1,index2,index3);
    end
end

Node_new = Node_new1;

for ilevel = 1:bf_ct
    level_b = bf(ilevel,2);
    Pmlevel = Pm(level_b).pts;
    index = bf(ilevel,1);
    Pmlevel(index,:) = ACP_new(ilevel,1:3);
    Pm(level_b).pts = Pmlevel;
end

Pm_new = Pm;

end