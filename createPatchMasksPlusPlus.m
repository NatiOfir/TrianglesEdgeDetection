function createPatchMasksPlusPlus(param)
    close all;
    index = 1;
    prm = getPrm(param);
    S = prm.minLength;
    w = prm.filterWidth;
    
    numOfEdges = (2+1)*S^2;
    maskTableL0 = zeros((S+2*w)^2,numOfEdges);
    maskTableL1 = zeros((S+2*w)^2,numOfEdges);
    maskTableL2 = zeros((S+2*w)^2,numOfEdges);
    maskTableL3 = zeros((S+2*w)^2,numOfEdges);
    
    edgeTableL0 = zeros((S)^2,numOfEdges);
    edgeTableL1 = zeros((S)^2,numOfEdges);
    edgeTableL2 = zeros((S)^2,numOfEdges);
    edgeTableL3 = zeros((S)^2,numOfEdges);
    
    edgeIndx = zeros(S,numOfEdges)+NaN;
    edgeIndxL1 = zeros(S,numOfEdges)+NaN;
    edgeIndxL2 = zeros(S,numOfEdges)+NaN;
    edgeIndxL3 = zeros(S,numOfEdges)+NaN;
    edgeIndy = zeros(S,numOfEdges)+NaN;
    edgeIndyL1 = zeros(S,numOfEdges)+NaN;
    edgeIndyL2 = zeros(S,numOfEdges)+NaN;
    edgeIndyL3 = zeros(S,numOfEdges)+NaN;

    metaSize = 8;
    metaTableL0 = zeros(metaSize,numOfEdges);
    metaTableL1 = zeros(metaSize,numOfEdges);
    metaTableL2 = zeros(metaSize,numOfEdges);
    metaTableL3 = zeros(metaSize,numOfEdges);
    gap = w+1;
    Snew = S+2*w;
    
    for edge1 = 0:1
        for edge2 = edge1+1:2
            inTriangle = getInTriangle(edge1,edge2,S,w);
            for v1 = 0:S-1
                for v2 = 0:S-1
                    [x0,y0] = TrianglesTree.getXYfromIndex(edge1,v1,S);
                    [x1,y1] = TrianglesTree.getXYfromIndex(edge2,v2,S);
                    ind = TrianglesTree.getIndexFromEdge(edge1,v1,edge2,v2,S,S);
                    L = max(abs(x1-x0),abs(y1-y0));
                    v = [x1-x0,y1-y0];
                    length = max(v);
                    theta = atan(v(2)/v(1));                                           
                    theta = 180*theta/pi;
                                        
                    if v(1)<0
                        theta = theta+180;
                    elseif theta<0
                        theta = theta+360;
                    end
                    if L == 0 
                        continue;
                    end
                    
                    curMask = zeros(Snew);
                    Ll = 0;
                    Lr = 0;
                    if edge1 == 0 && edge2 == 1
                        if 0%v2 == 0 || v1 == 4;
                        else
                            for k = 1:w
                                dy = k; %k*1/cost;
                                dx = k; %k*1/sint;
                                [l,maskRight] = getMask(x0+dx+gap,y0+gap,x1+gap,y1-dy+gap,Snew);
                                [r,maskLeft] = getMask(x0-dx+gap,y0+gap,x1+gap,y1+dy+gap,Snew);
                                Ll = Ll+l;
                                Lr = Lr+r;
                                curMask = curMask+maskLeft-maskRight;
                            end
                        end
                    elseif edge1 == 0 && edge2 == 2
                        if 0%v1 == 0 || v2 ==4;
                        else
                            for k = 1:w
                                dy = k; %k*1/cost;
                                dx = k; %k*1/sint;
                                dhx = k; %k*cosa/sqrt(2);
                                dhy = k; %k*cosa/sqrt(2);
                                [l,maskRight] = getMask(x0+dx+gap,y0+gap,x1+dhx+gap,y1+dhy+gap,Snew);
                                [r,maskLeft] = getMask(x0-dx+gap,y0+gap,x1-dhx+gap,y1-dhy+gap,Snew);
                                Ll = Ll+l;
                                Lr = Lr+r;
                                curMask = curMask+maskLeft-maskRight;
                            end
                        end
                    elseif edge1 == 1 && edge2 == 2
                        if 0%v1 == 4 || v2 == 0;
                        else
                            
                            for k = 1:w
                                dy = k; %k*1/cost;
                                dx = k; %k*1/sint;
                                dhx = k; %k*cosa;
                                dhy = k; %k*sina;
                                [l,maskLeft] = getMask(x0+gap,y0-dy+gap,x1-dhx+gap,y1-dhy+gap,Snew);
                                [r,maskRight] = getMask(x0+gap,y0+dy+gap,x1+dhx+gap,y1+dhy+gap,Snew);
                                Ll = Ll+l;
                                Lr = Lr+r;
                                curMask = curMask+maskLeft-maskRight;
                            end
                        end
                    end
                    
                    line1 = Line(x0+gap,y0+gap,x1+gap,y1+gap,1);
                    line1 = line1.calcPixels();
                    curEdgeMask = line1.getLineImage(Snew,Snew);
                    %[l,curEdgeMask] = getMask(x0+gap,y0+gap,x1+gap,y1+gap,Snew);
                    
                    %curMask(curEdgeMask == 1 | inTriangle~=1) = 0;
                    curMask(curEdgeMask == 1) = 0;
                    left = curMask>0.0;
                    right = curMask<0.0;
                    %left = curMask(curMask>0.0);
                    %right = curMask(curMask<0.0);
                    %bound = abs(curMask) == 0.5;
                    curMask = left-right;
                    %curMask(bound) = curMask(bound)*0.5;
                    
                    Ll = abs(sum(left(:)));
                    Lr = abs(sum(right(:)));
                                                          
                    if Lr<=0 || Ll<=0
                        curMask = zeros(S+2*w);
                        Lr = nan;
                        Ll = nan;
                    else
                        subplot(15,15,index),imshow(curMask/2+0.5);
                        index = index+1;
                    end
                    curMaskL1 = imrotate(curMask,90);
                    curMaskL2 = imrotate(curMask,180);
                    curMaskL3 = imrotate(curMask,270);
                    curMask = reshape(curMask,(S+2*w)^2,1);
                    curMaskL1 = reshape(curMaskL1,(S+2*w)^2,1);
                    curMaskL2 = reshape(curMaskL2,(S+2*w)^2,1);
                    curMaskL3 = reshape(curMaskL3,(S+2*w)^2,1);

                    line1 = Line(x0+1,y0+1,x1+1,y1+1,1);
                    line1 = line1.calcPixels();
                    curEdge = line1.getLineImage(S,S);
                    %curEdge
                    curEdgeL1 = imrotate(curEdge,90);
                    curEdgeL2 = imrotate(curEdge,180);
                    curEdgeL3 = imrotate(curEdge,270);
                    
                    [indx,indy] = find(curEdge);
                    indx = indx-1;
                    indy = indy-1;
                    %[indyL1,indxL1] = find(curEdgeL1);
                    indxL1 = -indy;
                    indyL1 = indx;
                    %[indxL2,indyL2] = find(curEdgeL2);
                    indxL2 = -indx;
                    indyL2 = -indy;
                    %[indyL3,indxL3] = find(curEdgeL3);
                    indxL3 = indy;
                    indyL3 = -indx;
                    
                    curEdge = reshape(curEdge,(S)^2,1);
                    curEdgeL1 = reshape(curEdgeL1,(S)^2,1);
                    curEdgeL2 = reshape(curEdgeL2,(S)^2,1);
                    curEdgeL3 = reshape(curEdgeL3,(S)^2,1); 

                    maskTableL0(:,ind) = curMask(:);
                    maskTableL1(:,ind) = curMaskL1(:);
                    maskTableL2(:,ind) = curMaskL2(:);
                    maskTableL3(:,ind) = curMaskL3(:);
                    
                    edgeTableL0(:,ind) = curEdge(:);
                    edgeTableL1(:,ind) = curEdgeL1(:);
                    edgeTableL2(:,ind) = curEdgeL2(:);
                    edgeTableL3(:,ind) = curEdgeL3(:);
                    
                    edgeIndx(1:numel(indx),ind) = indx(:);
                    edgeIndxL1(1:numel(indxL1),ind) = indxL1(:);
                    edgeIndxL2(1:numel(indxL2),ind) = indxL2(:);
                    edgeIndxL3(1:numel(indxL3),ind) = indxL3(:);
                    
                    edgeIndy(1:numel(indy),ind) = indy(:);
                    edgeIndyL1(1:numel(indyL1),ind) = indyL1(:);
                    edgeIndyL2(1:numel(indyL2),ind) = indyL2(:);
                    edgeIndyL3(1:numel(indyL3),ind) = indyL3(:);
                                        
                    metaTableL0(:,ind) = [Lr Ll length theta x0 y0 x1 y1]';
                    metaTableL1(:,ind) = [Lr Ll length theta+90 -y0 x0 -y1 x1]';
                    metaTableL2(:,ind) = [Lr Ll length theta+180 -x0 -y0 -x1 -y1]';
                    metaTableL3(:,ind) = [Lr Ll length theta+270 y0 -x0 y1 -x1]';
                end
            end    
        end
    end
    
    l0.maskTable = maskTableL0;
    l0.metaTable = metaTableL0;
    l0.edgeTable = edgeTableL0;
    l0.edgeIndx = edgeIndx;
    l0.edgeIndy = edgeIndy;
    
    l1.maskTable = maskTableL1;
    l1.metaTable = metaTableL1;
    l1.edgeTable = edgeTableL1;
    l1.edgeIndx = edgeIndxL1;
    l1.edgeIndy = edgeIndyL1;
    
    l2.maskTable = maskTableL2;
    l2.metaTable = metaTableL2;
    l2.edgeTable = edgeTableL2;
    l2.edgeIndx = edgeIndxL2;
    l2.edgeIndy = edgeIndyL2;
 
    l3.maskTable = maskTableL3;
    l3.metaTable = metaTableL3;
    l3.edgeTable = edgeTableL3;
    l3.edgeIndx = edgeIndxL3;
    l3.edgeIndy = edgeIndyL3;
 
    cd Files;
    save('l0.mat', '-struct', 'l0');
    save('l1.mat', '-struct', 'l1');
    save('l2.mat', '-struct', 'l2');
    save('l3.mat', '-struct', 'l3');
    cd ..;
end

function [L,mask] = getMask(x0,y0,x1,y1,S)
    mask = zeros(S);
    if notInRange(x0,1,S) || notInRange(y0,1,S) || notInRange(x1,1,S) || notInRange(y1,1,S) 
        L =0;
        return;
    end
    
    if x0 == x1 && y0 == y1
        mask(x0,y0) = 1;
        L = 1;
        return;
    end
    
    v = [x1-x0;y1-y0];
    v = v./max(abs(v));

    for i = 1:20 
        x(i) = x0+v(1)*(i-1); 
        y(i) = y0+v(2)*(i-1);
        
        if norm([x(i)-x1,y(i)-y1])<0.1
            x(i) = x1;
            y(i) = y1;
            break;
        end
    end

    xFloor = floor(x);
    yFloor = floor(y);
    xCeil = ceil(x);
    yCeil = ceil(y);
    xWeight = x-xFloor;
    yWeight = y-yFloor;
    
    if isnan(x)
        L =0;
        return;
    end
    
    L = length(x);
    
    for i=2:L-1
        if i==1 || i==L
            val = 1;
        else
            val = 1;
        end
        
        %mask(round(x(i)),round(y(i))) = 1;
        
        mask(xFloor(i),yFloor(i)) = mask(xFloor(i),yFloor(i))+val*(1-xWeight(i))*(1-yWeight(i));
        if yCeil(i)> yFloor(i)
           mask(xFloor(i),yCeil(i)) = mask(xFloor(i),yCeil(i))+val*(1-xWeight(i))*(yWeight(i));
        end
        if xCeil(i)> xFloor(i)
           mask(xCeil(i),yFloor(i)) = mask(xCeil(i),yFloor(i))+val*(xWeight(i))*(1-yWeight(i));
        end
        if yCeil(i)> yFloor(i) && xCeil(i)> xFloor(i) 
           mask(xCeil(i),yCeil(i)) = mask(xCeil(i),yCeil(i))+val*(xWeight(i))*(yWeight(i));
        end
    end
    L = L-1;
    val = 0.5;
    mask(round(x(1)),round(y(1))) = val;
    mask(round(x(end)),round(y(end))) = val;
end

function b = notInRange(v,a,b)
    b = v>=a && v<=b;
    b = ~b;
end

function inTriangle = getInTriangle(edge1,edge2,S,w)
    L = S+2*w;
    inTriangle = zeros(L);
    
    for i = 1:L
        for j = 1:L
            if edge1 == 0 && edge2 == 1 && j>=w+1 && i<=L-w
                inTriangle(i,j) = 1;
            elseif edge1 == 0 && edge2 == 2 && i>=j && j>=w+1
                inTriangle(i,j) = 1;
            elseif edge1 == 1 && edge2 == 2 && i>=j && i<=L-w
                inTriangle(i,j) = 1;
            end
        end
    end
end

