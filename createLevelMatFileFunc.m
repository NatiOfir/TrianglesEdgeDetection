function [ res ] = createLevelMatFileFunc( length,lengthLeft,lenghtRight,lengthHypo,lengthLeftHypo,lenghtRightHypo,numOfEdges )
    curSize = lengthLeft*2;
    res = zeros(curSize,numOfEdges); 

    % first edges pair
    edge1 = 0;
    edge2 = 1;
    sign1 = -1;
    el1 = 0;
    el2 = 2;
    sign2 = 1;
    er1 = 1;
    er2 = 2;
    len1 = length;
    len2 = length;
    subLenLeft = lengthLeft-1;
    subLenRight = lenghtRight-1;

    for v1 = 0:len1-1
        for v2 = 0:len2-1
            vl = v1;
            vr = v2;
            ind = getIndexFromEdge(edge1,v1,edge2,v2,length,lengthHypo);
            arr = zeros(curSize,1);
            for vmid = 0:subLenLeft
                indl = getIndexFromEdge(el1,subLenLeft-vmid,el2,vl,subLenLeft+1,lengthLeftHypo);
                indr = getIndexFromEdge(er1,vmid,er2,vr,subLenRight+1,lenghtRightHypo);
                %res{ind} = [res{ind} sign1*indl sign2*indr];
                arr(vmid*2+1:vmid*2+2,1) = [sign1*indl ;sign2*indr];
            end
            res(:,ind) = arr;
        end
    end

    % second edges pair
    edge1 = 0;
    edge2 = 2;
    sign1 = -1;
    el1 = 0;
    el2 = 2;
    sign2 = -1;
    er1 = 0;
    er2 = 1;
    len1 = length;
    len2 = lenghtRight;
    gap = floor(lengthHypo/2);
    subLenLeft = lengthLeft-1;
    subLenRight = lenghtRight-1;

    for v1 = 0:len1-1
        for v2 = 0:len2-1
            ind = getIndexFromEdge(edge1,v1,edge2,v2,length,lengthHypo);
            arr = zeros(curSize,1);
            vl = v1;
            vr = v2;                   
            for vmid = 0:subLenLeft
                indl = getIndexFromEdge(el1,vmid,el2,vl,subLenLeft+1,lengthLeftHypo);
                indr = getIndexFromEdge(er1,vr,er2,subLenRight-vmid,subLenLeft+1,lenghtRightHypo);
                %res{ind} = [res{ind} sign1*indl sign2*indr];
                arr(vmid*2+1:vmid*2+2,1) = [sign1*indl; sign2*indr];
            end
            res(:,ind) = arr;
        end

        for v2 = len2:lengthHypo-1
            ind = getIndexFromEdge(edge1,v1,edge2,v2,length,lengthHypo);
            indl = getIndexFromEdge(1,v2-gap,2,v1,subLenLeft+1,lengthLeftHypo);
            %res{ind} = [res{ind} sign1*indl nan];
            res(1:2,ind) = [sign1*indl; nan];
        end
    end

    % third edges pair
    edge1 = 1;
    edge2 = 2;
    sign1 = 1;
    el1 = 0;
    el2 = 1;
    sign2 = -1;
    er1 = 1;
    er2 = 2;
    len1 = length;
    len2 = lengthLeft;
    gap = floor(lengthHypo/2);
    subLenLeft = lengthLeft-1;
    subLenRight = lenghtRight-1;

    for v1 = 0:len1-1
        for v2 = gap:len2-1+gap
            ind = getIndexFromEdge(edge1,v1,edge2,v2,length,lengthHypo);
            arr = zeros(curSize,1);
            vr = v1;
            vl = v2-gap;                   
            for vmid = 0:subLenLeft
                indl = getIndexFromEdge(el1,vmid,el2,vl,subLenLeft+1,lengthLeftHypo);
                indr = getIndexFromEdge(er1,subLenRight-vmid,er2,vr,subLenRight+1,lenghtRightHypo);
                %res{ind} = [res{ind} sign1*indl sign2*indr];
                arr(vmid*2+1:vmid*2+2,1) = [sign1*indl ; sign2*indr];
            end
            res(:,ind) = arr;
        end

        for v2 = 0:gap-1
            ind = getIndexFromEdge(edge1,v1,edge2,v2,length,lengthHypo);
            indr = getIndexFromEdge(0,v2,2,v1,subLenRight+1,lenghtRightHypo);
            %res{ind} = [res{ind} nan sign2*indr];
            res(1:2,ind) = [ nan; sign2*indr];
        end
    end

end

function ind = getIndexFromEdge(edgeNumber1,vertexNumber1,edgeNumber2,vertexNumber2,edgeSize,hypoSize)
    subIndEdge = vertexNumber1*edgeSize+vertexNumber2;
    subIndHypo = vertexNumber1*hypoSize+vertexNumber2;
    if edgeNumber1 == 0 && edgeNumber2 == 1
        ind = subIndEdge;
    elseif edgeNumber1 == 0 && edgeNumber2 == 2
        ind = edgeSize^2+subIndHypo;
    elseif edgeNumber1 == 1 && edgeNumber2 == 2
        ind = (edgeSize^2)+edgeSize*hypoSize+subIndHypo;
    end
    ind = ind+1;
end
