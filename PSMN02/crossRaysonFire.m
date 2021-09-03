%% crossRaysonFire
function [crossP,D] = crossRaysonFire(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype)
D = 'nan';

[P1,V1]=findRaysDarcy02(CalibFileCam1,x_pxC1,y_pxC1,Ttype);
[P2,V2]=findRaysDarcy02(CalibFileCam2,x_pxC2,y_pxC2,Ttype);


if size(P1,1) == 0
    crossP = [];
elseif size(P2,1) == 0
    crossP = [];
else
    
    if size(P1,1) == 3
        P1 = P1';
    end
    if size(P2,1) == 3
        P2 = P2';
    end
    
    if isempty(P1)
        %break
    elseif isempty(P2)
        %break
    end
    
    
    clear lineA0 lineA1 lineB0 lineB1
    lineA0 = P1;
    lineA1 = (P1+V1);
    lineB0 = P2;
    lineB1 = (P2+V2);
    [D,Xcp,Ycp,Zcp,Xcq,Ycq,Zcq,Dmin,imin,jmin]= ll_dist3d(lineA0,lineA1,lineB0,lineB1);
    crossP = ([Xcp,Ycp,Zcp]+[Xcq,Ycq,Zcq])/2; % crossing oping
    
end
end

