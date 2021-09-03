%% DARCY02_matchingTracks
function [itraj2,dtraj,listPotTracks,prelist] = DARCY02_matchingTracks(itrajCam0,trajArray_CAM1,trajArray_CAM2RAW,tform1)
% We indicate the trajectory in camera 0,
% it finds the trajectory in camera 1
%
minTintersect = 10; % necessary overlapping time
distThresh = 10;
dstimestep = 5;

tminCAM01 = min(trajArray_CAM1(itrajCam0).track(:,3));
tmaxCAM01 = max(trajArray_CAM1(itrajCam0).track(:,3));
tmean = round(     length(trajArray_CAM1(itrajCam0).track(:,3))/2     );
xcam0 = trajArray_CAM1(itrajCam0).track(tmean,1);
ycam0 = trajArray_CAM1(itrajCam0).track(tmean,2);
[xcam1,ycam1] = transformPointsInverse(tform1,xcam0,ycam0);


% build the list of potential matching trajectory
% criteria for matching:
% 1/ having points at the same time
% 2/ being in the same region
dgi2cam1 = nan(1,length(trajArray_CAM2RAW));
listPotTracks = struct(); % list of potential tracks
ipt = 0; % i possible tracks

% a sort of prematching
for ic = 1 : length(trajArray_CAM2RAW)
    tmean = round( length(trajArray_CAM2RAW(ic).track(:,3))/2 );
    xcam1ic = trajArray_CAM2RAW(ic).track(tmean,1);
    ycam1ic = trajArray_CAM2RAW(ic).track(tmean,2);
    dtraj(ic) = sqrt((xcam1-xcam1ic)^2 + (ycam1-ycam1ic)^2);
end
prelist = find(dtraj<distThresh);

for iic = 1 : length(prelist) %1 : length(trajArray_CAM2RAW)
   ic = prelist(iic);
    tminCAM02 = min(trajArray_CAM2RAW(ic).track(:,3));
    tmaxCAM02 = max(trajArray_CAM2RAW(ic).track(:,3));
    clear A B C
    [A,tcam01,tcam02] = intersect([tminCAM01:tmaxCAM01],[tminCAM02:tmaxCAM02]);

    if length(A) > minTintersect % 1/ having points at the same time
        % 2/ being in the same region
        clear xC2 yC2 dd
        xC2 = [trajArray_CAM2RAW(ic).track(:,1)];
        yC2 = [trajArray_CAM2RAW(ic).track(:,2)];
        dd = sqrt((xC2-xcam1).^2 + (yC2-ycam1).^2);
        dgi2cam1(ic) = min(dd);
        if min(dd) < distThresh
            ipt = ipt + 1;
            listPotTracks(ipt).itr = ic;
            
            clear dsCam0 dsCam1
            dsCam0 = 0; dsCam1 = 0;
            for ids = 1 : dstimestep : length(A)-dstimestep
                xc0i = trajArray_CAM1(itrajCam0).track(tcam01(ids),1);
                yc0i = trajArray_CAM1(itrajCam0).track(tcam01(ids),2);
                xc0f = trajArray_CAM1(itrajCam0).track(tcam01(ids+dstimestep),1);
                yc0f = trajArray_CAM1(itrajCam0).track(tcam01(ids+dstimestep),2);
                xc1i = trajArray_CAM2RAW(ic).track(tcam02(ids),1);
                yc1i = trajArray_CAM2RAW(ic).track(tcam02(ids),2);
                xc1f = trajArray_CAM2RAW(ic).track(tcam02(ids+dstimestep),1);
                yc1f = trajArray_CAM2RAW(ic).track(tcam02(ids+dstimestep),2);
                dsCam0 = dsCam0 + sqrt((xc0i-xc0f)^2+(yc0i-yc0f)^2);
                dsCam1 = dsCam1 + sqrt((xc1i-xc1f)^2+(yc1i-yc1f)^2);
            end
            listPotTracks(ipt).dsCam0 = dsCam0;
            listPotTracks(ipt).dsCam1 = dsCam1;
            
            listPotTracks(ipt).dsRatio = dsCam0 / dsCam1;
            listPotTracks(ipt).distance = min(dd);
        end
    end
end

itraj2 = [];
if ipt > 0
    [~,iitraj2] = min([listPotTracks.distance]);
    if listPotTracks(iitraj2).dsRatio < 1.5 && listPotTracks(iitraj2).dsRatio > 0.5
        itraj2 = listPotTracks(iitraj2).itr;
    else
        itraj2 = [];
    end
end

end

