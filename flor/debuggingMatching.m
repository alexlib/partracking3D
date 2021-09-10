cd('C:\Users\darcy\Desktop\git\Robust-Estimation')
load('all_IFPEN_DARCY02_experiments.mat')

iexpe = 8;
folderScriptshell = allExpeStrct(iexpe).analysisFolder;
folderExperiment  = folderScriptshell;
nameAllTraj = 'alltraj_2021_08_15_electroVavle_at_40percent.mat';
iplane = 31;

allresults = struct();

% load centers
cd(folderExperiment)
allTrajLOAD = load(nameAllTraj);
allTraj = allTrajLOAD.allTraj;

iSeqa = iplane*2-1;
iSeqb = iplane*2;

ImMean = allTraj(iSeqa).ImMean;
him = size(ImMean,1); % 1152;
wim = size(ImMean,2); % 1152;
clear CCtemp CC1 CC2 totalnFrames
CC1 = allTraj(iSeqa).CC; % CCtemp.CC;
CC2 = allTraj(iSeqb).CC; % CCtemp.CC;
totalnFrames = size(CC1,2);
    %  normxcorr2 - we build the images
ACC1 = zeros(him,wim,'uint8');
ACC2 = zeros(him,wim,'uint8');
for it = 1 : totalnFrames
    for ip = 1 : length(CC1(it).X)
        xim1 = round(CC1(it).X(ip));
        yim1 = round(CC1(it).Y(ip));
        ACC1(yim1,xim1) = ACC1(yim1,xim1) + 255;
    end
    for ip = 1 : length(CC2(it).X)
        xim2 = round(CC2(it).X(ip));
        yim2 = round(CC2(it).Y(ip));
        ACC2(yim2,xim2) = ACC1(yim2,xim2) + 255;
    end
end

hcam01 = figure;
sgtitle('Raw images')
subplot(1,2,1)
imagesc(20*ACC1)%, colormap gray
title('CAM1')
subplot(1,2,2)
imagesc(20*ACC2)%, colormap gray
title('CAM2')
%%
filterOrder = 1;
figure;
sgtitle('Filtered images, order=10')
subplot(1,2,1)
imagesc(20*imgaussfilt(ACC1,filterOrder))%, colormap gray
title('CAM1')
subplot(1,2,2)
imagesc(20*imgaussfilt(ACC2,filterOrder))%, colormap gray
title('CAM2')



%%
filterOrder = 10;

% first pass
xm = 00+round(wim/2);
ym = 00+round(him/2);
wsub = 200;%250 %round(0.25*mean(xm,ym)); % width correlation template image

[xoffSet,yoffSet] = imageCorrelation(xm,ym,ACC1,ACC2,wsub,filterOrder);

dxPass01 =   xoffSet-xm;
dyPass01 =   yoffSet-ym;
R = (dxPass01^2+dyPass01^2)^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROBUST ESTIMATION PART 1.3 normxcorr2 pass 02 (on a small window)


wti = 100; % width template images
wstep = 80; % step for sampling the image
nPartMin = 100; % minimum number of particles to calculate the correlation
tmpl_IM_tStr = struct(); % structure storing information on template images
hcam01 = figure('defaultAxesFontSize',20);
imagesc(20*ACC1)%, colormap gray
title('Camera1'), hold on
clear nCol nRow
nCol = wim / wstep;
nLin = him / wstep;
iti = 0;
% cut the image in a lot of small images
clear nCol nRow
nCol = wim / wstep;
nLin = him / wstep;
iti = 0;
for iCol = 1 : nCol
    for iLin = 1 : nLin
        clear xc yc
        xc = round(1 + round(wstep/2) + (iCol-1)*wstep*nCol/floor(nCol));
        yc = round(1 + round(wstep/2) + (iLin-1)*wstep*nLin/floor(nLin));
        if xc-wti/2 < 0 || yc-wti/2 < 0 || xc+wti/2 > wim || yc+wti/2 > him
            continue
        end
        iti = iti + 1;
        tmpl_IM_tStr(iti).x = xc;
        tmpl_IM_tStr(iti).y = yc;
        clear subIm
        subIm =  ACC1(yc-wti/2:yc+wti/2,xc-wti/2:xc+wti/2);
        tmpl_IM_tStr(iti).subIm = subIm;
        tmpl_IM_tStr(iti).meanSubIm = mean(subIm(:));
        if tmpl_IM_tStr(iti).meanSubIm*(101*101)/255 > nPartMin  && ...
                (1.5*dxPass01) + xc + wti/2 > 0 && ...
                (1.5*dyPass01) + yc + wti/2 > 0
            tmpl_IM_tStr(iti).correlable = 1;
            pcol='g';
        else
            tmpl_IM_tStr(iti).correlable = 0;
            pcol='b';
        end
        figure(hcam01)
        clear xp yp
        xp = .5*[-1  1  1 -1 -1]*wti+tmpl_IM_tStr(iti).x;
        yp = .5*[-1 -1  1  1 -1]*wti+tmpl_IM_tStr(iti).y;
        patch('xdata',xp,'ydata',yp,'faceColor','none','faceAlpha',.3,'edgeColor',pcol)
        pause(.2)
        if tmpl_IM_tStr(iti).correlable == 1
            clear xm ym xoffSet yoffSet
            xm = tmpl_IM_tStr(iti).x;
            ym = tmpl_IM_tStr(iti).y;
            [xoffSet,yoffSet] = imageCorrelation(xm,ym,ACC1,ACC2,...
                round(wti/2),filterOrder,'cleanC',dxPass01,dyPass01,150);
            tmpl_IM_tStr(iti).xoffSet = xoffSet;
            tmpl_IM_tStr(iti).yoffSet = yoffSet;

            figure(hcam01)
            if abs(xoffSet-xm- dxPass01)<100 && abs(yoffSet-ym- dyPass01)<100
                quiver(xm,ym,xoffSet-xm,yoffSet-ym,'-r','lineWidth',2)
            else
                tmpl_IM_tStr(iti).correlable = 0;
                quiver(xm,ym,xoffSet-xm,yoffSet-ym,'--r','lineWidth',1)
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build tform1
clear fixedPoints movingPoints
corOK = ([tmpl_IM_tStr.correlable] == 1);
fixedPoints  = [[tmpl_IM_tStr(corOK).x]',      [tmpl_IM_tStr(corOK).y]'];
movingPoints = [[tmpl_IM_tStr(corOK).xoffSet]',[tmpl_IM_tStr(corOK).yoffSet]']; 
transformationType = 'affine'; % 'nonreflectivesimilarity' OR 'similarity'
tform1 = fitgeotrans(movingPoints,fixedPoints,transformationType);

%% check tform1 is OK

ACC2T = transformPointsForward(tform1,ACC2,ACC2);
figure;
imshow(ACC1)
plot(ACC2T(:,1),ACC2T(:,2),'.')
legend('C1','C2on1')


%%

% from BLP TRAJECTOIRE 2D
clear part_cam1 part_cam2 part_cam2RAW
for it = 1 : size(CC1,2)
    part_cam1(it).pos(:,1) = [CC1(it).X]; % out_CAM1(:,1);
    part_cam1(it).pos(:,2) = [CC1(it).Y]; % out_CAM1(:,2);
    part_cam1(it).pos(:,3) = ones(length([CC1(it).X]),1)*it;
    part_cam1(it).intensity = 0; %mI;

    clear cam2X cam2Y
    part_cam2RAW(it).pos(:,1) = [CC2(it).X]; % out_CAM1(:,1);
    part_cam2RAW(it).pos(:,2) = [CC2(it).Y]; % out_CAM1(:,2);
    part_cam2RAW(it).pos(:,3) = ones(length([CC2(it).X]),1)*it;
    part_cam2RAW(it).intensity = 0; %mI;
end

maxdist = 3;
longmin = 8;
[trajArray_CAM1,tracks_CAM1]          = TAN_track2d(part_cam1,maxdist,longmin);
[trajArray_CAM2RAW,tracks_CAM2RAW]    = TAN_track2d(part_cam2RAW,maxdist,longmin);

lcrossStitchTHRSHLD = 4;
itStep = 4;
timeShift = 5;  % 5; % frames
rA = 20;
trajArray_CAM1    = Darcy02stitching(trajArray_CAM1,   lcrossStitchTHRSHLD,itStep,timeShift,rA);
trajArray_CAM2RAW = Darcy02stitching(trajArray_CAM2RAW,lcrossStitchTHRSHLD,itStep,timeShift,rA);

allTraj(iSeqa).trajArray = trajArray_CAM1;
allTraj(iSeqa).tracks    = tracks_CAM1;
allTraj(iSeqb).trajArray = trajArray_CAM2RAW;
allTraj(iSeqb).tracks    = tracks_CAM2RAW;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keep only long trajectories  %  unuseful because trajectories already
% have a minimal length
ltraj = 8; % = 10
clear ikill01 ikill02
for itraj1 = 1 : length(trajArray_CAM1)
    lt1(itraj1) = size(trajArray_CAM1(itraj1).track,1) ;
end
ikill01 = find(lt1<ltraj);

for itraj2 = 1 : length(trajArray_CAM2RAW)
    lt2(itraj2) = size(trajArray_CAM2RAW(itraj2).track,1) ;
end
ikill02 = find(lt2<ltraj);

trajArray_CAM1(ikill01)=[];
trajArray_CAM2RAW(ikill02)=[];


%%%%% %%%%% %%%%% %%%%% %%%%%
%%%%% %%%%% %%%%% %%%%% %%%%%
% STEP 5 - associate trajectories
listMatchedTracks = struct(); % list of potential tracks
ilist = 0;
for itrajCam0 = 1 : length(trajArray_CAM1)
%         if (~mod(itrajCam0,100) == 1) || (itrajCam0==1)
%             fprintf('index trajectory: %0.0f / %0.0f \n',itrajCam0,length(trajArray_CAM1))
%         end
    [itrajCam1,dtraj,listPotTracks,prelist] = DARCY02_matchingTracks(itrajCam0,trajArray_CAM1,trajArray_CAM2RAW,tform1);
    if itrajCam1
        ilist = ilist +1;
        listMatchedTracks(ilist).trajcam0 = itrajCam0;
        listMatchedTracks(ilist).trajcam1 = itrajCam1;
    end
end
%%%%% %%%%% %%%%% %%%%% %%%%%
%%%%% %%%%% %%%%% %%%%% %%%%%
% cross rays with trajectories found with DARCY02_matchingTracks
someTrajectories = struct();
Ttype = 'T1';%T1

%% imageCorrelation
% debug the cleaning of C (varargin in function)
function [xoffSet,yoffSet] = imageCorrelation(xc,yc,ACC1,ACC2,w,filterOrder,varargin)
% varargin: ,'cleanC',dxPass01,dyPass01,R);
ACC1sub = zeros(w+1,w+1,'uint8');
ACC1sub = ACC1(yc-w:yc+w,xc-w:xc+w);
C = normxcorr2(imgaussfilt(ACC1sub,filterOrder),imgaussfilt(ACC2,filterOrder));

% set C to zero above a predefined radius
% Checking varargin structure
if ( length(varargin) > 1 )
    %fprintf('cleaning C \n')
    dxPass01 = double(varargin{:,2});
    dyPass01 = double(varargin{:,3});
    R = double(varargin{:,4});
    x0 = round(xc+dxPass01 + size(ACC1sub,1)/2);
    y0 = round(yc+dyPass01 + size(ACC1sub,2)/2);
    x = 1:size(C,2);
    y = 1:size(C,1);
    [xx,yy] = meshgrid(x,y);
    %     figure
    %     imagesc(C)
    C(((xx-x0).^2+(yy-y0).^2) > R^2)=0;
    %     figure
    %     imagesc(C)
end

[ypeak,xpeak] = find(C==max(C(:)));
yoffSet = ypeak-size(ACC1sub,1) + w;
xoffSet = xpeak-size(ACC1sub,2) + w;
end

%% Darcy02stitching
function trajArray_CAM1_sttchd = Darcy02stitching(trajArray_CAM1,lcrossStitchTHRSHLD,itStep,timeShift,rA)
% lcrossStitchTHRSHLD = 4;
% itStep = 4;
% timeShift = 5;  % 5; % frames
% rA = 20;
%
% trajArray_CAM1_sttchd stitched trajectories

%%% Part 01 simplify tracks
%%%% %%%% %%%%
%%%% %%%% %%%% calculate ds
%%%% %%%% %%%%
for itraj = 1 : size(trajArray_CAM1,2)
    for ittime = 2 : size(trajArray_CAM1(itraj).track,1)
        xi = trajArray_CAM1(itraj).track(ittime-1,1);
        yi = trajArray_CAM1(itraj).track(ittime-1,2);
        xf = trajArray_CAM1(itraj).track(ittime,1);
        yf = trajArray_CAM1(itraj).track(ittime,2);
        trajArray_CAM1(itraj).track(ittime,6) = sqrt((xf-xi)^2+(yf-yi)^2);
    end
    % subset the trajectory with a point every 5 steps in time
    iit = 0;
    % prepare timeSubSample
    iti = 1+itStep/2;
    itf = size(trajArray_CAM1(itraj).track,1)-itStep/2;
    timeSubSample_i = iti :  itStep : round(itf/2);
    timeSubSample_f = itf : -itStep : round(itf/2);
    if timeSubSample_i(end) == timeSubSample_f(end)
        timeSubSample_f(end) = [];
        timeSubSample = [timeSubSample_i,flip(timeSubSample_f)];
    elseif (timeSubSample_f(end)-timeSubSample_i(end)) > 5
        timeSubSample_i = [timeSubSample_i,round((timeSubSample_i(end)+timeSubSample_f(end))/2)];
        timeSubSample   = [timeSubSample_i,flip(timeSubSample_f)];
    else
        timeSubSample = [timeSubSample_i,flip(timeSubSample_f)];
    end
    
    for iittime = 1 : length(timeSubSample) 
        iit = iit + 1;
        ittime = timeSubSample(iittime);
        tmean_i = ittime-itStep/2;
        tmean_f = ittime+itStep/2;
        trajArray_CAM1(itraj).smplTrack(iit,1) = ...
            mean( trajArray_CAM1(itraj).track(tmean_i:tmean_f,1));
        trajArray_CAM1(itraj).smplTrack(iit,2) = ...
            mean( trajArray_CAM1(itraj).track(tmean_i:tmean_f,2)); 
        trajArray_CAM1(itraj).smplTrack(iit,3) = trajArray_CAM1(itraj).track(ittime,3); 
        trajArray_CAM1(itraj).smplTrack(iit,4) = trajArray_CAM1(itraj).track(ittime,4); 
        trajArray_CAM1(itraj).smplTrack(iit,5) = trajArray_CAM1(itraj).track(ittime,5); 
    end
    
    for ittime = 2 : size(trajArray_CAM1(itraj).smplTrack,1)
        xi = trajArray_CAM1(itraj).smplTrack(ittime-1,1);
        yi = trajArray_CAM1(itraj).smplTrack(ittime-1,2);
        xf = trajArray_CAM1(itraj).smplTrack(ittime,1);
        yf = trajArray_CAM1(itraj).smplTrack(ittime,2);
        trajArray_CAM1(itraj).smplTrack(ittime,6) = sqrt((xf-xi)^2+(yf-yi)^2);
    end
    
    trajArray_CAM1(itraj).dsSUM    = sum([trajArray_CAM1(itraj).track(:,6)]);
    trajArray_CAM1(itraj).dsSUMwindow = sum([trajArray_CAM1(itraj).smplTrack(:,6)]);
end
%%%% %%%% %%%%
%%%% %%%% %%%% calculate ds END
%%%% %%%% %%%%



%%% Part 02 stitch
continue2stich = 'on';
conversionsSTR = struct();
icSTR = 0;
while strcmp(continue2stich,'on') % as long as we can stitch we continue to stitch
    itA = 0;
    conversions = 0 ;
    while(1)% trouver une trach qui peut stitcher
        itA = itA+1;
        if itA > size(trajArray_CAM1,2)
        %fprintf('stitched %0.0f trajectories \n',conversions)
            if conversions == 0
                continue2stich = 'off';
            end
            break
        end

    tmaxA = max(trajArray_CAM1(itA).track(:,3));
    itBcandidates = [];
    dABall = [];
    dCrossingCandidates = [];
    for itB = 1 : length(trajArray_CAM1)
        tminB = min(trajArray_CAM1(itB).track(:,3));
        if (tminB - tmaxA) < timeShift && (tminB - tmaxA) > 0
            clear xA yA xB yB dAB
            xA = trajArray_CAM1(itA).track(end,1);
            yA = trajArray_CAM1(itA).track(end,2);
            xB = trajArray_CAM1(itB).track(1,1);
            yB = trajArray_CAM1(itB).track(1,2);
            dAB = sqrt((xA-xB)^2+(yA-yB)^2);
            if dAB < rA
                dABall = [dABall,dAB];
                itBcandidates = [itBcandidates,itB];
                % extrapolate the position of traj A and B and show where the tracer would be
                tA2B = (tmaxA+tminB)/2;
                Dt = (trajArray_CAM1(itA).smplTrack(end,3)-trajArray_CAM1(itA).smplTrack(end-1,3));
                Dx = (trajArray_CAM1(itA).smplTrack(end,1)-trajArray_CAM1(itA).smplTrack(end-1,1));
                Dy = (trajArray_CAM1(itA).smplTrack(end,2)-trajArray_CAM1(itA).smplTrack(end-1,2));
                vAsmplX = Dx / Dt;
                vAsmplY = Dy / Dt;
                xA_extra = trajArray_CAM1(itA).smplTrack(end,1) + ...
                    (tA2B - trajArray_CAM1(itA).smplTrack(end,3)) * vAsmplX;
                yA_extra = trajArray_CAM1(itA).smplTrack(end,2) + ...
                    (tA2B - trajArray_CAM1(itA).smplTrack(end,3)) * vAsmplY;
               
                Dt = (trajArray_CAM1(itB).smplTrack(2,3)-trajArray_CAM1(itB).smplTrack(1,3));
                Dx = (trajArray_CAM1(itB).smplTrack(1,1)-trajArray_CAM1(itB).smplTrack(2,1));
                Dy = (trajArray_CAM1(itB).smplTrack(1,2)-trajArray_CAM1(itB).smplTrack(2,2));
                
                vBsmplX = Dx/Dt;
                vBsmplY = Dy/Dt;
                xB_extra = trajArray_CAM1(itB).smplTrack(1,1) + ...
                    (- tA2B + trajArray_CAM1(itB).smplTrack(1,3)) * vBsmplX;
                yB_extra = trajArray_CAM1(itB).smplTrack(1,2) + ...
                    (- tA2B + trajArray_CAM1(itB).smplTrack(1,3)) * vBsmplY;
                
                dCrossingCandidates = [dCrossingCandidates,sqrt((xB_extra-xA_extra)^2+(yB_extra-yA_extra)^2)];
            end
        end
    end
    
    % stitch the best candidate if it is possible
    [mindist,itBstitch] = min(dCrossingCandidates);
    if mindist < lcrossStitchTHRSHLD        
        itB = itBcandidates(itBstitch);
        
        conversions = conversions + 1;
        
        icSTR = icSTR + 1;
        conversionsSTR(icSTR).mindist = mindist;
        conversionsSTR(icSTR).tmaxA = tmaxA; 
        conversionsSTR(icSTR).tminB = trajArray_CAM1(itB).track(1,3); 
        conversionsSTR(icSTR).itA = itA;
        conversionsSTR(icSTR).itB = trajArray_CAM1(itB).track(1,3); 
        conversionsSTR(icSTR).Ax = trajArray_CAM1(itA).track(:,1);
        conversionsSTR(icSTR).Ay = trajArray_CAM1(itA).track(:,2);
        conversionsSTR(icSTR).At = trajArray_CAM1(itA).track(:,3);
        conversionsSTR(icSTR).Bx = trajArray_CAM1(itB).track(:,1);
        conversionsSTR(icSTR).By = trajArray_CAM1(itB).track(:,2);
        conversionsSTR(icSTR).Bt = trajArray_CAM1(itB).track(:,3);
        xA = trajArray_CAM1(itA).track(end,1);
        yA = trajArray_CAM1(itA).track(end,2);
        xB = trajArray_CAM1(itB).track(1,1);
        yB = trajArray_CAM1(itB).track(1,2);
        dAB = sqrt((xA-xB)^2+(yA-yB)^2);
        conversionsSTR(icSTR).dAB = dAB;
        conversionsSTR(icSTR).xA_extra = xA_extra;
        conversionsSTR(icSTR).yA_extra = yA_extra;
        conversionsSTR(icSTR).xB_extra = xB_extra;
        conversionsSTR(icSTR).yB_extra = yB_extra;
        
        % attach B to A
        xA_f = trajArray_CAM1(itA).track(end,1);
        yA_f = trajArray_CAM1(itA).track(end,2);
        xB_i = trajArray_CAM1(itB).track(1,1);
        yB_i = trajArray_CAM1(itB).track(1,2);
        tA_f = trajArray_CAM1(itA).track(end,3);
        tB_i = trajArray_CAM1(itB).track(1,3);
        for it = tA_f+1 : tB_i-1
            trajArray_CAM1(itA).track(end+1,1) =  xA_f + (xB_i - xA_f) * ((it-tA_f)/(tB_i-tA_f)) ;
            trajArray_CAM1(itA).track(end,2)   =  yA_f + (yB_i - yA_f) * ((it-tA_f)/(tB_i-tA_f)) ;
            trajArray_CAM1(itA).track(end,3)   =  it;
            trajArray_CAM1(itA).track(end,4)   =  0;
        end
        LitA = size(trajArray_CAM1(itA).track  ,1);
        for it = 1 : size(trajArray_CAM1(itB).track  ,1)
            trajArray_CAM1(itA).track(LitA+it,1) =  trajArray_CAM1(itB).track(it,1);
            trajArray_CAM1(itA).track(LitA+it,2) =  trajArray_CAM1(itB).track(it,2);
            trajArray_CAM1(itA).track(LitA+it,3) =  trajArray_CAM1(itB).track(it,3);
            trajArray_CAM1(itA).track(LitA+it,4) =  trajArray_CAM1(itB).track(it,4);
        end
        
        % kill B
        trajArray_CAM1(itB) = [];
        
        % recalculate smplTrack
        trajArray_CAM1(itA).smplTrack  = [];
        itraj = itA;
        iti = 1+itStep/2;
        itf = size(trajArray_CAM1(itraj).track,1)-itStep/2;
        timeSubSample_i = iti :  itStep : round(itf/2);
        timeSubSample_f = itf : -itStep : round(itf/2);
        if timeSubSample_i(end) == timeSubSample_f(end)
            timeSubSample_f(end) = [];
            timeSubSample = [timeSubSample_i,flip(timeSubSample_f)];
        elseif (timeSubSample_f(end)-timeSubSample_i(end)) > 5
            timeSubSample_i = [timeSubSample_i,round((timeSubSample_i(end)+timeSubSample_f(end))/2)];
            timeSubSample = [timeSubSample_i,flip(timeSubSample_f)];
        else
            timeSubSample = [timeSubSample_i,flip(timeSubSample_f)];
        end
        
        iit = 0;
        for iittime = 1 : length(timeSubSample) % 1 : 5 : size(trajArray_CAM1(itraj).track,1)
            iit = iit + 1;
            ittime = timeSubSample(iittime);
            tmean_i = ittime-itStep/2;
            tmean_f = ittime+itStep/2;
            trajArray_CAM1(itraj).smplTrack(iit,1) = ...
                mean( trajArray_CAM1(itraj).track(tmean_i:tmean_f,1));
            trajArray_CAM1(itraj).smplTrack(iit,2) = ...
                mean( trajArray_CAM1(itraj).track(tmean_i:tmean_f,2));
            trajArray_CAM1(itraj).smplTrack(iit,3) = trajArray_CAM1(itraj).track(ittime,3);
            trajArray_CAM1(itraj).smplTrack(iit,4) = trajArray_CAM1(itraj).track(ittime,4);
            trajArray_CAM1(itraj).smplTrack(iit,5) = trajArray_CAM1(itraj).track(ittime,5);
        end
        
        for ittime = 2 : size(trajArray_CAM1(itraj).smplTrack,1)
            xi = trajArray_CAM1(itraj).smplTrack(ittime-1,1);
            yi = trajArray_CAM1(itraj).smplTrack(ittime-1,2);
            xf = trajArray_CAM1(itraj).smplTrack(ittime,1);
            yf = trajArray_CAM1(itraj).smplTrack(ittime,2);
            trajArray_CAM1(itraj).smplTrack(ittime,6) = sqrt((xf-xi)^2+(yf-yi)^2);
        end
        
        trajArray_CAM1(itraj).dsSUM    = sum([trajArray_CAM1(itraj).track(:,6)]);
        trajArray_CAM1(itraj).dsSUMwindow = sum([trajArray_CAM1(itraj).smplTrack(:,6)]);
    end
    
    end
end

trajArray_CAM1_sttchd = trajArray_CAM1;

end
%% matching
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

