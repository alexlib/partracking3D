%% partracking3D

% two problems
% 1/ points de la mires mal croisés en 3D. (crossrays.m ??)
% 2/ tform1 weird with xz
% 3/ is matching of the trajectories OK?

if strcmp(getenv('COMPUTERNAME'),'DESKTOP-3ONLTD9')
    cd('C:\Users\Lenovo\Jottacloud\RECHERCHE\Projets\21_IFPEN\git\partracking3D')
else
    cd('C:\Users\darcy\Desktop\git\Robust-Estimation')
end

load('all_IFPEN_DARCY02_experiments.mat')

% ho to load calib file in the structure
% allExpeStrct(7).calib = load('D:\IFPEN\analysisExperiments\calibFiles\calib_2021_08_13\calib.mat');

iexpe = 8; 

if strcmp(getenv('COMPUTERNAME'),'DESKTOP-3ONLTD9')
    if iexpe == 8
        allExpeStrct(iexpe).inputFolder = 'C:\Users\Lenovo\Jottacloud\RECHERCHE\Projets\21_IFPEN\manips\expe_2021_08_15_40percent';
        allExpeStrct(iexpe).analysisFolder = 'C:\Users\Lenovo\Jottacloud\RECHERCHE\Projets\21_IFPEN\analysis\analysis_2021_08_15';
    end
end
%% find Centers 
allresults = struct();
cd(allExpeStrct(iexpe).analysisFolder)
file_log_ID = fopen(sprintf('log_%s.txt',allExpeStrct(iexpe).name), 'a');
 
planeI = 13;
planeF = 13;

cIN = clock;
allTraj = struct();

for iplane = planeI : planeF
    
    fprintf('plane: %3.0d / %3.0d \n',iplane,planeF)
    iSeqa = iplane*2-1;
    iSeqb = iplane*2;
    
    cPlane_i = clock;
    c1i = clock; fprintf('starts looking for trajectories at %0.2dh%0.2dm\n',c1i(4),c1i(5))
    
    maxdist = allExpeStrct(iexpe).maxdist;
    longmin = allExpeStrct(iexpe).longmin;
    for iSeq = iSeqa:iSeqb % loop on images sequences
        fprintf('working on sequence: %4.0f \n',iSeq)
        clear trajArray_loc tracks_loc CCout
        [CCall,CCout,ImMean,ImMax,filename] = ...,
            DARCY02_findCC(allExpeStrct,iexpe,iSeq,maxdist,longmin,'figures','no');
        %allTraj(iSeq).trajArray = trajArray_loc;
        %allTraj(iSeq).tracks    = tracks_loc;
        allTraj(iSeq).CCall    = CCall;
        allTraj(iSeq).CC       = CCout;
        allTraj(iSeq).ImMean   = ImMean;
        allTraj(iSeq).ImMax    = ImMax;
        allTraj(iSeq).filename = filename;
    end
    
    
end
cOUT = clock;

%%
%% FLOR%% see crossRaysonFire with calibration files
iexpe=8;
calib = allExpeStrct(iexpe).calib;

CalibFileCam1 = calib(:,1);
CalibFileCam2 = calib(:,2);
Ttype= 'T1';

x01 = struct();
y01 = struct();
r3D = struct();
for i = 1:length(CalibFileCam1)   
   % for ipos3D = 1 : length( CalibFileCam1(i).pos3D  )
        % finish the work here
        xy01(i).x = CalibFileCam1(i).pimg(:,1);
        xy01(i).y = CalibFileCam1(i).pimg(:,2);
        xy02(i).x = CalibFileCam2(i).pimg(:,1);
        xy02(i).y = CalibFileCam2(i).pimg(:,2);
   % end
end
colorTime = jet(length(CalibFileCam1));
clear r3D D x3D y3D z3D x_pxC1 y_pxC1 x_pxC2 y_pxC2
figure, hold on, box on
for ipoints = 1:length(CalibFileCam1)
    clear x3D y3D z3D
    for ixy = 1 : length(xy01(ipoints).x)
        x_pxC1 = xy01(ipoints).x(ixy);
        y_pxC1 = xy01(ipoints).y(ixy);
        x_pxC2 = xy02(ipoints).x(ixy);
        y_pxC2 = xy02(ipoints).y(ixy);
        
        [crossP,D] = crossRays(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype);
        if length(crossP)>0
            r3D(ipoints).x(ixy) = crossP(1);
            r3D(ipoints).y(ixy) = crossP(2);
            r3D(ipoints).z(ixy) = crossP(3);
            D(ixy)   = D;
        end
    end
    %plot3(x3D,y3D,z3D,'Color',colorTime(ipoints))
    plot3(r3D(ipoints).x,r3D(ipoints).y,r3D(ipoints).z,'o')
end
xlabel('x')
ylabel('y')
zlabel('z')

%% visualise alltraj
% choose time
iplane = planeI*2; % scan

figure('defaultAxesFontSize',20)
imagesc(allTraj(iplane).ImMax )
colormap gray
clear X Y 
X = allTraj(iplane).CCall(:,1);
Y = allTraj(iplane).CCall(:,2);
hold on
plot(X,Y,'.b')
%set(gca,'ydir','reverse')
hold on

%% save allTraj
cd(allExpeStrct(iexpe).analysisFolder)
%nameMat = 'alltraj_2021_08_15_electroVavle_at_40percent.mat';
nameAllTraj = 'alltraj_2021_05_05.mat';
save(nameAllTraj,'allTraj','-v7.3')
% %%
% planeI = 13;    
% planeF = 13;
%% matching and crossing the rays
planeI=1;
planeF=40;
folderExperiment = allExpeStrct(iexpe).analysisFolder;
folderScriptshell = folderExperiment;
allresults = PSMN_DARCY02(folderScriptshell,folderExperiment,nameAllTraj,iexpe,planeI,planeF);

%% Display results

% all tracks
iplane = planeI*2;


for itrck = 1 : length(allresults(iplane).trajArray_CAM1)
    sizeTRCK1(itrck) = length(allresults(iplane).trajArray_CAM1(itrck).track(:,1));
end
ltrckXY1 = max(sizeTRCK1);
for itrck = 1 : length(allresults(iplane).trajArray_CAM2RAW)
    sizeTRCK2(itrck) = length(allresults(iplane).trajArray_CAM2RAW(itrck).track(:,1));
end
ltrckXY2 = max(sizeTRCK2);

figure('defaultAxesFontSize',20), hold on
clear X1 Y1 X2 Y2
X1 = NaN * ones(ltrckXY1,length(allresults(iplane).trajArray_CAM1));
Y1 = NaN * ones(ltrckXY1,length(allresults(iplane).trajArray_CAM1));
for itrck = 1 : length(allresults(iplane).trajArray_CAM1) 
    X1(1:sizeTRCK1(itrck),itrck) = allresults(iplane).trajArray_CAM1(itrck).track(1:sizeTRCK1(itrck),1);
    Y1(1:sizeTRCK1(itrck),itrck) = allresults(iplane).trajArray_CAM1(itrck).track(1:sizeTRCK1(itrck),2);
end
plot(X1,Y1,'-b','lineWidth',4)
X2 = NaN * ones(ltrckXY2,length(allresults(iplane).trajArray_CAM2RAW));
Y2 = NaN * ones(ltrckXY2,length(allresults(iplane).trajArray_CAM2RAW));
for itrck = 1 : length(allresults(iplane).trajArray_CAM2RAW)
   X2(1:sizeTRCK2(itrck),itrck) = allresults(iplane).trajArray_CAM2RAW(itrck).track(1:sizeTRCK2(itrck),1);
   Y2(1:sizeTRCK2(itrck),itrck) = allresults(iplane).trajArray_CAM2RAW(itrck).track(1:sizeTRCK2(itrck),2);
end
plot(X2,Y2,'-r','lineWidth',4)
%% find track number
xblue = 1976;
yblue =  716;
figure('defaultAxesFontSize',20), hold on, box on
%set(gca,'ydir','reverse')
clear dloc
for itrck = 1 : length(allresults(iplane).trajArray_CAM1) 
   clear X Y 
   X = allresults(iplane).trajArray_CAM1(itrck).track(:,1);
   Y = allresults(iplane).trajArray_CAM1(itrck).track(:,2);
   dloc(itrck) = min(sqrt((X-xblue).^2+(Y-yblue).^2));
end
[~,it1] = min(dloc);
Xb = allresults(iplane).trajArray_CAM1(it1).track(:,1);
Yb = allresults(iplane).trajArray_CAM1(it1).track(:,2);

xblue = 1991;
yblue =  681;
clear dloc
for itrck = 1 : length(allresults(iplane).trajArray_CAM2RAW) 
   clear X Y 
   X = allresults(iplane).trajArray_CAM2RAW(itrck).track(:,1);
   Y = allresults(iplane).trajArray_CAM2RAW(itrck).track(:,2);
   dloc(itrck) = min(sqrt((X-xblue).^2+(Y-yblue).^2));
end
[~,it2] = min(dloc);
clear Xr Yr
Xr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,1);
Yr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,2);

plot(Xb,Yb,'-b','lineWidth',4)
plot(Xr,Yr,'-r','lineWidth',4)

% identify paired tracks and non-paired tracks

%%
calib = allExpeStrct(iexpe).calib;

CalibFileCam1 = calib(:,1);
CalibFileCam2 = calib(:,2);
Ttype= 'T3';

clear x3D y3D z3D
for ixy = 1 : length(x01)
    x_pxC1 = x01(ixy);
    y_pxC1 = y01(ixy);
    x_pxC2 = x02(ixy);
    y_pxC2 = y02(ixy);
    
    [crossP,D] = crossRays(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype);
    if length(crossP)>0
        x3D(ixy) = crossP(1);
        y3D(ixy) = crossP(2);
        z3D(ixy) = crossP(3);
        D(ixy)   = D;
        
        
    end
end
%%
figure('defaultAxesFontSize',20), box on
plot3(x3D,y3D,z3D)


%% matched tracks
figure, hold on, box on, view(3)
xlabel('x')
ylabel('y')
zlabel('z')
iplane = 13;
 for itrck = 1 : length(allresults(iplane).someTrajectories)
    clear X Y
    X = allresults(iplane).someTrajectories(itrck).x3D;
    Y = allresults(iplane).someTrajectories(itrck).y3D;
    Z = allresults(iplane).someTrajectories(itrck).z3D;
    plot3(X,Y,Z,'lineWidth',4)
 end
%% 

%% matched tracks
% test by Flor. Will filter out the tracks that fall 
% outside the deviation

figure, hold on, box on, view(3)
xlabel('x')
ylabel('y')
zlabel('z')
iplane = 33;

 for itrck = 1 : length(allresults(iplane).someTrajectories)
     min = mean(allresults(iplane).hist3D)-std(allresults(iplane).hist3D)
     max = mean(allresults(iplane).hist3D)+std(allresults(iplane).hist3D)
    if mean(allresults(iplane).someTrajectories(itrck).z3D) > min && mean(allresults(iplane).someTrajectories(itrck).z3D) < max
        clear X Y
        X = allresults(iplane).someTrajectories(itrck).x3D;
        Y = allresults(iplane).someTrajectories(itrck).y3D;
        Z = allresults(iplane).someTrajectories(itrck).z3D;
        plot3(X,Y,Z,'lineWidth',4)
    end
 end

%% showing all the planes

%% matched tracks

colorsL = jet(length(allresults));

figure, hold on, box on, view(3)
xlabel('x')
ylabel('y')
zlabel('z')
for  iplane = 1 : 1 : length(allresults)
    clear X Y Z
    for itrck = 1 : length(allresults(iplane).someTrajectories)
        d(itrck) = length(allresults(iplane).someTrajectories(itrck).x3D);
    end
    X = false(max(d),length(allresults(iplane).someTrajectories));
    X = NaN * X;
    Y = false(max(d),length(allresults(iplane).someTrajectories));
    Y = NaN * Y;
    Z = false(max(d),length(allresults(iplane).someTrajectories));
    Z = NaN * Z;
    for itrck = 1 : length(allresults(iplane).someTrajectories)
        X(1:d(itrck),itrck) = allresults(iplane).someTrajectories(itrck).x3D;
        Y(1:d(itrck),itrck) = allresults(iplane).someTrajectories(itrck).y3D;
        Z(1:d(itrck),itrck) = allresults(iplane).someTrajectories(itrck).z3D;
        plot3(X,Y,Z,'lineWidth',4,'color',colorsL(iplane,:))
    end
end

%% DARCY02_findTracks - function part

function [CCall,CCout,ImMean,ImMax,filename] = DARCY02_findCC(allExpeStrct,iexpe,ifile,maxdist,longmin,varargin)

% 1. load image
% 2. subtract mean of the image sequence
% 3. find particles positions on all images

dofigures = 'no';
if numel(varargin)
    dofigures = 'no';
    if strcmp(varargin{2},'yes')
        dofigures = varargin{2};
    end
end


inputFolder = allExpeStrct(iexpe).inputFolder;

cd(inputFolder)
listMcin2 = dir('*.mcin2');
filename  = listMcin2(ifile).name;

cd(inputFolder)
[~,~,params] = mCINREAD2(filename,1,1);
totalnFrames = params.total_nframes;

cd(inputFolder)
[M,~,params]=mCINREAD2(filename,1,totalnFrames);

% calculate mean image
ImMean = uint8(mean(M,3));
ImMax = max(M,[],3);
% subtract
Im01 = M - ImMean;

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% determine particules positions at pixel precision - then refine at
% subpixel precision
th = allExpeStrct(iexpe).centerFinding_th;
sz = allExpeStrct(iexpe).centerFinding_sz;
Nwidth = 1;
CCRAW = struct();

switch dofigures
    case 'yes'
        hmax = figure('defaultAxesFontSize',20);
        imshow(imadjust(ImMax,[1 15]/255))
        colorTime = jet(size(M,3));
end

for it = 1 : size(M,3)
    Im = zeros(size(Im01,1),size(Im01,2),class(Im01));
    Im(:,:) = Im01(:,:,it);
    
    CCRAW(it).xyRAW = pkfnd(Im01(:,:,it),th,sz);
    
    %     figure
    %     imagesc(Im01(:,:,it))
    %     hold on
    %     CCooo = CCRAW(it).xyRAW;
    %     plot(CCooo(:,1),CCooo(:,2),'or')
    %     pause(.5)
    %     close all
    
    
    %refine at subpixel precision
    for ixy = 1 : size(CCRAW(it).xyRAW,1)
        clear xpkfnd ypkfnd Ip
        Ip = zeros(2*Nwidth+1,2*Nwidth+1,'double');
        
        xpkfnd = CCRAW(it).xyRAW(ixy,1);
        ypkfnd = CCRAW(it).xyRAW(ixy,2);
        Ip = double(Im(ypkfnd-Nwidth:ypkfnd+Nwidth,xpkfnd-Nwidth:xpkfnd+Nwidth));
        % replace all 0 Ip values by a tiny value
        list0 = find(Ip==0);
        Ip(list0) = 1e-6;
        
        CC(it).xy(ixy,1) = xpkfnd + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
        CC(it).xy(ixy,2) = ypkfnd + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
    end
    switch dofigures
        case 'yes'
    if size(CCRAW(it).xyRAW,1) >0
        CCwrk = CC(it).xy;
        figure(hmax), title(sprintf('time: %4.0f - nparts: %4.0f',it,length(CCwrk)))
        hold on
        hp = plot(CCwrk(:,1),CCwrk(:,2),'o','markerFaceColor',colorTime(it,:),...
            'markerEdgeColor','none','markerSize',5);
    end
    end
end
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% remove the NaNs for all t
if exist('CC')
    for it = 1 : size(CC,2)
        
        if size(CCRAW(it).xyRAW,1) >0
            clear ikill CCX CCY
            ikill = [];
            for ip = 1 : size(CC(it).xy,1)
                CCX(ip) = CC(it).xy(ip,1);
                CCY(ip) = CC(it).xy(ip,2);
                if isnan(CC(it).xy(ip,1)) || isnan(CC(it).xy(ip,2))
                    ikill = [ikill,ip];
                end
            end
            CCX(ikill) = [];
            CCY(ikill) = [];
            CC(it).xy = [];
            for ip = 1 : length(CCX)
                CC(it).xy(ip,1) = CCX(ip);
                CC(it).xy(ip,2) = CCY(ip);
            end
            
            % fprintf('ifile: %4.0f - it: %4.0f - n part found: %4.0f n_found_RAW: %4.0f \n',...
            %     ifile,it,length(CCX),size(CC(it).xy,1))
        end
        
        
    end
end
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% put all positions in only one variable called CCall
% also in CCout, an output of the function for later calculations
clear CCall CCout
if exist('CC')
for it = 1 : size(M,3)
    %fprintf('ifile is %4.0f and it is: %4.0f\n',ifile,it)
    if  ((length(CC) >= it)) && size(CC(it).xy,1) > 0
        X = CC(it).xy(:,1);
        Y = CC(it).xy(:,2);
        CCout(it).X = X;
        CCout(it).Y = Y;
        T = it * ones(1,length(X));
        %if it == 1
        if ~exist('CCall')
            CCall = [X,Y];
            CCall(:,3) = [T];
        else
            CCtemp = [X,Y];
            CCtemp(:,3) = [T];
            CCall = [CCall;CCtemp];
        end
    end
end
else
    CCall = NaN;
    CCout = NaN;
end
end

%%

%% PSMN_DARCY02
function allresults = PSMN_DARCY02(folderScriptshell,folderExperiment,nameAllTraj,iexpe,planeI,planeF)
%
% folderScriptshell
% folderExperiment = 'D:\IFPEN\analysisExperiments\analysis_expe_2021_06_17\run02';
% nameAllTraj
% iexpe = 6 / 1 / 2 / 3
% planeI = 1;
% planeF = 1;
%
%
%

fprintf('start \n' );

cd(folderScriptshell)
load('all_IFPEN_DARCY02_experiments.mat')

allresults = struct();
% CHECK WITH FLOR : cd(allExpeStrct(iexpe).analysisFolder)
% cd(folderExperiment)
% file_log_ID = fopen(sprintf('log_%s.txt',allExpeStrct(iexpe).name), 'a');

% load centers
cd(folderExperiment)
allTrajLOAD = load(nameAllTraj);
allTraj = allTrajLOAD.allTraj;

for iplane = planeI : planeF
    
    fprintf('plane %0.3d - start \n', iplane  );
    
    cPlane_i = clock;
    
    iSeqa = iplane*2-1;
    iSeqb = iplane*2;

    ImMean = allTraj(iSeqa).ImMean;
    him = size(ImMean,1); % 1152;
    wim = size(ImMean,2); % 1152;
    clear CCtemp CC1 CC2 totalnFrames
    CC1 = allTraj(iSeqa).CC; % CCtemp.CC;
    CC2 = allTraj(iSeqb).CC; % CCtemp.CC;
    totalnFrames = size(CC1,2);
    
    
fprintf('plane %0.3d - normxcorr2 \n', iplane  );
ci = clock;
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
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normxcorr2 pass 01 (on a large window)
    % xm,ym : fixed points in camera 1
    filterOrder = 10;
    
    % first pass
    xm = 00+round(wim/2);
    ym = 00+round(him/2);
    wsub = 250; %round(0.25*mean(xm,ym)); % width correlation template image
    
    [xoffSet,yoffSet] = imageCorrelation(xm,ym,ACC1,ACC2,wsub,filterOrder);
    
    dxPass01 =   xoffSet-xm;
    dyPass01 =   yoffSet-ym;
    R = (dxPass01^2+dyPass01^2)^(1/2);

    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ROBUST ESTIMATION PART 1.3 normxcorr2 pass 02 (on a small window)
    
    fprintf('plane %0.3d - find corresponding points between the 2 cameras for tform1 \n', iplane  );
    ci = clock;
    
    wti = 200; % width template images
    wstep = 100; % step for sampling the image
    nPartMin = 100; % minimum number of particles to calculate the correlation
    tmpl_IM_tStr = struct(); % structure storing information on template images
    
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
            else
                tmpl_IM_tStr(iti).correlable = 0;
            end
            
            if tmpl_IM_tStr(iti).correlable == 1
                clear xm ym xoffSet yoffSet
                xm = tmpl_IM_tStr(iti).x;
                ym = tmpl_IM_tStr(iti).y;
                [xoffSet,yoffSet] = imageCorrelation(xm,ym,ACC1,ACC2,...
                    round(wti/2),filterOrder,'cleanC',dxPass01,dyPass01,150);
                tmpl_IM_tStr(iti).xoffSet = xoffSet;
                tmpl_IM_tStr(iti).yoffSet = yoffSet;
            end
        end
    end
        cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    fprintf('plane %0.3d - tform1 \n', iplane  );
    ci = clock;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build tform1
    clear fixedPoints movingPoints
    corOK = ([tmpl_IM_tStr.correlable] == 1);
    fixedPoints  = [[tmpl_IM_tStr(corOK).x]',      [tmpl_IM_tStr(corOK).y]'];
    movingPoints = [[tmpl_IM_tStr(corOK).xoffSet]',[tmpl_IM_tStr(corOK).yoffSet]']; 
    transformationType = 'nonreflectivesimilarity';
    tform1 = fitgeotrans(movingPoints,fixedPoints,transformationType);
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    fprintf('plane %0.3d - find trajectories \n', iplane  );
    ci = clock;
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
    
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    fprintf('plane %0.3d - associate trajectories \n', iplane  );
    ci = clock;
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
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    fprintf('plane %0.3d - cross rays \n', iplane  );
    ci = clock;
    %%%%% %%%%% %%%%% %%%%% %%%%%
    %%%%% %%%%% %%%%% %%%%% %%%%%
    % cross rays with trajectories found with DARCY02_matchingTracks
    someTrajectories = struct();
    Ttype = 'T1';
    
    % CHECK WITH FLOR:
    %     CalibFile = allExpeStrct(iexpe).CalibFile;
    %     %calibTemp = load(CalibFile,'calib'); calib = calibTemp.calib;
    %     calib = allExpeStrct(iexpe).calib;
    calib = allExpeStrct(iexpe).calib;
    
    CalibFileCam1 = calib(:,1);
    CalibFileCam2 = calib(:,2);
    
    for iselTraj = 1 : size(listMatchedTracks,2)
        itraj1 = listMatchedTracks(iselTraj).trajcam0;
        itraj2 = listMatchedTracks(iselTraj).trajcam1;

        someTrajectories(iselTraj).itraj1 = itraj1;
        someTrajectories(iselTraj).itraj2 = itraj2;
        
        % cross the two choosen rays
        clear x01 y01 x02 y02 x02incam01 y02incam01
        tminCAM01 = min(trajArray_CAM1(itraj1).track(:,3));
        tmaxCAM01 = max(trajArray_CAM1(itraj1).track(:,3));
        tminCAM02 = min(trajArray_CAM2RAW(itraj2).track(:,3));
        tmaxCAM02 = max(trajArray_CAM2RAW(itraj2).track(:,3));
        [A,B,C] = intersect([tminCAM01:tmaxCAM01],[tminCAM02:tmaxCAM02]);

        if A
            clear x01 y01 x02 y02
            x01 = trajArray_CAM1(itraj1).track(min(B):max(B),1);
            y01 = trajArray_CAM1(itraj1).track(min(B):max(B),2);
            x02 = trajArray_CAM2RAW(itraj2).track(min(C):max(C),1);
            y02 = trajArray_CAM2RAW(itraj2).track(min(C):max(C),2);
            
            clear x_pxC1 y_pxC1 x_pxC2 y_pxC2
            for ixy = 1 : length(x01)
                
                x_pxC1 = x01(ixy);
                y_pxC1 = y01(ixy);
                x_pxC2 = x02(ixy);
                y_pxC2 = y02(ixy);
                
                [crossP,D] = crossRays(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype);
                if length(crossP)>0
                    someTrajectories(iselTraj).x3D(ixy) = crossP(1);
                    someTrajectories(iselTraj).y3D(ixy) = crossP(2);
                    someTrajectories(iselTraj).z3D(ixy) = crossP(3);
                    %someTrajectories(iselTraj).t(ixy)   = t;         % image number
                    %someTrajectories(iselTraj).t(ixy)   = 20210617T161636-20210617T154756;         % absolute time
                    %someTrajectories(iselTraj).t(ixy)   = t/fps;         % relative time
                    someTrajectories(iselTraj).D(ixy)   = D;
                end
            end
        end
        
        
    end
    
    %%%%%%%%%
    hist3D = [];
    histD  = [];
    % plot the result
    for itraj3D = 1 : size(someTrajectories,2)
        hist3D = [hist3D,[someTrajectories(itraj3D).z3D]]; % histogram z of particles
        histD  = [histD, [someTrajectories(itraj3D).D]];   % histogram distance between the rays
    end
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    
    fprintf('plane %0.3d - save data \n', iplane  );
    ci = clock;
    
    allresults(iplane).someTrajectories = someTrajectories;
    allresults(iplane).hist3D = hist3D;
    allresults(iplane).histD = histD;
    allresults(iplane).trajArray_CAM1 = trajArray_CAM1;
    allresults(iplane).trajArray_CAM2RAW = trajArray_CAM2RAW;
    allresults(iplane).listMatchedTracks = listMatchedTracks;
    allresults(iplane).tform1 = tform1;
    
    cd(folderExperiment)
    save('allResults_auFilDeLEau.mat','allresults')
    
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    cPlane_f = clock;
    % fprintf(file_log_ID, 'just run plane %0.3d in %0.0f s \n', iplane , etime(cPlane_f,cPlane_i) );
    fprintf('just run plane %0.3d in %0.0f s \n', iplane , etime(cPlane_f,cPlane_i) );
    
end

fprintf('end \n')
% fclose(file_log_ID)

end

%% FUNCTIONS

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

%% TAN_track2d
function [traj,tracks]=TAN_track2d(pos,maxdist,longmin)

% nearest neighbor particle tracking algo
% pos is a structure with field .pos : pos.pos(:,1)=x1 ; pos.pos(:,2)=x2
% pos.pos(:,3)=framenumber
% frame number must an integer from 1->N
% maxdist=1; disc radius in where we are supposed to find a particle in
% next frame
% longmin : minimum length trajectory
%
% in case there are few frames bit many particles this code could be
% improved using say 100^2 vertex (only try particle in nearest vertices)
%
% example:
% load tracking_rom_positions.mat
% tracks=track2d_rom(positions,1,50);
%
% use plot_traj to display results
%
% written by r. volk 09/2014 (modified 01/2020)

tic;

tracks = pos;
for ii = 1:size(tracks,2)
    % frame number in the trajectory (from 1 ->p) for a trajectory of length p
    tracks(ii).pos(:,4) = zeros(size(tracks(ii).pos,1),1);
    % 5th column 1 if particle is free to be linked to a new trajectory, 0 if no, 2 if
    % linked to 2 or more trajectories.
    tracks(ii).pos(:,5) = ones(size(tracks(ii).pos,1),1);
    % we don't create trajectories, we only write numbers in the 4th and 5th
    % columns
end


% number of active tracks, we strat with frame 1
ind_actif = (1:size(tracks(1).pos,1));

tracks(1).pos(:,4) = (1:size(tracks(1).pos,1));
tracks(1).pos(:,5) = zeros(size(tracks(1).pos,1),1);

% number of trajectories created at this step
% will increase each time we create a new trajectory
numtraj=size(tracks(1).pos,1);

% loop over frames
%tic
for kk=2:size(tracks,2)
    % frame number we are looking at
    %numframe=kk;
    % indices of those paricles
    %ind_new=find(tracks(:,3)==numframe);

    % loop over active particles in previous frame (kk-1)
    for ll=1:length(ind_actif)
        % position of particle ll in frame kk-1
        actx = tracks(kk-1).pos(ind_actif(ll),1);
        acty = tracks(kk-1).pos(ind_actif(ll),2);

        % trajectory number of the active particle (frame kk-1)
        actnum = tracks(kk-1).pos(ind_actif(ll),4);

        % could add a tag: frame number in this trajectory
        % si tag<=2 : rien
        % si tag==3 actx=actx+vx*dt avec dt=1 ici (3 frames best estimate)
        % si tag==4 actx

        % new particle positions in frame kk
        newx = tracks(kk).pos(:,1);
        newy = tracks(kk).pos(:,2);


        % compute distance
        dist = sqrt((actx-newx).^2+(acty-newy).^2);
        % take the min
        [dmin,ind_min]=min(dist);

        % test with maxdist criterion
        if dmin < maxdist
            dispo = tracks(kk).pos(ind_min,5);

            if dispo==1
                % part is dispo=free we change dispo into 0
                tracks(kk).pos(ind_min,5) = 0;

                % we link the particle to the active particle set
                % trajectory number equal to the one of the active
                % particle
                tracks(kk).pos(ind_min,4) = actnum;

            elseif dispo==0
                % the part already linked, change dispo into 2
                % can't be linked to 2 trajectories
                tracks(kk).pos(ind_min,5) = 2;

                % and we set its trajectory number to zero
                % will be rejected at the end
                tracks(kk).pos(ind_min,4) = 0;

            end
        end
    end

        % define particles to be tracked
        % keep particles found only one time, and the non found particles
        % those will create new trajectories
        ind_actif = find(tracks(kk).pos(:,5)==0);

        % new (not found) particles are given a new trajectory number
        ind_new_part = find(tracks(kk).pos(:,5)==1);

        % if there are new particles
        if isempty(ind_new_part)==0
            % loop of new part -> increase numtraj
            for mm=1:length(ind_new_part)
                numtraj = numtraj + 1;
                tracks(kk).pos(ind_new_part(mm),4) = numtraj;
            end
        end
        ind_actif=[ind_actif;ind_new_part];
    %toc
end

% write trajectories in right order.
% reject all particles found 2 times
% keep only trajectories longer than longmin

tracks_array = cat(1,tracks.pos); %put all traj in an array
tracks_array = sortrows(tracks_array,4); %sort traj in ascendent order

tracks_array = tracks_array(tracks_array(:,4)~=0,:); %kick numtraj == 0 
fstart = zeros(length(tracks_array),1);
fstart(1) = 1; fstart(2:end) = diff(tracks_array(:,4)); %make diff -> 0 when two 2 succesive lines belong to the same traj, 1 if not
fstart = find(fstart==1); %find indices corresponding to the start of trajectories

flong = diff(fstart); %find the length of trajectories (except the last one)
last_traj = find(tracks_array(:,4)==length(fstart)); %find the length of the last trajectory
flong(end+1) = length(last_traj);

ftrackselec = find(flong>=longmin); %select traj longer than longmin

traj = struct(); %rearange in a structure
for kk = 1 :length(ftrackselec)
    ll  = flong(ftrackselec(kk));
    deb = fstart(ftrackselec(kk));
    traj(kk).track = tracks_array(deb:deb+ll-1,:);
end
end

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

%% crossRays
function [crossP,D] = crossRays(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype)
D = 'nan';

 [P1,V1]=findRaysDarcy02(CalibFileCam1,x_pxC1,y_pxC1,Ttype);
 [P2,V2]=findRaysDarcy02(CalibFileCam2,x_pxC2,y_pxC2,Ttype);

%[P1,V1]=findRaysDarcy02_smallTarget(CalibFileCam1,x_pxC1,y_pxC1,Ttype);
%[P2,V2]=findRaysDarcy02_smallTarget(CalibFileCam2,x_pxC2,y_pxC2,Ttype);

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

%% findRaysDarcy02
function [P,V,XYZ]=findRaysDarcy02(calib,x_px,y_px,Ttype)
%%% calib : calibration data for this camera
%%% x_px  : x coordinates in px,
%%% y_px  : y coordinates in px,
%%% Ttype : type of the transformation to use (T1=Linear, T3=Cubic).

% calibTemp = load(CalibFile,'calib'); calib = calibTemp.calib;

Npart = numel(x_px);
Nplans = numel(calib);

XYZ = zeros(numel(calib),3,numel(x_px));

% for kplan = 1:Nplans
%     I = inpolygon(x_px,y_px,calib(kplan).pimg(calib(kplan).cHull,1),calib(kplan).pimg(calib(kplan).cHull,2));
%     if max(I)>0
%         if Ttype=='T1'
%             [Xtmp,Ytmp]=transformPointsInverse((calib(kplan).T1px2rw),x_px(I==1),y_px(I==1));
%         elseif Ttype=='T3'
%             [Xtmp,Ytmp]=transformPointsInverse((calib(kplan).T3px2rw),x_px(I==1),y_px(I==1));
%         end
%         
%         XYZ(kplan,1,I==1)=Xtmp;
%         XYZ(kplan,2,I==1)=Ytmp;
%         XYZ(kplan,3,I==1)=calib(kplan).posPlane;
%     end
%     
%     XYZ(kplan,1,I==0) = NaN;
%     XYZ(kplan,2,I==0) = NaN;
%     XYZ(kplan,3,I==0) = NaN;
% end
for kplan = 1:Nplans
    I = inpolygon(x_px,y_px,calib(kplan).pimg(calib(kplan).cHull,1),calib(kplan).pimg(calib(kplan).cHull,2));
    if max(I)>0
        if Ttype=='T1'
            [Xtmp,Ytmp]=transformPointsForward((calib(kplan).T1px2rw),x_px(I==1),y_px(I==1));
        elseif Ttype=='T3'
            [Xtmp,Ytmp]=transformPointsForward((calib(kplan).T3px2rw),x_px(I==1),y_px(I==1));
        end
        
        XYZ(kplan,1,I==1)=Xtmp;
        XYZ(kplan,2,I==1)=Ytmp;
        XYZ(kplan,3,I==1)=calib(kplan).posPlane;
    end
    
    XYZ(kplan,1,I==0) = NaN;
    XYZ(kplan,2,I==0) = NaN;
    XYZ(kplan,3,I==0) = NaN;
end



[P, V] = fit3Dline(XYZ);


end

%% findRaysDarcy02_smallTarget

function [P,V,XYZ]=findRaysDarcy02_smallTarget(calib,x_px,y_px,Ttype)
% the line  I = inpolygon(x_px,y_px,calib(kplan).pimg(calib(kplan).cHull,1),calib(kplan).pimg(calib(kplan).cHull,2));
% is removed for situations when the target really doesn't fill the field
% of view of the cameras

Npart = numel(x_px);
Nplans = numel(calib);

XYZ = zeros(numel(calib),3,numel(x_px));

for kplan = 1:Nplans
    if Ttype=='T1'
        [Xtmp,Ytmp]=transformPointsInverse((calib(kplan).T1px2rw),x_px,y_px);
    elseif Ttype=='T3'
        [Xtmp,Ytmp]=transformPointsInverse((calib(kplan).T3px2rw),x_px,y_px);
    end
    XYZ(kplan,1)=Xtmp;
    XYZ(kplan,2)=Ytmp;
    XYZ(kplan,3)=calib(kplan).posPlane;
end
[P, V] = fit3Dline(XYZ);
end



%% fit3Dline
function [xyz0,direction] = fit3Dline(XYZ)

if max(max(max(isnan(XYZ)))) ==0
    [xyz0,direction] = fit3Dline_nonan(XYZ);
else
    [P V] = arrayfun(@(I)(fit3Dline_nan(XYZ(:,:,I))),1:size(XYZ,3),'UniformOutput',false);
    xyz0 = (cell2mat(P'));
    direction = (cell2mat(V'));
    
    xyz0(isnan(xyz0)) = [];
    direction(isnan(direction)) = [];
end

end

%% fit3Dline_nan
function [xyz0,direction]=fit3Dline_nan(XYZ)
%%% [xyz0,direction]=fit3Dline_jv(XYZ)
%
% @MBourgoin 01/2019

I = find(isnan(XYZ(:,1)));
XYZ(I,:)=[];

if size(XYZ,1)>2
    
    xyz0=mean(XYZ);
    %xyz0=cell2mat(arrayfun(@(x) mean(x.CCrw),Proj,'UniformOutput',false));
    
    A=bsxfun(@minus,XYZ,xyz0); %center the data
    
    % xyz0=XYZ(3,:);
    % A= XYZ;
    
    % xyz0=XYZ(plan_centre,:);
    % A=bsxfun(@minus,XYZ,xyz0); %center the data
    
    %[U,S,V]=svd(A);
    [Uac Sac Vac]=arrayfun(@(kkk) svd(A(:,:,kkk)),[1:size(A,3)],'UniformOutput',false);
    Ua=cat(3,Uac{:});
    Sa=cat(3,Sac{:});
    Va=cat(3,Vac{:}); clear Uac Sac Vac;
    
    %direction=cross(V(:,end),V(:,end-1));
    dd=arrayfun(@(x) cross(Va(:,end,x),Va(:,end-1,x)),[1:size(Va,3)],'UniformOutput',false);
    direction=cat(3,dd{:})';  clear dd;
else
    %xyz0 = [NaN NaN NaN];
    %direction = [NaN NaN NaN];
    xyz0=[];
    direction=[];
end

%line = [xyz0'  direction];
end

%% fit3Dline_nonan
function [xyz0,direction]=fit3Dline_nonan(XYZ)

% @JVessaire 01/2019

xyz0=mean(XYZ,1);
Aa=bsxfun(@minus,XYZ,xyz0); %center the data
xyz0=squeeze(xyz0)';

%Aa=permute(A,[3 2 1]);

[~, ~, Vac]=arrayfun(@(kkk) svd(Aa(:,:,kkk)),[1:size(Aa,3)],'UniformOutput',false);
Va=cat(3,Vac{:});

dd=arrayfun(@(x) cross(Va(:,end,x),Va(:,end-1,x)),[1:size(Va,3)],'UniformOutput',false);
direction=cat(2,dd{:})'; clear dd Vac A;

end

%% ll_dist3d
function [D,Xcp,Ycp,Zcp,Xcq,Ycq,Zcq,Dmin,imin,jmin]= ll_dist3d(P0,P1,Q0,Q1)
%ll_dist3d - Find the distances between each pair of straight 3D lines in
% two sets. Find the closest points on each pair, and the pair with minimum
% distance. Each line is defined by two distinct points.
%
% Input:
% P0 - array of first points of the first set (m X 3), where m is the
% number of lines in the first set. P0(j,1), P0(j,2), P0(j,3) are X, Y
% and X coordinates, accordingly, of point j.
% Pl - array of second points of the first set (m X 3), where m is the
% number of lines in the first set. P1(j,1), Pl(j,2), Pl(j,3) are X, Y
% and X coordinates, accordingly, of point j.
% Q0 - array of first points of the second set (n % 3), where n is the
% number of lines in the second set. Q0(k,1), Q0(k,2), Q0(k,3) are X, Y
% and X coordinates, accordingly, of point k.
% Ql - array of second points of the second set (n % 3), where n is the
% number of lines in the second set. Q0(k,1), Q0(k,2), Q0(k,3) are X, Y
% and X coordinates accordingly of point k.
% Output:
% D - array of distances between line pairs (m X n). D(j,k) is the
% distance between line j from the first (P) set, and line k from the
% second (Q) set.
% Xcp - array of X coordinates of closest points belonging to the first
% (P) set (m X n). Xcp(j,k) is an % coordinate of the closest point on a
% line j defined by P0(j,:) and P1(j,:), computed to the line k defined
% by Q0(k,:) and Q1(k,:).
% Ycp - array of Y coordinates of closest points belonging to the first
% (P) set (m X n). See Xcp definition.
% Zcp - array of Y coordinates of closest points belonging to the first
% (P) set (m X n). See Xcp definition.
% Xcq - array of X coordinates of closest points belonging to the second
% (Q) set (m X n). Xcq(j,k) is an % coordinate of the closest point on a
% line k defined by Q0(k,:) and Q1(k,:), computed to the line j defined
% by P0(j,:) and P1(1,:).
% Ycq - array of Y coordinates of closest points belonging to the second
% (Q) set (m X n). See Xcq definition.
% Zcq - array of % coordinates of closest points belonging to the second
% (Q) set (m X n). See Xcq definition.
%
% Remarks:
% Below is a simple unit test for this function. The test creates
% 2 sets of random 3D lines, finds the distances between each pair of
% lines, and plots the pair with shortest distance
% To run the test, uncommnent the following lines:
%
% n1 = 4; % number of lines in first set
% n2 = 2; % number of lines in first set
% P0 = rand(n1,3); P1 = rand(n1,3); Q0 = rand(n2,3); Q1 = rand(n2,3);
% [D,Xcp,Ycp,Zcp,Xcq,Ycq,Zcq,Dmin,imin,jmin] = ll_dist3d(P0, P1, Q0, Q1);
% t = (-2:0.01:2);
% Tp = repmat(t(:), 1, size(P0,1));
% Tq = repmat(t(:), 1, size(Q0,1));
% Xp = repmat(P0(:,1)',[size(t,2), 1]) + Tp.*(repmat(P1(:,1)',[size(t,2),1])-...
% repmat(P0(:,1)', size(t,2), 1));
% Yp = repmat(P0(:,2)',[size(t,2), 1]) + Tp.*(repmat(P1(:,2)',[size(t,2),1])-...
% repmat(P0(:,2)', size(t,2), 1));
% Zp = repmat(P0(:,3)',[size(t,2), 1]) + Tp.*(repmat(P1(:,3)',[size(t,2),1])-...
% repmat(P0(:,3)', size(t,2), 1));
% Xq = repmat(Q0(:,1)', size(t,2), 1) + Tq.*(repmat(Q1(:,1)',size(t,2),1)-...
% repmat(Q0(:,1)', size(t,2), 1));
% Yq = repmat(Q0(:,2)',size(t,2), 1) + Tq.*(repmat(Q1(:,2)',size(t,2),1)-...
% repmat(Q0(:,2)', size(t,2), 1));
% Zq = repmat(Q0(:,3)',size(t,2), 1) + Tq.*(repmat(Q1(:,3)',size(t,2),1)-...
% repmat(Q0(:,3)', size(t,2), 1));
% figure;
% plot3(Xp(:,imin),Yp(:,imin),Zp(:,imin),Xq(:,jmin),Yq(:,jmin),Zq(:,jmin));
% hold on
% plot3(Xcp(imin,jmin),Ycp(imin,jmin),Zcp(imin,jmin),'ro',Xcq(imin,jmin),Ycq(imin,jmin),Zcq(imin,jmin),'mo');
% axis equal
% grid on
% xlabel('X'); ylabel('Y'); zlabel('Z');
%
% Revision history:
% March 03, 2016 - created (Michael Yoshpe)
%**************************************************************************
% check inputs validity
[mp0, np0] = size(P0);
if(np0 ~=3 )
    error('Array P0 should of size (m X 3)');
end
[mpl, npl] = size(P1);
if((mpl ~= mp0) || (npl ~= np0))
    error('P0 and Pl arrays must be of same size');
end
[mq0, nq0] = size(Q0);
if(nq0 ~= 3)
    error('Array Q0 should of size (n X 3)');
end
[mq1, nq1] = size(Q1);
if((mq1 ~= mq0) || (nq1 ~= nq0))
    error('Q0 and Ql arrays must be of same size');
end
u = P1 - P0; % vectors from P0 to P1
uu = repmat(u,[1,1,mq0]);
v = Q1 - Q0; % vectors from Q0 to Q1
vv = permute(repmat(v,[1,1,mp0]), [3 2 1]);
PP0 = repmat(P0,[1,1,mq0]);
QQ0 = permute(repmat(Q0,[1,1,mp0]), [3 2 1]);
w0 = PP0 - QQ0;
aa = dot(uu,uu,2);
bb = dot(uu,vv,2);
cc = dot(vv,vv,2);
dd = dot(uu,w0,2);
ee = dot(vv,w0,2);
ff = aa.*cc - bb.*bb;
idx_par = (ff < 5*eps); % indices of parallel lines
idx_nonpar = ~idx_par; % indices of non-parallel lines
sc = NaN(mp0,1,mq0);
tc = NaN(mp0,1,mq0);
sc(idx_nonpar) = (bb(idx_nonpar).*ee(idx_nonpar) - ...
    cc(idx_nonpar).*dd(idx_nonpar))./ff(idx_nonpar);
tc(idx_nonpar) = (aa(idx_nonpar).*ee(idx_nonpar) - ...
    bb(idx_nonpar).*dd(idx_nonpar))./ff(idx_nonpar);
PPc = PP0 + repmat(sc, [1,3,1]).*uu;
QQc = QQ0 + repmat(tc, [1,3,1]).*vv;
Xcp = permute(PPc(:,1,:), [1 3 2]);
Ycp = permute(PPc(:,2,:), [1 3 2]);
Zcp = permute(PPc(:,3,:), [1 3 2]);
Xcq = permute(QQc(:,1,:), [1 3 2]);
Ycq = permute(QQc(:,2,:), [1 3 2]);
Zcq = permute(QQc(:,3,:), [1 3 2]);
% If there are parallel lines, find the distances  between them
% Note, that for parallel lines, the closest points will be undefined
% (will contain NaN's)
if(any(idx_par))
    idx_par3 = repmat(idx_par, [1,3,1]); % logical indices
    PPc(idx_par3) = PP0(idx_par3);
    tmpl = repmat(dd(idx_par)./bb(idx_par), [1, 3, 1]);
    tmp2 = vv(find(idx_par3));
    
    QQc(idx_par3) = QQ0(idx_par3) + tmpl(:).*tmp2;
end
PQc = (PPc - QQc);
D = permute(sqrt(dot(PQc,PQc,2)), [1 3 2]);
[Dmin, idx_min] = min(D(:));
[imin,jmin] = ind2sub(size(D), idx_min);
end

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
