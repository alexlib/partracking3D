%% workflow_mcin2_prepare_for_PSMN

cd('C:\Users\darcy\Desktop\git\Robust-Estimation')
load('all_IFPEN_DARCY02_experiments.mat')

% ho to load calib file in the structure
% allExpeStrct(7).calib = load('D:\IFPEN\analysisExperiments\calibFiles\calib_2021_08_13\calib.mat');
iexpe = 9; 

%% find Centers 
allresults = struct();
cd(allExpeStrct(iexpe).analysisFolder)
file_log_ID = fopen(sprintf('log_%s.txt',allExpeStrct(iexpe).name), 'a');
 
planeI = 001;
planeF = 060;

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
            DARCY02_findCC(allExpeStrct,iexpe,iSeq,maxdist,longmin,'figures','yes');
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

%% visualise alltraj
% choose time
planeI=34
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
nameMat = 'alltraj_2021_08_15_electroVavle_at_40percent.mat';
%nameMat = 'alltraj_2021_05_05.mat';
save(nameMat,'allTraj','-v7.3')

%% FLOR
%modified lines 145-185
%created new variable CCRAW_loc (localized)
%in which we only see a selected area (see in figs)
%we can choose a better criteria

%% DARCY02_findTracks - function part

function [CCall,CCout,ImMean,ImMax,filename] = DARCY02_findCC(allExpeStrct,iexpe,ifile,maxdist,longmin,varargin)

% 1. load image
% 2. subtract mean of the image sequence
% 3. find particles positions on all images

dofigures = 'yes';
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
Nwidth = 1.;
CCRAW = struct();
CCRAW_loc = struct();

switch dofigures
    case 'yes'
        hmax = figure('defaultAxesFontSize',20);
        imshow(imadjust(ImMax,[1 15]/255))
        colorTime = jet(size(M,3));
end

sizeX = size(Im01,2);
sizeY = size(Im01,1);
centerX = sizeX/2;
centerY = sizeY/2;

for it = 1 : size(M,3)
    Im = zeros(sizeY,sizeX,class(Im01));
    Im(:,:) = Im01(:,:,it);
    
    CCRAW(it).xyRAW = pkfnd(Im01(:,:,it),th,sz);
    
    minx = centerX-1*std(CCRAW(it).xyRAW(:,1));
    miny = centerY-1*std(CCRAW(it).xyRAW(:,2));
    maxx = centerX+1*std(CCRAW(it).xyRAW(:,1));
    maxy = centerY+1*std(CCRAW(it).xyRAW(:,2));
    
    condXmin = CCRAW(it).xyRAW(:,1)>minx;
    condXmax = CCRAW(it).xyRAW(:,1)<maxx;
    condYmin = CCRAW(it).xyRAW(:,2)>miny;
    condYmax = CCRAW(it).xyRAW(:,2)<maxy;
    CCRAW_loc(it).xy(:,1) = CCRAW(it).xyRAW(condXmin&condXmax&condYmin&condYmax,1);
    CCRAW_loc(it).xy(:,2) = CCRAW(it).xyRAW(condXmin&condXmax&condYmin&condYmax,2);

    
    %refine at subpixel precision
    for ixy = 1 : size(CCRAW_loc(it).xy,1)
        clear xpkfnd ypkfnd Ip
        Ip = zeros(2*Nwidth+1,2*Nwidth+1,'double');
        xpkfnd = CCRAW_loc(it).xy(ixy,1);
        ypkfnd = CCRAW_loc(it).xy(ixy,2);
        Ip = double(Im(ypkfnd-Nwidth:ypkfnd+Nwidth,xpkfnd-Nwidth:xpkfnd+Nwidth));
        % replace all 0 Ip values by a tiny value
        list0 = find(Ip==0);
        Ip(list0) = 1e-6;
        
        CC(it).xy(ixy,1) = xpkfnd + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
        CC(it).xy(ixy,2) = ypkfnd + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
    end
    


    switch dofigures
        case 'yes'
    if size(CCRAW_loc(it).xy,1) >0
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
        
        if size(CCRAW_loc(it).xy,1) >0
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

