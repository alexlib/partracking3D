
ifile = 1;

% 1. load image
% 2. subtract mean of the image sequence
% 3. find particles positions on all images

dofigures = 'no';

iexpe = 2;
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

%%
figure
imshow(imadjust(M(:,:,1),[1 10]/255))

figure
imshow(imadjust(ImMax,[1 15]/255))

figure
imshow(imadjust(Im01(:,:,1),[1 10]/255))
%%
clear CC
CC = struct();
hmax = figure('defaultAxesFontSize',20,'position',[ 473 172  1051 794]);
imshow(imadjust(ImMax,[1 15]/255))

colorTime = jet(250);
for it = 1 : 64
    
    if exist('hp')
        if isvalid(hp)
            hp.MarkerFaceColor = 'none';
            hp.MarkerEdgeColor = colorTime(it,:);
        end
    end
    Im = zeros(size(Im01,1),size(Im01,2),class(Im01));
    Im(:,:) = Im01(:,:,it);
    
    CCRAW(it).xyRAW = pkfnd(Im01(:,:,it),th,sz);
    
    % fix the problem when there are 0 values in image Ip !!!
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
    CCwrk = CC(it).xy;
    figure(hmax), title(sprintf('time: %4.0f - nparts: %4.0f',it,length(CCwrk)))
    hold on
    hp = plot(CCwrk(:,1),CCwrk(:,2),'o','markerFaceColor',colorTime(it,:),...
        'markerEdgeColor','none','markerSize',5);
   
    
end

%%
figure
imshow(imadjust(M(:,:,it),[1 10]/255))
 
%%
iexpe = 7;
ifile = 390;
clear CC

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
th = allExpeStrct(iexpe).centerFinding_th
sz = allExpeStrct(iexpe).centerFinding_sz
Nwidth = 1;
CCRAW = struct();

hmax = figure('defaultAxesFontSize',20);
imshow(imadjust(ImMax,[1 15]/255))
colorTime = jet(size(M,3));

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
    if size(CCRAW(it).xyRAW,1) >0
        CCwrk = CC(it).xy;
        figure(hmax), title(sprintf('time: %4.0f - nparts: %4.0f',it,length(CCwrk)))
        hold on
        hp = plot(CCwrk(:,1),CCwrk(:,2),'o','markerFaceColor',colorTime(it,:),...
            'markerEdgeColor','none','markerSize',5);
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
            
            fprintf('ifile: %4.0f - it: %4.0f - n part found: %4.0f n_found_RAW: %4.0f \n',...
                ifile,it,length(CCX),size(CC(it).xy,1))
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
for it = 1 : size(M,3)
    fprintf('ifile is %4.0f and it is: %4.0f\n',ifile,it)
    if ((length(CC) >= it)) && size(CC(it).xy,1) > 0
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