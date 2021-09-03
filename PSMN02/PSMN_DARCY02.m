
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


load(strcat(folderScriptshell,'all_IFPEN_DARCY02_experiments.mat'))

allresults = struct();
% CHECK WITH FLOR : cd(allExpeStrct(iexpe).analysisFolder)
% cd(folderExperiment)
% file_log_ID = fopen(sprintf('log_%s.txt',allExpeStrct(iexpe).name), 'a');

% load centers
allTrajLOAD = load(strcat(folderExperiment,nameAllTraj));
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
    
%     % STEP 3 - removing the NaNs for all t
%     for it = 1 : size(CC1,2)
%         ikill = [];
%         for ip = 1 : size(CC1(it).X,1)
%             if isnan(CC1(it).X(ip)) || isnan(CC1(it).Y(ip))
%                 ikill = [ikill,ip];
%             end
%         end
%         CC1(it).X(ikill) = [];
%         CC1(it).Y(ikill) = [];
%         clear ikill
%         ikill = [];
%         for ip = 1 : size(CC2(it).X,1)
%             if isnan(CC2(it).X(ip)) || isnan(CC2(it).Y(ip))
%                 ikill = [ikill,ip];
%             end
%         end
%         CC2(it).X(ikill) = [];
%         CC2(it).Y(ikill) = [];
%     end
    
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
    
    wti = 300; % width template images
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
    transformationType = 'affine';
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
    Ttype = 'T3';
    
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
                
                [crossP,D] = crossRaysonFire(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype);
                if length(crossP)>0
                    someTrajectories(iselTraj).x3D(ixy) = crossP(1);
                    someTrajectories(iselTraj).y3D(ixy) = crossP(2);
                    someTrajectories(iselTraj).z3D(ixy) = crossP(3);
                    someTrajectories(iselTraj).D(ixy) = D;
                end
            end
        end
        
        
    end
    
    %%%%%%%%%
    hist3D = [];
    histD  = [];
    % plot the result
    for itraj3D = 1 : size(someTrajectories,2)
        hist3D = [hist3D,[someTrajectories(itraj3D).z3D]];
        histD  = [histD, [someTrajectories(itraj3D).D]];
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
    
    save(strcat(folderExperiment,'allResults_auFilDeLEau.mat'),'allresults')
    
    cf = clock;
    fprintf('\t \t \t it took %4.0f s \n',etime(cf,ci))
    
    cPlane_f = clock;
    % fprintf(file_log_ID, 'just run plane %0.3d in %0.0f s \n', iplane , etime(cPlane_f,cPlane_i) );
    fprintf('just run plane %0.3d in %0.0f s \n', iplane , etime(cPlane_f,cPlane_i) );
    
end

fprintf('end \n')
% fclose(file_log_ID)

end

