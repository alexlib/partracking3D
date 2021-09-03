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

