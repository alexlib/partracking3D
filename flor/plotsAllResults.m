cd('C:\Users\darcy\Desktop\git\Robust-Estimation')
load('all_IFPEN_DARCY02_experiments.mat')

iexpe = 8;
folderScriptshell = allExpeStrct(iexpe).analysisFolder;
folderExperiment  = folderScriptshell;
nameAllTraj = 'alltraj_2021_08_15_electroVavle_at_40percent.mat';

cd(folderExperiment)
nameResults = 'allResults_auFilDeLEau.mat';
load(nameResults)

planeI = 1;    
planeF = length(allresults);

%% Display results
% all tracks
iplane = 1;

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
title('All points')
xlabel('x')
ylabel('y')

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
Xb = allresults(iplane).trajArray_CAM1(it2).track(:,1);
Yb = allresults(iplane).trajArray_CAM1(it2).track(:,2);


xred = 1970;
yred = 700;
clear dloc
for itrck = 1 : length(allresults(iplane).trajArray_CAM2RAW) 
   clear X Y 
   X = allresults(iplane).trajArray_CAM2RAW(itrck).track(:,1);
   Y = allresults(iplane).trajArray_CAM2RAW(itrck).track(:,2);
   dloc(itrck) = min(sqrt((X-xred).^2+(Y-yred).^2));
end
[~,it2] = min(dloc);
clear Xr Yr
Xr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,1);
Yr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,2);

plot(Xb,Yb,'-b','lineWidth',4)
plot(Xr,Yr,'-r','lineWidth',4)
title('Two matched trajectories')
% identify paired tracks and non-paired tracks


%% matched tracks
figure, hold on, box on, view(3)
xlabel('x')
ylabel('y')
zlabel('z')
for iplane = 31%planeI : planeF
 for itrck = 1 : length(allresults(iplane).someTrajectories)
    clear X Y
    X = allresults(iplane).someTrajectories(itrck).x3D;
    Y = allresults(iplane).someTrajectories(itrck).y3D;
    Z = allresults(iplane).someTrajectories(itrck).z3D;
    plot3(X,Y,Z,'lineWidth',4)
 end
end
title('All matched tracks')
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
    xy01(i).x = CalibFileCam1(i).pos3D(:,1);
    xy01(i).y = CalibFileCam1(i).pos3D(:,2);
    xy02(i).x = CalibFileCam2(i).pos3D(:,1);
    xy02(i).y = CalibFileCam2(i).pos3D(:,2);
end

colorTime = jet(length(CalibFileCam1));
clear x3D y3D z3D
%figure('defaultAxesFontSize',20), box on
%xlabel('x')
%ylabel('y')
%zlabel('z')
%hold on
for ipoints = 2%:length(CalibFileCam1)
    clear x3D y3D z3D
for ixy = 1 : length(xy01(ipoints).x)
    x_pxC1 = xy01(ipoints).x(ixy);
    y_pxC1 = xy01(ipoints).y(ixy);
    x_pxC2 = xy02(ipoints).x(ixy);
    y_pxC2 = xy02(ipoints).y(ixy);

    [crossP,D] = crossRaysonFire(CalibFileCam1,CalibFileCam2,x_pxC1,y_pxC1,x_pxC2,y_pxC2,Ttype);
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


%% matched tracks
% test by Flor. Will filter out the tracks that fall 
% outside the deviation

figure, hold on, box on, view(3)
xlabel('x')
ylabel('y')
zlabel('z')
iplane = 31;

 for itrck = 1 : length(allresults(iplane).someTrajectories)
     minz = mean(allresults(iplane).hist3D)-2.335*std(allresults(iplane).hist3D);
     maxz = mean(allresults(iplane).hist3D)+2.335*std(allresults(iplane).hist3D);
     minx = mean(allresults(iplane).histx)-2.335*std(allresults(iplane).histx);
     maxx = mean(allresults(iplane).histx)+2.335*std(allresults(iplane).histx);
     miny = mean(allresults(iplane).histy)-2.335*std(allresults(iplane).histy);
     maxy = mean(allresults(iplane).histy)+2.335*std(allresults(iplane).histy);
    if mean(allresults(iplane).someTrajectories(itrck).z3D) > minz && mean(allresults(iplane).someTrajectories(itrck).z3D) < maxz ...
            && mean(allresults(iplane).someTrajectories(itrck).y3D) > miny && mean(allresults(iplane).someTrajectories(itrck).y3D) < maxy ...
            && mean(allresults(iplane).someTrajectories(itrck).x3D) > minx && mean(allresults(iplane).someTrajectories(itrck).x3D) < maxx
        clear X Y
        X = allresults(iplane).someTrajectories(itrck).x3D;
        Y = allresults(iplane).someTrajectories(itrck).y3D;
        Z = allresults(iplane).someTrajectories(itrck).z3D;
        plot3(X,Y,Z,'lineWidth',4)
    end
 end
