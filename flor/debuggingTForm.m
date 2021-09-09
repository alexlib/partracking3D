%"bad" trajectory
xPoint = 15.47;
yPoint = -10.16;
zPoint = 26;

%"good" traj
xPoint = 3.32;
yPoint = -1.35;
zPoint = 18.86;

xPoint = -7.4;
yPoint = -17.82;
zPoint = 16.;

%these dont look to be the same
xPoint = -24.0888;
yPoint = 11.8934;
zPoint = 21.3539;

%seems to be good
xPoint = -5.5315;
yPoint = -4.24088;
zPoint = 18.3831;

%bad
xPoint = -4.82721;
yPoint = -5.14635;
zPoint = 18.5414;

%bad
xPoint = 1.93945;
yPoint = 9.45962;
zPoint = 23.073;

for itrck = 1 : length(allresults(iplane).someTrajectories) 
   clear X Y Z
   X = allresults(iplane).someTrajectories(itrck).x3D;
   Y = allresults(iplane).someTrajectories(itrck).y3D;
   Z = allresults(iplane).someTrajectories(itrck).z3D;
   dloc(itrck) = min(sqrt((X-xPoint).^2+(Y-yPoint).^2+(Z-zPoint).^2));
end
[~,trck] = min(dloc);

it1 = allresults(iplane).someTrajectories(trck).itraj1;
it2 = allresults(iplane).someTrajectories(trck).itraj2;

Xb = allresults(iplane).trajArray_CAM1(it1).track(:,1);
Yb = allresults(iplane).trajArray_CAM1(it1).track(:,2);

Xr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,1);
Yr = allresults(iplane).trajArray_CAM2RAW(it2).track(:,2);

figure;
sgtitle(sprintf('Trajectory %0.3d',trck))
subplot(1,3,1)
plot(Xb,Yb,'-b','lineWidth',4)
title('CAM1')
xlabel('x')
ylabel('y')
subplot(1,3,2)
plot(Xr,Yr,'-r','lineWidth',4)
title('CAM2')
xlabel('x')
ylabel('y')
subplot(1,3,3)
box on
view(3)
X = allresults(iplane).someTrajectories(trck).x3D;
Y = allresults(iplane).someTrajectories(trck).y3D;
Z = allresults(iplane).someTrajectories(trck).z3D;
plot3(X,Y,Z,'lineWidth',4)
xlabel('x')
ylabel('y')
zlabel('z')
title('Matched trajs')

%%
%I choose some by hand
% iexpe=8; iplane=31

badTracks = [112,132,218];
goodTracks = [];

badTransf = struct();
goodTransf = struct();

j=1
for i = badTracks
    j = j+1
    badTransf(j) = allresults(iplane).tform1.T
end

j=1
for i = goodTracks
    j = j+1
    goodTransf(j) = allresults(iplane).tform1.T
end

 

