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

%"good" traj
xPoint = 3.32;
yPoint = -1.35;
zPoint = 18.86;

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

%plot variability in each tranf coordinate
for i = 1:length(allresults)
xx(i) = allresults(i).tform1.T(1,1);
yy(i) = allresults(i).tform1.T(2,2);
zz(i) = allresults(i).tform1.T(3,3);
xy(i) = allresults(i).tform1.T(1,2);
%xz(i) = allresults(i).tform1.T(1,3); 
yx(i) = allresults(i).tform1.T(2,1);
%yz(i) = allresults(i).tform1.T(2,3);
zx(i) = allresults(i).tform1.T(3,1); 
zy(i) = allresults(i).tform1.T(3,2); 
end

colors=jet(7);
figure;
hold on
plot(xx,'-o','color',colors(1,:),'markerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:))
plot(yy,'-o','color',colors(2,:),'markerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))
plot(zz,'-o','color',colors(3,:),'markerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:))
plot(xy,'-o','color',colors(4,:),'markerFaceColor',colors(4,:),'MarkerEdgeColor',colors(4,:))
%plot(xz,'-o','color',colors(5,:),'markerFaceColor',colors(5,:),'MarkerEdgeColor',colors(5,:))
plot(yx,'-o','color',colors(5,:),'markerFaceColor',colors(5,:),'MarkerEdgeColor',colors(5,:))
%plot(yz,'-o','color',colors(7,:),'markerFaceColor',colors(7,:),'MarkerEdgeColor',colors(7,:))
plot(zx,'-o','color',colors(6,:),'markerFaceColor',colors(6,:),'MarkerEdgeColor',colors(6,:))
plot(xz,'-o','color',colors(7,:),'markerFaceColor',colors(7,:),'MarkerEdgeColor',colors(7,:))
legend('xx','yy','zz','xy','yx','zx','zy')
xlabel('Plane')
ylabel('Transformation')



 

