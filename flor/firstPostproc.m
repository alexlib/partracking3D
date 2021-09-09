%%
%Getting familiarized with the trajectories
%%
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
%%
fps = 200; 
dt = 1/fps; %difference between two successive images, in s
%%
index = 1;
oneTraj = allresults(iplane).someTrajectories;
trj = 36;
x = oneTraj(trck).x3D; %mm
y = oneTraj(trck).y3D; %mm
z = oneTraj(trck).z3D; %mm


figure; hold on
plot3(x,y,z,'.-')
plot3(medfilt1(x),medfilt1(y),medfilt1(z),'.-')
legend('raw','filtered')
xlabel('x')
ylabel('y')
zlabel('z')
title('trajectory')
%% DIFFERENT FILTERED TRAJECTORIES

%% DIFFERENT VELOCITY INCREMENTS

%% 2nd order accu
for t = 1:(length(x)-1)
    vxF(t) = (x(t+1)-x(t))/dt;
    vyF(t) = (y(t+1)-y(t))/dt;
    vzF(t) = (z(t+1)-z(t))/dt;
end

%% 2nd order accu. 
%THIS ONE WORKS BETTER
for t = 2:(length(x)-1)
    vxC(t) = (x(t+1)-x(t-1))/(2*dt);
    vyC(t) = (y(t+1)-y(t-1))/(2*dt);
    vzC(t) = (z(t+1)-z(t-1))/(2*dt);
end

%% 3rd order accu
for t = 3:(length(x)-1)
    vxT(t) = (2*x(t+1)+3*x(t)-6*x(t-1)+x(t-2))/(6*dt);
    vyT(t) = (2*y(t+1)+3*y(t)-6*y(t-1)+y(t-2))/(6*dt);
    vzT(t) = (2*z(t+1)+3*z(t)-6*z(t-1)+z(t-2))/(6*dt);
end

%% 
figure; hold on
%plot3(vxF,vyF,vzF,'.-')
%plot3(vxC,vyC,vzC,'.-')
plot(vxF,'.-')
plot(vxC,'.-')
plot(vxT,'.-')
legend('F','C','T')
xlabel('x')
ylabel('y')
zlabel('z')