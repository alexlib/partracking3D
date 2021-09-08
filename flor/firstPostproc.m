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
dt = 1/fps; %difference between two successive images
%%
index = 30;
oneTraj = allresults(index).someTrajectories;
trj = 18;
x = oneTraj(trj).x3D;
y = oneTraj(trj).y3D;
z = oneTraj(trj).z3D;
%% diferencias forward
for t = 1:(length(x)-1)
    vxF(t) = (x(t+1)-x(t))/dt;
    vyF(t) = (y(t+1)-y(t))/dt;
    vzF(t) = (z(t+1)-z(t))/dt;
end

%% diferencia centrada
for t = 2:(length(x)-1)
    vxC(t) = (x(t+1)-x(t-1))/dt;
    vyC(t) = (y(t+1)-y(t-1))/dt;
    vzC(t) = (z(t+1)-z(t-1))/dt;
end
%% 
figure; hold on
plot3(vxF,vyF,vzF)
plot3(vxC,vyC,vzC)
legend('F','C')
