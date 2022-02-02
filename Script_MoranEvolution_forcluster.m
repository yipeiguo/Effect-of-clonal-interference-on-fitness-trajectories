% Script to carry out full population simulation in a Moran process.
% This uses the function 'EvolveMoran_Gillespie.m' (or other versions)


close all; clear variables;
rng('shuffle')

tic

%% Define parameters
% model parameters:
N = 1e6; modelparams.N = N;
mub0 = 1e-8; %5e-5;
WTF = 1;
alpha0 = 1/0.05; % initial inverse mean fitness effect
k = 1; % 1 for moran, 2 for WF, 2*log(D)/(D-1) for dilution protocol
A = N*log(N)*k;
modelparams.A = A;
modelparams.k = k;

% functions defining macroscopic epistasis
% possible options include house of cards('HoC'), non-epistatic('NEPI'),
% stairway to heaven('STH'), diminishing returns('DR'), 
% replenishing/releasing mutations ('RM') and others ('others') in which 
% case we will need to manually specify the functional form.
landscape = 'DR'; 
if strcmp(landscape,'HoC')
    alphafunc = @(x) alpha0.*x;
    mubfunc = @(x) mub0.*exp(-(x-1));
elseif strcmp(landscape,'NEPI')
    alphafunc = @(x) alpha0.*x;
    mubfunc = @(x) mub0;
elseif strcmp(landscape,'STH')
    alphafunc = @(x) alpha0;
    mubfunc = @(x) mub0;
elseif strcmp(landscape,'DR')
    g = 5;
    alphafunc = @(x) alpha0.*x.^g;
    mubfunc = @(x) mub0;
elseif strcmp(landscape,'RM')
    alphafunc = @(x) alpha0.*exp(x-1);
    mubfunc = @(x) mub0.*exp(x-1);
end
modelparams.alphafunc = alphafunc;
modelparams.mubfunc = mubfunc;

% experimental/simulation set-up
% tstoreVec = 1:1:10;
% tstoreVec = [1:1:10,20:10:100,200:100:1e3,2e3:1e3:1e4];
% tstoreVec = [1:1:10,20:20:1e4];
% tstoreVec = [1:1:10,20:20:1000,1200:200:1e5];
tstoreVec = [20:20:100,200:200:1e3,2e3:2e3:1e6];
numreplicates = 1;  % number of replicates of each experiment

disp('Parameters defined.');
toc

%% Compile into structure initialpop

initialpop.gvec = WTF;
initialpop.Nvec = N;
initialpop.smat = [];
initialpop.sEtimemat = [];
initialpop.NumMutvec = 0;

otherparams.tstoreVec = tstoreVec;


%% file name for storage
runindx = 1;
fn = strcat('MoranEvolSimv2Data_N',num2str(N),'_mub0_',num2str(mub0),'_alpha0_',...
    num2str(alpha0),'_maxT_',num2str(max(tstoreVec)),'_landscape_',landscape,...
    '_numreplicates',num2str(numreplicates),'_run',num2str(runindx));
fn = strrep(fn,'.','pt');

basenamestr_temp = strcat('MoranEvolSimv2Data_N',num2str(N),'_mub0_',num2str(mub0),'_alpha0_',...
    num2str(alpha0),'_maxT_',num2str(max(tstoreVec)),'_landscape_',landscape);
basenamestr_temp = strrep(basenamestr_temp,'.','pt');

%% Define Storage vectors
gavtrajmat = ones(numreplicates,length(tstoreVec));
numgrpsmat = ones(numreplicates,length(tstoreVec));
nummutavtrajmat = zeros(numreplicates,length(tstoreVec));
popTrajsCell = cell(1,numreplicates);
finalPopCells = cell(1,numreplicates);

%% Main Simulation
tic
for repindx = 1:numreplicates
    fprintf('Trial index =%d \n',repindx)
    namestr_temp = strcat(basenamestr_temp,'_rep',num2str(repindx));
    
    [gavTraj,nummutavTraj,numgrpsTraj,allpopsCell,finalpop] = ...
        EvolveMoran_Gillispie_v2(initialpop,modelparams,otherparams,namestr_temp);
    gavtrajmat(repindx,:) = gavTraj;
    numgrpsmat(repindx,:) = numgrpsTraj;
    nummutavtrajmat(repindx,:) = nummutavTraj;
    popTrajsCell{repindx} = allpopsCell;
    finalPopCells{repindx} = finalpop;
    save(fn);
    toc

end
toc

save(fn);













