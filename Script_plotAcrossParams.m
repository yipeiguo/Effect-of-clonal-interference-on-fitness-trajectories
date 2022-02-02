% This is a script to plot trajectories for different sets of parameters on
% the same plot

close all; 
clear variables;

%% load data
whichlandscape = 'HoC'; % 'DR' or 'HoC' or 'RM'
pathname = strcat('./',whichlandscape,'/');
Sdir = dir(strcat(pathname,'*.mat'));
filenames = {Sdir.name};
numfiles = length(filenames);

% specify parameters common to all files
alpha0 = 100;

% Extract parameters of different data sets & the actual data for plotting
% storage arrays for parameters:
N_allsets = zeros(1,numfiles);
mub0_allsets = zeros(1,numfiles);
maxT_allsets = zeros(1,numfiles);
numtrials_allsets = zeros(1,numfiles);
% storage arrays for simulation data:
tstoreVec_cell = cell(1,numfiles);
gavTraj_cell = cell(1,numfiles);
gavTraj_ub_cell = cell(1,numfiles);
gavTraj_lb_cell = cell(1,numfiles);
gTrajMat_cell = cell(1,numfiles);
mavTraj_cell = cell(1,numfiles);
mavTraj_ub_cell = cell(1,numfiles);
mavTraj_lb_cell = cell(1,numfiles);
mTrajMat_cell = cell(1,numfiles);
% storage arrays for analytical expression:
tvec_cell = cell(1,numfiles);
Fvec_cell = cell(1,numfiles);
Svec_cell = cell(1,numfiles);
tvec_noCI_cell = cell(1,numfiles);
Fvec_noCI_cell = cell(1,numfiles);
Svec_noCI_cell = cell(1,numfiles);
% storage arrays for other useful quantities:
CIintegralF0_allsets = zeros(1,numfiles); % for calculating initial Fgrad
CIintegralS0_allsets = zeros(1,numfiles); % for calculating initial Sgrad
for fileIndx = 1:numfiles
    fn = filenames{fileIndx};
    N_pos = strfind(fn,'N');
    mu_pos = strfind(fn,'mub0');
    alpha_pos = strfind(fn,'alpha');
    maxT_pos = strfind(fn,'maxT');
    landscape_pos = strfind(fn,'landscape');
    numtrials_pos = strfind(fn,'numtrials');
    set_pos = strfind(fn,'set');
    % Extract parameters
    N_allsets(fileIndx) = str2double(fn(N_pos+1:mu_pos-2));
    mub0_allsets(fileIndx) = str2double(fn(mu_pos+5:alpha_pos-2));
    maxT_allsets(fileIndx) = str2double(fn(maxT_pos+5:landscape_pos-2));
    numtrials_allsets(fileIndx) = str2double(fn(numtrials_pos+9:set_pos-2));
    
    % Extract data
    load(strcat(pathname,fn),'tstoreVec_CI','gavtraj_CI_mean','gavtraj_CI_ub','gavtraj_CI_lb',...
        'nummutavtraj_CI_mean','nummutavtraj_CI_ub','nummutavtraj_CI_lb',...
        'tvec_CI','Fvec_CI','Svec_CI','tvec_noCI','Fvec_noCI','Svec_noCI',...
        'CIintegral_Finit','CIintegral_Sinit','gavtrajmat_CI','nummutavtrajmat_CI');
    tstoreVec_cell{fileIndx} = tstoreVec_CI;
    gavTraj_cell{fileIndx} = gavtraj_CI_mean;
    gavTraj_ub_cell{fileIndx} = gavtraj_CI_ub;
    gavTraj_lb_cell{fileIndx} = gavtraj_CI_lb;
    gTrajMat_cell{fileIndx} = gavtrajmat_CI;
    mavTraj_cell{fileIndx} = nummutavtraj_CI_mean;
    mavTraj_ub_cell{fileIndx} = nummutavtraj_CI_ub;
    mavTraj_lb_cell{fileIndx} = nummutavtraj_CI_lb;
    mTrajMat_cell{fileIndx} = nummutavtrajmat_CI;
    tvec_cell{fileIndx} = tvec_CI;
    Fvec_cell{fileIndx} = Fvec_CI;
    Svec_cell{fileIndx} = Svec_CI;
    tvec_noCI_cell{fileIndx} = tvec_noCI;
    Fvec_noCI_cell{fileIndx} = Fvec_noCI;
    Svec_noCI_cell{fileIndx} = Svec_noCI;
    CIintegralF0_allsets(fileIndx) = CIintegral_Finit;
    CIintegralS0_allsets(fileIndx) = CIintegral_Sinit;
end

% calculate theoretical initial gradients
FgradInit_noCI_allsets = 2.*N_allsets.*mub0_allsets./(alpha0^2);
FgradInit_CI_allsets = CIintegralF0_allsets.*N_allsets.*mub0_allsets./(alpha0^2);
SgradInit_noCI_allsets = N_allsets.*mub0_allsets./alpha0;
SgradInit_CI_allsets = CIintegralS0_allsets.*N_allsets.*mub0_allsets./alpha0;

%% Extract sets with a given N 
NoI = 1e6;
NindsOI = find(N_allsets == NoI);
mu_rel = mub0_allsets(NindsOI);
[mu_sorted,mu_I] = sort(mu_rel);
setInds_Nfix = NindsOI(mu_I);
numsets_Nfix = length(NindsOI);
if numsets_Nfix > 1
    % extract common time points between all sets
    tstored_common = intersect(tstoreVec_cell{setInds_Nfix(1)},tstoreVec_cell{setInds_Nfix(2)});
    for kk = 3:numsets_Nfix
        tstored_common = intersect(tstored_common,tstoreVec_cell{setInds_Nfix(kk)});
    end

    % coarsely-sample these trajectories and scale time accordingly
    deltat = 100; % sample at most once in each deltat time interval
    tCoarse_Nfix = cell(1,numsets_Nfix);
    gavTrajCoarse_Nfix = cell(1,numsets_Nfix);
    gubCoarse_Nfix = cell(1,numsets_Nfix);
    glbCoarse_Nfix = cell(1,numsets_Nfix);
    mavTrajCoarse_Nfix = cell(1,numsets_Nfix);
    mubCoarse_Nfix = cell(1,numsets_Nfix);
    mlbCoarse_Nfix = cell(1,numsets_Nfix);
    for jj = 1:numsets_Nfix
        tstoreVec = tstoreVec_cell{setInds_Nfix(jj)};
        [~,iA] = unique(floor(tstoreVec./deltat));
        tstoreNew = [0,tstoreVec(iA(2:end))];
        tCoarse_Nfix{jj} = tstoreNew;
        gavTraj = gavTraj_cell{setInds_Nfix(jj)};
        gavTrajCoarse_Nfix{jj} = [1,gavTraj(iA(2:end))];
        gub = gavTraj_ub_cell{setInds_Nfix(jj)};
        glb = gavTraj_lb_cell{setInds_Nfix(jj)};
        gubCoarse_Nfix{jj} = [1,gub(iA(2:end))];
        glbCoarse_Nfix{jj} = [1,glb(iA(2:end))];
        mavTraj = mavTraj_cell{setInds_Nfix(jj)};
        mavTrajCoarse_Nfix{jj} = [0,mavTraj(iA(2:end))];
        mub = mavTraj_ub_cell{setInds_Nfix(jj)};
        mlb = mavTraj_lb_cell{setInds_Nfix(jj)};
        mubCoarse_Nfix{jj} = [0,mub(iA(2:end))];
        mlbCoarse_Nfix{jj} = [0,mlb(iA(2:end))];
    end

    % obtain F gradient based on when trajectories first hit some target
    % fitness
    Fhit_scan = [1.01,1.05,1.1,1.2,1.3,1.4];
    initFgradMat_Nfix = zeros(numsets_Nfix,length(Fhit_scan));
    for jj = 1:numsets_Nfix
        tstoreVec = tCoarse_Nfix{jj};
        gavTraj = gavTrajCoarse_Nfix{jj};
        for kk = 1:length(Fhit_scan)
            Fhit = Fhit_scan(kk);
            tIndx = find(gavTraj>=Fhit,1);
            initgrad = (gavTraj(tIndx)-1)/tstoreVec(tIndx);
            initFgradMat_Nfix(jj,kk) = initgrad;
        end
    end
    
    % obtain S gradient based on when trajectories first hit some target
    % fitness
    Shit_scan = [5,10,15,20];
    initSgradMat_Nfix = zeros(numsets_Nfix,length(Shit_scan));
    for jj = 1:numsets_Nfix
        tstoreVec = tCoarse_Nfix{jj};
        mavTraj = mavTrajCoarse_Nfix{jj};
        for kk = 1:length(Shit_scan)
            Shit = Shit_scan(kk);
            tIndx = find(mavTraj>=Shit,1);
            initgrad = mavTraj(tIndx)/tstoreVec(tIndx);
            initSgradMat_Nfix(jj,kk) = initgrad;
        end
    end
end

% plot scaled trajectories
if numsets_Nfix > 1
    colormat_Nconstant = rand(numsets_Nfix,3);
    namecell_Nconstant = cell(1,numsets_Nfix);
    numrows = ceil(sqrt(length(Fhit_scan)));
    numcols = ceil(length(Fhit_scan)/numrows);
    figure;
    for kk = 1:length(Fhit_scan)
        subplot(numrows,numcols,kk);
        for jj = 1:numsets_Nfix
            mub0 = mu_sorted(jj);
            tstoreVec = tCoarse_Nfix{jj};
            initgrad = initFgradMat_Nfix(jj,kk);
            plot(tstoreVec.*initgrad,gavTrajCoarse_Nfix{jj},'-.',...
                'color',colormat_Nconstant(jj,:),'LineWidth',1.5);
            hold on
            namecell_Nconstant{jj} = strcat('\mu = ',num2str(mub0));
        end
        for jj = 1:numsets_Nfix
            tstoreVec = tCoarse_Nfix{jj};
            initgrad = initFgradMat_Nfix(jj,kk);
            fill([tstoreVec.*initgrad, fliplr(tstoreVec.*initgrad)],...
                [gubCoarse_Nfix{jj},fliplr(glbCoarse_Nfix{jj})],...
                colormat_Nconstant(jj,:),'linestyle', 'none','FaceAlpha',0.4);
            hold on
        end
%         for jj = 1:numsets_Nfix
%             setIndx = setInds_Nfix(jj);
%             initgrad = initFgradMat_Nfix(jj,kk);
%             plot(tvec_cell{setIndx}.*initgrad,Fvec_cell{setIndx},...
%                 'color',colormat_Nconstant(jj,:),'LineWidth',1.5);
%             hold on
%         end
%         xlim([0 1]);
        xlabel('t scaled');
        ylabel('F');
        title(strcat('Fhit = ',num2str(Fhit_scan(kk))));
        if kk == 1
            legend(namecell_Nconstant,'location','best');
        end        
    end
    
    numrows = ceil(sqrt(length(Shit_scan)));
    numcols = ceil(length(Shit_scan)/numrows);
    figure;
    for kk = 1:length(Shit_scan)
        subplot(numrows,numcols,kk);
        for jj = 1:numsets_Nfix
            mub0 = mu_sorted(jj);
            tstoreVec = tCoarse_Nfix{jj};
            initgrad = initSgradMat_Nfix(jj,kk);
            plot(tstoreVec.*initgrad,mavTrajCoarse_Nfix{jj},'-.',...
                'color',colormat_Nconstant(jj,:),'LineWidth',1.5);
            hold on
            namecell_Nconstant{jj} = strcat('\mu = ',num2str(mub0));
        end
        for jj = 1:numsets_Nfix
            tstoreVec = tCoarse_Nfix{jj};
            initgrad = initSgradMat_Nfix(jj,kk);
            fill([tstoreVec.*initgrad, fliplr(tstoreVec.*initgrad)],...
                [mubCoarse_Nfix{jj},fliplr(mlbCoarse_Nfix{jj})],...
                colormat_Nconstant(jj,:),'linestyle', 'none','FaceAlpha',0.4);
            hold on
        end
        for jj = 1:numsets_Nfix
            setIndx = setInds_Nfix(jj);
            initgrad = initSgradMat_Nfix(jj,kk);
            plot(tvec_cell{setIndx}.*initgrad,Svec_cell{setIndx},...
                'color',colormat_Nconstant(jj,:),'LineWidth',1.5);
            hold on
        end
        xlim([0 200]);
        xlabel('t scaled');
        ylabel('S');
        title(strcat('Shit = ',num2str(Shit_scan(kk))));
        if kk == 1
            legend(namecell_Nconstant,'location','best');
        end
        
    end
end

%% Extract sets with a given mu
muOI = 1e-5;
muIndsOI = find(mub0_allsets == muOI);
N_rel = N_allsets(muIndsOI);
[N_sorted,N_I] = sort(N_rel);
setInds_mufix = muIndsOI(N_I);
numsets_mufix = length(muIndsOI);

if numsets_mufix > 1
    % coarsely-sample these trajectories and scale time accordingly
    deltat = 200; % sample at most once in each deltat time interval
    tCoarse_mufix = cell(1,numsets_mufix);
    gavTrajCoarse_mufix = cell(1,numsets_mufix);
    gubCoarse_mufix = cell(1,numsets_mufix);
    glbCoarse_mufix = cell(1,numsets_mufix);
    mavTrajCoarse_mufix = cell(1,numsets_mufix);
    mubCoarse_mufix = cell(1,numsets_mufix);
    mlbCoarse_mufix = cell(1,numsets_mufix);
    for jj = 1:numsets_mufix
        tstoreVec = tstoreVec_cell{setInds_mufix(jj)};
        [~,iA] = unique(floor(tstoreVec./deltat));
        tstoreNew = [0,tstoreVec(iA(2:end))];
        tCoarse_mufix{jj} = tstoreNew;
        gavTraj = gavTraj_cell{setInds_mufix(jj)};
        gavTrajCoarse_mufix{jj} = [1,gavTraj(iA(2:end))];
        gub = gavTraj_ub_cell{setInds_mufix(jj)};
        glb = gavTraj_lb_cell{setInds_mufix(jj)};
        gubCoarse_mufix{jj} = [1,gub(iA(2:end))];
        glbCoarse_mufix{jj} = [1,glb(iA(2:end))];
        mavTraj = mavTraj_cell{setInds_mufix(jj)};
        mavTrajCoarse_mufix{jj} = [1,mavTraj(iA(2:end))];
        mub = mavTraj_ub_cell{setInds_mufix(jj)};
        mlb = mavTraj_lb_cell{setInds_mufix(jj)};
        mubCoarse_mufix{jj} = [1,mub(iA(2:end))];
        mlbCoarse_mufix{jj} = [1,mlb(iA(2:end))];
    end

    % obtain F gradient based on when trajectories first hit some target
    % fitness
    Fhit_scan = [1.01,1.05,1.1,1.2,1.3,1.4];
    initFgradMat_mufix = zeros(numsets_mufix,length(Fhit_scan));
    for jj = 1:numsets_mufix
        tstoreVec = tCoarse_mufix{jj};
        gavTraj = gavTrajCoarse_mufix{jj};
        for kk = 1:length(Fhit_scan)
            Fhit = Fhit_scan(kk);
            tIndx = find(gavTraj>=Fhit,1);
            initgrad = (gavTraj(tIndx)-1)/tstoreVec(tIndx);
            initFgradMat_mufix(jj,kk) = initgrad;
        end
    end
    
    % obtain S gradient based on when trajectories first hit some target
    % fitness
    Shit_scan = [2,5,10,15];
    initSgradMat_mufix = zeros(numsets_Nfix,length(Shit_scan));
    for jj = 1:numsets_mufix
        tstoreVec = tCoarse_mufix{jj};
        mavTraj = mavTrajCoarse_mufix{jj};
        for kk = 1:length(Shit_scan)
            Shit = Shit_scan(kk);
            tIndx = find(mavTraj>=Shit,1);
            initgrad = mavTraj(tIndx)/tstoreVec(tIndx);
            initSgradMat_mufix(jj,kk) = initgrad;
        end
    end
end

% plot scaled trajectories
if numsets_mufix > 1
    colormat_mufix = rand(numsets_mufix,3);
    namecell_mufix = cell(1,numsets_mufix);
    numrows = ceil(sqrt(length(Fhit_scan)));
    numcols = ceil(length(Fhit_scan)/numrows);
    figure;
    for kk = 1:length(Fhit_scan)
        subplot(numrows,numcols,kk);
        for jj = 1:numsets_mufix
            N = N_sorted(jj);
            tstoreVec = tCoarse_mufix{jj};
            initgrad = initFgradMat_mufix(jj,kk);
            plot(tstoreVec.*initgrad,gavTrajCoarse_mufix{jj},'-.',...
                'color',colormat_mufix(jj,:),'LineWidth',1.5);
            hold on
            namecell_mufix{jj} = strcat('N = ',num2str(N));
        end
        for jj = 1:numsets_mufix
            tstoreVec = tCoarse_mufix{jj};
            initgrad = initFgradMat_mufix(jj,kk);
            fill([tstoreVec.*initgrad, fliplr(tstoreVec.*initgrad)],...
                [gubCoarse_mufix{jj},fliplr(glbCoarse_mufix{jj})],...
                colormat_mufix(jj,:),'linestyle', 'none','FaceAlpha',0.4);
            hold on
        end
%         for jj = 1:numsets_mufix
%             setIndx = setInds_mufix(jj);
%             initgrad = initFgradMat_mufix(jj,kk);
%             plot(tvec_cell{setIndx}.*initgrad,Fvec_cell{setIndx},...
%                 'color',colormat_mufix(jj,:),'LineWidth',1.5);
%             hold on
%         end        
%         xlim([0 3]);
        xlabel('t scaled');
        ylabel('F');
        title(strcat('Fhit = ',num2str(Fhit_scan(kk))));
        if kk == 1
            legend(namecell_mufix,'location','best');
        end
        
    end
    
    numrows = ceil(sqrt(length(Shit_scan)));
    numcols = ceil(length(Shit_scan)/numrows);
    figure;
    for kk = 1:length(Shit_scan)
        subplot(numrows,numcols,kk);
        for jj = 1:numsets_mufix
            N = N_sorted(jj);
            tstoreVec = tCoarse_mufix{jj};
            initgrad = initSgradMat_mufix(jj,kk);
            plot(tstoreVec.*initgrad,mavTrajCoarse_mufix{jj},'-.',...
                'color',colormat_mufix(jj,:),'LineWidth',1.5);
            hold on
            namecell_mufix{jj} = strcat('N = ',num2str(N));
        end
        for jj = 1:numsets_mufix
            tstoreVec = tCoarse_mufix{jj};
            initgrad = initSgradMat_mufix(jj,kk);
            fill([tstoreVec.*initgrad, fliplr(tstoreVec.*initgrad)],...
                [mubCoarse_mufix{jj},fliplr(mlbCoarse_mufix{jj})],...
                colormat_mufix(jj,:),'linestyle', 'none','FaceAlpha',0.4);
            hold on
        end
        for jj = 1:numsets_mufix
            setIndx = setInds_mufix(jj);
            initgrad = initSgradMat_mufix(jj,kk);
            plot(tvec_cell{setIndx}.*initgrad,Svec_cell{setIndx},...
                'color',colormat_mufix(jj,:),'LineWidth',1.5);
            hold on
        end    
        xlim([0 100]);
        xlabel('t scaled');
        ylabel('S');
        title(strcat('Shit = ',num2str(Shit_scan(kk))));
        if kk == 1
            legend(namecell_mufix,'location','best');
        end
        
    end
end

%% Extract sets with a given N*mu

%% Plot all sets
% close all;
colormat = rand(numfiles,3);
namecell = cell(1,numfiles);
figure;
subplot(2,1,1);
for fileIndx = 1:numfiles
    N = N_allsets(fileIndx);
    mub0 = mub0_allsets(fileIndx);
    plot(tvec_cell{fileIndx}.*FgradInit_CI_allsets(fileIndx),Fvec_cell{fileIndx},...
        'color',colormat(fileIndx,:),'LineWidth',1.5);
%     plot(tvec_cell{fileIndx}.*(N*mub0),Fvec_cell{fileIndx},...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
    hold on
    namecell{fileIndx} = strcat('N = ',num2str(N),', mu = ',num2str(mub0));
end
for fileIndx = 1:numfiles
    N = N_allsets(fileIndx);
    mub0 = mub0_allsets(fileIndx);
    tstoreVec_CI = tstoreVec_cell{fileIndx};
    gTrajMat = gTrajMat_cell{fileIndx};
%     for trialindx = 1:numtrials_allsets(fileIndx)
%         plot(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx),...
%             gTrajMat(trialindx,:),'-.','color',colormat(fileIndx,:),...
%             'LineWidth',1.2);
%         hold on
%     end
%     plot(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx),gavTraj_cell{fileIndx},'-.',...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
%     plot(tstoreVec_CI.*(N*mub0),gavTraj_cell{fileIndx},'-.',...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
    hold on
    fill([tstoreVec_CI.*FgradInit_CI_allsets(fileIndx), ...
        fliplr(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx))],...
        [gavTraj_ub_cell{fileIndx},fliplr(gavTraj_lb_cell{fileIndx})],...
        colormat(fileIndx,:),'linestyle', 'none','FaceAlpha',0.4);
%     fill([tstoreVec_CI.*(N*mub0), fliplr(tstoreVec_CI.*(N*mub0))],...
%         [gavTraj_ub_cell{fileIndx},fliplr(gavTraj_lb_cell{fileIndx})],...
%         colormat(fileIndx,:),'linestyle', 'none','FaceAlpha',0.4);
    hold on
end
% xlim([0 1e5]);
xlabel('tscaled');
% xlabel('t*(N \mu)');
ylabel('F');
legend(namecell,'location','best');

subplot(2,1,2);
for fileIndx = 1:numfiles
    N = N_allsets(fileIndx);
    mub0 = mub0_allsets(fileIndx);
    plot(tvec_cell{fileIndx}.*FgradInit_CI_allsets(fileIndx),...
        Svec_cell{fileIndx}.*FgradInit_CI_allsets(fileIndx)./SgradInit_CI_allsets(fileIndx),...
        'color',colormat(fileIndx,:),'LineWidth',1.5);
%     plot(tvec_cell{fileIndx}.*(N*mub0),Svec_cell{fileIndx},...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
    hold on    
end
for fileIndx = 1:numfiles
    N = N_allsets(fileIndx);
    mub0 = mub0_allsets(fileIndx);
    tstoreVec_CI = tstoreVec_cell{fileIndx};
    mTrajMat = mTrajMat_cell{fileIndx};
%     for trialindx = 1:numtrials_allsets(fileIndx)
%         plot(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx),...
%             mTrajMat(trialindx,:).*FgradInit_CI_allsets(fileIndx)./...
%             SgradInit_CI_allsets(fileIndx),'-.','color',colormat(fileIndx,:),...
%             'LineWidth',1.2);
%         hold on
%     end
%     plot(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx),...
%         mavTraj_cell{fileIndx}.*FgradInit_CI_allsets(fileIndx)./...
%         SgradInit_CI_allsets(fileIndx),'-.',...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
%     plot(tstoreVec_CI.*(N*mub0),mavTraj_cell{fileIndx},'-.',...
%         'color',colormat(fileIndx,:),'LineWidth',1.5);
    hold on
    fill([tstoreVec_CI.*FgradInit_CI_allsets(fileIndx), ...
        fliplr(tstoreVec_CI.*FgradInit_CI_allsets(fileIndx))],...
        [mavTraj_ub_cell{fileIndx}.*FgradInit_CI_allsets(fileIndx)./...
        SgradInit_CI_allsets(fileIndx),fliplr(mavTraj_lb_cell{fileIndx}.*...
        FgradInit_CI_allsets(fileIndx)./SgradInit_CI_allsets(fileIndx))],...
        colormat(fileIndx,:),'linestyle', 'none','FaceAlpha',0.4);
%     fill([tstoreVec_CI.*(N*mub0), fliplr(tstoreVec_CI.*(N*mub0))],...
%         [mavTraj_ub_cell{fileIndx},fliplr(mavTraj_lb_cell{fileIndx})],...
%         colormat(fileIndx,:),'linestyle', 'none','FaceAlpha',0.4);
    hold on
end
% xlim([0 1e5]);
xlabel('tscaled');
% xlabel('t*(N \mu)');
ylabel('Sscaled');
% ylabel('S');
legend(namecell,'location','best');
    
%% save
versionIndx = 3;
fn2save = strcat('AllSimData_',whichlandscape,'_v',num2str(versionIndx));
save(fn2save)

