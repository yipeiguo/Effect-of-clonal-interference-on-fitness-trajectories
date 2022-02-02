% This is a script to plot simulation data from 'EvolveMoran_Gillespie and
% 'Script_MoranEvolution_forcluster.m'

% In this v2, we compare simulation results to numerics obtained from
% odesolver.

close all; 
clear variables;
addpath('../');
addpath('../../');

whichlandscape = 'HoC'; % 'DR' or 'HoC' or 'RM'

Ssolve = true; % if true, we solve for both S and F. otherwise, find F only

timeunit = 'real'; % either 'generation' or 'real'

%% load data (CI regime)
if strcmp(whichlandscape,'HoC')
    runInds = 1:1:10;
    numruns = length(runInds);
    numrepsVec = ones(1,numruns); % vector of length numruns indicating the number of replicates in each run
    totnumtrials = sum(numrepsVec);

    pathname = './fromCluster/HoClandscape_N1e6_mu1eminus5/';
    fn = 'MoranEvolSimv2Data_N1000000_mub0_1e-05_alpha0_100_maxT_10000_landscape_HoC_numreplicates1_run';

elseif strcmp(whichlandscape,'DR')
    runInds = 1:1:10;
    numruns = length(runInds);
    numrepsVec = ones(1,numruns);
    totnumtrials = sum(numrepsVec);

    pathname = './fromCluster/DRlandscape_N1e6_mu1eminus5/';
    fn = 'MoranEvolSimv2Data_N1000000_mub0_1e-05_alpha0_20_maxT_10000_landscape_DR_numreplicates1_run';

elseif strcmp(whichlandscape,'RM')
   runInds = 1:1:10;
   numruns = length(runInds);
   numrepsVec = ones(1,numruns);
   totnumtrials = sum(numrepsVec);

   pathname = './fromCluster/RMlandscape_N1e6_mu1eminus5_set3/';
   fn = 'MoranEvolSimv2Data_N1000000_mub0_1e-05_alpha0_50_maxT_3000_landscape_RM_numreplicates1_run';

end

load(strcat(pathname,fn,num2str(runInds(1))),'tstoreVec','N','mub0','alpha0',...
    'modelparams');
if ~exist('mub0','var')
    mub0 = modelparams.mubfunc(1);
    alpha0 = modelparams.alphafunc(1);
end
tstoreVec_CI = tstoreVec;
N_CI = N;
mub0_CI = mub0;
alpha0_CI = alpha0;
modelparams_CI = modelparams;

gavtrajmat_CI = ones(totnumtrials,length(tstoreVec));
numgrpsmat_CI = ones(totnumtrials,length(tstoreVec));
nummutavtrajmat_CI = zeros(totnumtrials,length(tstoreVec));
numfixtrajmat_CI = zeros(totnumtrials,length(tstoreVec));
popTrajsCell_CI = cell(1,totnumtrials);
% finalPopCells = cell(1,totnumtrials);

counter = 1;
for runkk = 1:numruns
    runindx = runInds(runkk);
    numreps = numrepsVec(runkk);
    load(strcat(pathname,fn,num2str(runindx)), ...
        'gavtrajmat','nummutavtrajmat','numgrpsmat','popTrajsCell');
    for repindx = 1:numreps
        gavTraj = gavtrajmat(repindx,:);
        numgrpsTraj = numgrpsmat(repindx,:);
        nummutavTraj = nummutavtrajmat(repindx,:);
        allpopsCell = popTrajsCell{repindx};
        gavtrajmat_CI(counter,:) = gavTraj;
        numgrpsmat_CI(counter,:) = numgrpsTraj;
        nummutavtrajmat_CI(counter,:) = nummutavTraj;
        popTrajsCell_CI{counter} = allpopsCell;
%         finalPopCells{counter} = finalpop;
        % Extract number of fixed mutations at every time point:
        numfixTraj = zeros(1,length(tstoreVec));
        for tindx = 1:length(tstoreVec)
            if ~isempty(allpopsCell{tindx})
                whichgrp = find(allpopsCell{tindx}.NumMutvec == min(allpopsCell{tindx}.NumMutvec),1);
                numfixmuts = 0;
                for mutindx = 1:allpopsCell{tindx}.NumMutvec(whichgrp)
                    if sum(allpopsCell{tindx}.smat(:) == allpopsCell{tindx}.smat(whichgrp,mutindx)) == numgrpsTraj(tindx)
                        numfixmuts = numfixmuts + 1;
                    end
                end
                numfixTraj(tindx) = numfixmuts;
            end
        end
        numfixtrajmat_CI(counter,:) = numfixTraj;
        counter = counter + 1;
    end
end

% calculate average:
gavtraj_CI_mean = mean(gavtrajmat_CI,1);
if Ssolve == true
    nummutavtraj_CI_mean = mean(nummutavtrajmat_CI,1);
    numfixtraj_CI_mean = mean(numfixtrajmat_CI,1);
end

% unbiased estimate s of population standard deviation
gavtraj_CI_sd = std(gavtrajmat_CI,[],1);
if Ssolve == true
    nummutavtraj_CI_sd = std(nummutavtrajmat_CI,[],1);
    numfixtraj_CI_sd = std(numfixtrajmat_CI,[],1);
end

% calculate standard error of the mean (i.e. s/sqrt(N))
gavtraj_CI_sem = gavtraj_CI_sd./sqrt(totnumtrials);
if Ssolve == true
    nummutavtraj_CI_sem = nummutavtraj_CI_sd./sqrt(totnumtrials);
    numfixtraj_CI_sem = numfixtraj_CI_sd./sqrt(totnumtrials);
end

% 95percent confidence interval of the mean
gavtraj_CI_95CI = 1.96.*gavtraj_CI_sem;
nummutavtraj_CI_95CI = 1.96.*nummutavtraj_CI_sem;
numfixtraj_CI_95CI = 1.96.*numfixtraj_CI_sem;
    

gavtraj_CI_ub = gavtraj_CI_mean + gavtraj_CI_95CI;
gavtraj_CI_lb = gavtraj_CI_mean - gavtraj_CI_95CI;
% gavtraj_CI_ub = gavtraj_CI_mean + gavtraj_CI_sd;
% gavtraj_CI_lb = gavtraj_CI_mean - gavtraj_CI_sd;

nummutavtraj_CI_ub = nummutavtraj_CI_mean + nummutavtraj_CI_95CI;
nummutavtraj_CI_lb = nummutavtraj_CI_mean - nummutavtraj_CI_95CI;

numfixtraj_CI_ub = numfixtraj_CI_mean + numfixtraj_CI_95CI;
numfixtraj_CI_lb = numfixtraj_CI_mean - numfixtraj_CI_95CI;


%% numerical results using odesolver (with integral expression)
% solve ode here
modelparams.soversmean_max = 1000;
Tmax = max(tstoreVec_CI);
Finit = 1;

if strcmp(timeunit,'real')
    Fcommonfactorfunc = @(x) N*x^2*modelparams_CI.mubfunc(x)/modelparams_CI.alphafunc(x)^2;
elseif strcmp(timeunit,'generation')
    Fcommonfactorfunc = @(x) N*x*modelparams_CI.mubfunc(x)/modelparams_CI.alphafunc(x)^2;
end
% function representing degree of clonal interference
lambdafunc = @(x,stilde) (modelparams_CI.A*modelparams_CI.mubfunc(x)).*(1+1./stilde).*exp(-stilde);
Fintegfunc_CI = @(x,stilde) stilde.^2.*exp(-stilde).*exp(-lambdafunc(x,stilde)); 

modelparams.Fcommonfactorfunc = Fcommonfactorfunc;
modelparams.Fintegfunc_CI = Fintegfunc_CI;

if Ssolve == true
    Sinit = 0;
    if strcmp(timeunit,'real')
        Scommonfactorfunc = @(x) N*x*modelparams_CI.mubfunc(x)/modelparams_CI.alphafunc(x);
    elseif strcmp(timeunit,'generation')
        Scommonfactorfunc = @(x) N*modelparams_CI.mubfunc(x)/modelparams_CI.alphafunc(x);
    end
    Sintegfunc_CI = @(x,stilde) stilde.*exp(-stilde).*exp(-lambdafunc(x,stilde)); 
    modelparams.Scommonfactorfunc = Scommonfactorfunc;
    modelparams.Sintegfunc_CI = Sintegfunc_CI;

    modelparams.ifCI = true;
    [tvec_CI,FSvec_CI] = ode45(@(t,y) dSdtfunc_General(t,y,modelparams), [0 Tmax], [Finit;Sinit]);
    Fvec_CI = FSvec_CI(:,1);
    Svec_CI = FSvec_CI(:,2);
    modelparams.ifCI = false;
    [tvec_noCI,FSvec_noCI] = ode45(@(t,y) dSdtfunc_General(t,y,modelparams), [0 Tmax], [Finit;Sinit]);
    Fvec_noCI = FSvec_noCI(:,1);
    Svec_noCI = FSvec_noCI(:,2);
else
    modelparams.ifCI = true;
    [tvec_CI,Fvec_CI] = ode45(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], Finit);
%     [tvec_CI,Fvec_CI] = ode15s(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], Finit);
    modelparams.ifCI = false;
%     [tvec_noCI,Fvec_noCI] = ode15s(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], Finit);
    [tvec_noCI,Fvec_noCI] = ode45(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], Finit);
end



%% scale time if desired

if strcmp(timeunit,'generation')
    if fullsimCI_ifexist == true
        dtvec = tstoreVec_CI(2:end)-tstoreVec_CI(1:end-1);
        dtvec = gavtraj_CI_mean(1:end-1).*dtvec;
        tstoreVec_CI(2:end) = tstoreVec_CI(1:end-1) + dtvec;
    end
    if fullsimnoCI_ifexist == true
        dtvec = tstoreVec_noCI(2:end)-tstoreVec_noCI(1:end-1);
        dtvec = gavtraj_noCI_mean(1:end-1).*dtvec;
        tstoreVec_noCI(2:end) = tstoreVec_noCI(1:end-1) + dtvec;
    end
    tvec_noCI(2:end) = tvec_noCI(1:end-1) + ...
        Fvec_noCI(1:end-1).*(tvec_noCI(2:end)-tvec_noCI(1:end-1));
    tvec_CI(2:end) = tvec_CI(1:end-1) + ...
        Fvec_CI(1:end-1).*(tvec_CI(2:end)-tvec_CI(1:end-1));
end

%% Plot
if ~exist('mub0_noCI','var')
    mub0_noCI = 0.1/N_CI;
end
% tscale = 5e-5/1e-7;
tscale = mub0_CI/mub0_noCI;
CIintegral_Finit = integral(@(stilde) Fintegfunc_CI(Finit,stilde),...
    0,1000);
if Ssolve == true
    CIintegral_Sinit = integral(@(stilde) Sintegfunc_CI(Finit,stilde),...
        0,modelparams.soversmean_max);
end

% close all
figure;
plot(tstoreVec_CI,gavtraj_CI_mean,'r-.','LineWidth',1.5,...
    'DisplayName','full simulation (with CI)');
hold on
for trialindx = 1:totnumtrials
    plot(tstoreVec_CI,gavtrajmat_CI(trialindx,:),'r-.',...
        'DisplayName','sim data');
    hold on
end
fill([tstoreVec_CI, fliplr(tstoreVec_CI)],...
    [gavtraj_CI_ub,fliplr(gavtraj_CI_lb)],...
    'r','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data)');
hold on
% scaled data:
plot(tstoreVec_CI.*(CIintegral_Finit/2),gavtraj_CI_mean,'g-.','LineWidth',1.5,...
    'DisplayName','full simulation (with CI, scaled)');
hold on
fill([tstoreVec_CI.*(CIintegral_Finit/2), fliplr(tstoreVec_CI.*(CIintegral_Finit/2))],...
    [gavtraj_CI_ub,fliplr(gavtraj_CI_lb)],...
    'g','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data, scaled)');
hold on
plot(tvec_CI,Fvec_CI,'r-','DisplayName','with CI','LineWidth',1.5,...
    'DisplayName','Analytics (with CI)');
hold on
plot(tvec_noCI,Fvec_noCI,'b-','DisplayName','no CI','LineWidth',1.5,...
    'DisplayName','Analytics (no CI)');
hold on
% plot(tvec_lowmu_ode./tscale,Fvec_lowmu_ode,'b-','DisplayName','no CI','LineWidth',1.5);
% hold on
plot(tvec_CI.*(CIintegral_Finit/2),Fvec_CI,'g-',...
    'DisplayName','Analytics (with CI, scaled)','LineWidth',1.5);
hold on

% if strcmp(whichlandscape,'HoC')
% for trialindx = 1:numtrials_SSWM
%     plot(tvec_SSWM./tscale,gTrajMat_SSWM(trialindx,:),'b-.');
%     hold on
% end
% end
xlabel('t');
ylabel('<g>');
legend('location','best');
% legend('with CI','without CI');
title(strcat(whichlandscape, ' landscape'));
xlim([0 Tmax]);      

%%
if Ssolve == true
    figure;
    
    plot(tstoreVec_CI,nummutavtraj_CI_mean,'r-.','LineWidth',1.5,...
        'DisplayName','full simulation (with CI)');
    hold on
%         for trialindx = 1:totnumtrials
%             plot(tstoreVec_CI,nummutavtrajmat_CI(trialindx,:),'r-.',...
%                 'DisplayName','sim data');
%             hold on
%         end
    fill([tstoreVec_CI, fliplr(tstoreVec_CI)],...
        [nummutavtraj_CI_ub,fliplr(nummutavtraj_CI_lb)],...
        'r','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data)');
    hold on
    % scaled data:
    plot(tstoreVec_CI.*(CIintegral_Finit/2),...
        nummutavtraj_CI_mean.*(CIintegral_Finit/(2*CIintegral_Sinit)),...
        'g-.','LineWidth',1.5,...
        'DisplayName','full simulation (with CI, scaled)');
    hold on
    fill([tstoreVec_CI.*(CIintegral_Finit/2), fliplr(tstoreVec_CI.*(CIintegral_Finit/2))],...
        [nummutavtraj_CI_ub.*(CIintegral_Finit/(2*CIintegral_Sinit)),...
        fliplr(nummutavtraj_CI_lb.*(CIintegral_Finit/(2*CIintegral_Sinit)))],...
        'g','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data, scaled)');
    hold on

    
    plot(tvec_CI,Svec_CI,'r-','DisplayName','with CI','LineWidth',1.5,...
        'DisplayName','Analytics (with CI)');
    hold on
    plot(tvec_noCI,Svec_noCI,'b-','DisplayName','no CI','LineWidth',1.5,...
        'DisplayName','Analytics (no CI)');
    hold on
    plot(tvec_CI.*(CIintegral_Finit/2),...
            Svec_CI.*(CIintegral_Finit/(2*CIintegral_Sinit)),'g-',...
            'DisplayName','Analytics (with CI, scaled)','LineWidth',1.5);
    hold on
    % if strcmp(whichlandscape,'HoC')
    % for trialindx = 1:numtrials_SSWM
    %     plot(tvec_SSWM./tscale,gTrajMat_SSWM(trialindx,:),'b-.');
    %     hold on
    % end
    % end
    xlabel('t');
    ylabel('<m>');
    legend('location','best');
    title(strcat(whichlandscape, ' landscape'));
    % xlim([0 1e3]);      
    xlim([0 Tmax]);  
end

if Ssolve == true
    figure;
    
    plot(tstoreVec_CI,numfixtraj_CI_mean,'r-.','LineWidth',1.5,...
        'DisplayName','full simulation (with CI)');
    hold on
%         for trialindx = 1:totnumtrials
%             plot(tstoreVec_CI,numfixtrajmat_CI(trialindx,:),'r-.',...
%                 'DisplayName','sim data');
%             hold on
%         end
    fill([tstoreVec_CI, fliplr(tstoreVec_CI)],...
        [numfixtraj_CI_ub,fliplr(numfixtraj_CI_lb)],...
        'r','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data)');
    hold on        
    % scaled data:
    plot(tstoreVec_CI.*(CIintegral_Finit/2),...
        numfixtraj_CI_mean.*(CIintegral_Finit/(2*CIintegral_Sinit)),...
        'g-.','LineWidth',1.5,...
        'DisplayName','full simulation (with CI, scaled)');
    hold on
    fill([tstoreVec_CI.*(CIintegral_Finit/2), fliplr(tstoreVec_CI.*(CIintegral_Finit/2))],...
        [numfixtraj_CI_ub.*(CIintegral_Finit/(2*CIintegral_Sinit)),...
        fliplr(numfixtraj_CI_lb.*(CIintegral_Finit/(2*CIintegral_Sinit)))],...
        'g','linestyle', 'none','FaceAlpha',0.4,'DisplayName','error (sim data, scaled)');
    hold on

    plot(tvec_CI,Svec_CI,'r-','DisplayName','with CI','LineWidth',1.5,...
        'DisplayName','Analytics (with CI)');
    hold on
    plot(tvec_noCI,Svec_noCI,'b-','DisplayName','no CI','LineWidth',1.5,...
        'DisplayName','Analytics (no CI)');
    hold on
    plot(tvec_CI.*(CIintegral_Finit/2),...
            Svec_CI.*(CIintegral_Finit/(2*CIintegral_Sinit)),'g-',...
            'DisplayName','Analytics (with CI, scaled)','LineWidth',1.5);
    hold on
    % if strcmp(whichlandscape,'HoC')
    % for trialindx = 1:numtrials_SSWM
    %     plot(tvec_SSWM./tscale,gTrajMat_SSWM(trialindx,:),'b-.');
    %     hold on
    % end
    % end
    xlabel('t');
    ylabel('# fixations');
    legend('location','best');
    title(strcat(whichlandscape, ' landscape'));
    % xlim([0 1e3]);      
    xlim([0 Tmax]);  
end

%% Save data (for plotting in figure)
% setIndx = 2;
% fn_tosave = strcat('Data_CIeffectonTrajs_N',num2str(N),'_mub0_',num2str(mub0),'_alpha0_',...
%     num2str(alpha0),'_maxT_',num2str(max(tstoreVec)),'_landscape_',whichlandscape,...
%     '_numtrials',num2str(totnumtrials),'_set',num2str(setIndx));
% fn_tosave = strrep(fn_tosave,'.','pt');
% save(fn_tosave);

% save(strcat('Data_CIeffectonTrajs_landscape',whichlandscape));
