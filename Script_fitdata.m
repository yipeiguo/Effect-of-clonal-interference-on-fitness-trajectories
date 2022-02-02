% This is a script to fit experimental (/simulated) data to numerical
% solution.

close all;
clear variables

whichlandscape = 'DR';
addpath('../../.');

%% load data
pathname = './forPlotting_AnalyzedData/';
% fn = strcat('AllSimData_',whichlandscape);
fn = strcat('AllSimData_',whichlandscape,'_v2');
load(strcat(pathname,fn),'muOI','N_sorted','tCoarse_mufix',...
    'gavTrajCoarse_mufix','gubCoarse_mufix','glbCoarse_mufix',...
    'Fhit_scan','initFgradMat_mufix');

%% combine data into structure
numsets_all = length(N_sorted);
numsets2fit = min(3,numsets_all);
sets2fit = sort(randperm(numsets2fit));
% sets2fit = 1;
fn_sets2fit = strcat('fittedsets_',num2str(sets2fit(1)));
if numsets2fit > 1
    for kk = 2:numsets2fit
        fn_sets2fit = strcat(fn_sets2fit,'_',num2str(sets2fit(kk)));
    end
end

Nvec = N_sorted(sets2fit);
tscan_cell = tCoarse_mufix(sets2fit);
gavTraj_cell = gavTrajCoarse_mufix(sets2fit);
exptdata.Nvec = Nvec;
exptdata.tscan_cell = tscan_cell;
exptdata.gavTraj_cell = gavTraj_cell;

gubCoarse_cell = gubCoarse_mufix(sets2fit);
glbCoarse_cell = gubCoarse_mufix(sets2fit);
initFgradMat_rel = initFgradMat_mufix(sets2fit,:);

%% further coarse-grain data if desired
deltat = 100; % sample at most once in each deltat time interval
for jj = 1:numsets2fit
    tscan = tscan_cell{jj};
    [~,iA] = unique(floor(tscan./deltat));
    tstoreNew = [0,tscan(iA(2:end))];
    tscan_cell{jj} = tstoreNew;
    gavTraj = gavTraj_cell{jj};
    gavTraj_cell{jj} = [1,gavTraj(iA(2:end))];
    gub = gubCoarse_cell{jj};
    glb = glbCoarse_cell{jj};
    gubCoarse_cell{jj} = [1,gub(iA(2:end))];
    glbCoarse_cell{jj} = [1,glb(iA(2:end))];
end

% update gradients
for jj = 1:numsets2fit
    tscan = tscan_cell{jj};
    gavTraj = gavTraj_cell{jj};
    for kk = 1:length(Fhit_scan)
        Fhit = Fhit_scan(kk);
        tIndx = find(gavTraj>=Fhit,1);
        initgrad = (gavTraj(tIndx)-1)/tscan(tIndx);
        initFgradMat_rel(jj,kk) = initgrad;
    end
end


%% filename to save
versionIndx = 1;
if exist('deltat','var')
    fn2save = strcat('Datafitting_',whichlandscape,'_dt',num2str(deltat),...
        '_',fn_sets2fit,'_v',num2str(versionIndx));
else
    fn2save = strcat('Datafitting_',whichlandscape,'_',fn_sets2fit,...
        '_v',num2str(versionIndx));
end

%% Define other parameters
modelparams.ifCI = true;
modelparams.k = 1; % 1 for moran, 2 for WF, 2*log(D)/(D-1) for dilution protocol
modelparams.soversmean_max = 1000;

% grad_estimate = mean(initFgradMat_mufix,2);

%% scan parameters
alpha0_scan = 20:10:140;
alphaExp_scan = 0:1:10;
mub_scan = 10.^(-8:1:-4);
mubExp_scan = 0:0.1:1;
% alpha0_scan = 20:10:100;
% alphaExp_scan = 0;
% mub_scan = 10.^(-5);
% mubExp_scan = 1;

LSEscan = ones(length(alpha0_scan),length(alphaExp_scan),...
    length(mub_scan),length(mubExp_scan)).*100;
tic
for alpha0Indx = 1:length(alpha0_scan)
    alpha0 = alpha0_scan(alpha0Indx);
    fprintf('alpha = %d \n',alpha0);
    for alphaExpIndx = 1:length(alphaExp_scan)
        alphaExp = alphaExp_scan(alphaExpIndx);
        fprintf('alpha exponent = %d \n',alphaExp);
        for mubIndx = 1:length(mub_scan)
            mub0 = mub_scan(mubIndx);
            for mubExpIndx = 1:length(mubExp_scan)
                mubExp = mubExp_scan(mubExpIndx);
                params_all = [alpha0;alphaExp;mub0;mubExp];
                LSE = ObjFunc(params_all,exptdata,modelparams);
                LSEscan(alpha0Indx,alphaExpIndx,mubIndx,mubExpIndx) = LSE;
            end
        end
    end
    toc
end
toc

%% find the best parameters 
LSEscan(LSEscan==0) = inf;
[minLSE,linIndx] = min(LSEscan,[],'all','omitnan','linear');
[I_alpha,I_alphaExp,I_mub0,I_mubExp] = ind2sub(size(LSEscan),linIndx);
alpha0_opt = alpha0_scan(I_alpha);
alphaExp_opt = alphaExp_scan(I_alphaExp);
mub0_opt = mub_scan(I_mub0);
mubExp_opt = mubExp_scan(I_mubExp);

%% Estimate parameters
% alpha0_guess = 20;
% alphaExp_guess = 5; %0;
% mubExp_guess = 0;
% mub0_guess = 1e-5;
% mub0_guess = grad_estimate(1).*(alpha0_guess^2)./(2.*N_sorted(1));

xparams_guess = [alpha0_opt;alphaExp_opt;mub0_opt;mubExp_opt];
% xparams_guess = [alpha0_guess;alphaExp_guess;mub0_guess;mubExp_guess];
% lbvec = [0; -1; 0; -1];
lbvec = [0; 0; 0; 0];
ubvec = [1000;20;10*mub0_opt;inf];

%% carry out fitting
tic
[xparams_opt, fval, exitflag] = fmincon(@(x) ObjFunc(x,exptdata,modelparams),...
    xparams_guess,[],[],[],[],lbvec,ubvec);
toc

% if final solution is worse than initial guess, revert to initial guess
if fval > minLSE
    xparams_opt = xparams_guess;
end

%% numerical trajectory 
tic
[~, tvec_cell_opt, Fvec_cell_opt] = ObjFunc(xparams_opt,exptdata,modelparams);
toc

%% obtain initial gradient using fitted landscape parameters
modelparams_opt = modelparams;
modelparams_opt.alphafunc = @(x) xparams_opt(1).*x.^xparams_opt(2);
modelparams_opt.mubfunc = @(x) xparams_opt(3).*exp(-xparams_opt(4).*(x-1));

FgradInit_opt = zeros(numsets2fit,1);
for kk = 1:numsets2fit
    N = Nvec(kk);
    A = N*log(N)*modelparams_opt.k;
    modelparams_opt.A = A;

    % function representing degree of clonal interference
    lambdafunc = @(x,stilde) (modelparams_opt.A*modelparams_opt.mubfunc(x)).*(1+1./stilde).*exp(-stilde);
    Fintegfunc_CI = @(x,stilde) stilde.^2.*exp(-stilde).*exp(-lambdafunc(x,stilde)); 
    CIintegral_Finit = integral(@(stilde) Fintegfunc_CI(1,stilde),0,1000);    
    FgradInit_opt(kk) = CIintegral_Finit*N*xparams_opt(3)/(xparams_opt(1)^2);
end

%% Plot data and fitted trajectory
close all;
colormat = rand(numsets2fit,3);
namecell = cell(1,numsets2fit);  
% FhitIndxOI = min(3,size(initFgradMat_rel,2));
figure;
for kk = 1:numsets2fit
    N = Nvec(kk);
    namecell{kk} = strcat('N = ',num2str(N));
    tscan_data = tscan_cell{kk};
    gavTraj_data = gavTraj_cell{kk};
    gub_data = gubCoarse_cell{kk};
    glb_data = glbCoarse_cell{kk};
    
%     tscale = initFgradMat_rel(kk,FhitIndxOI)/initFgradMat_rel(1,FhitIndxOI);
    tscale = FgradInit_opt(kk)/FgradInit_opt(1);
    errorbar(tscan_data.*tscale,gavTraj_data,gavTraj_data-glb_data,...
        gub_data-gavTraj_data,'color',colormat(kk,:),'marker','o',...
        'MarkerSize',4,'MarkerEdgeColor',colormat(kk,:),...
        'MarkerFaceColor',colormat(kk,:),'LineStyle','none');
    hold on    
end
for kk = 1:numsets2fit
    tvec_fit = tvec_cell_opt{kk};
    Fvec_fit = Fvec_cell_opt{kk};
%     tscale = initFgradMat_rel(kk,FhitIndxOI)/initFgradMat_rel(1,FhitIndxOI);
    tscale = FgradInit_opt(kk)/FgradInit_opt(1);
    plot(tvec_fit.*tscale,Fvec_fit,'-','color','k','LineWidth',1.5);
%     plot(tvec_fit,Fvec_fit,'-','color',colormat(kk,:),'LineWidth',1.5);
    hold on    
end
% xlabel('t');
xlabel('tscaled');
ylabel('F');
legend(namecell,'location','best');
title(strcat('\alpha_0  = ',num2str(xparams_opt(1)),...
    ', g = ',num2str(xparams_opt(2)),', \mu_{b0} = ',num2str(xparams_opt(3)),...
    ', \gamma = ',num2str(xparams_opt(4))));

        

%% save
save(fn2save);






