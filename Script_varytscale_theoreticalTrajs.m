% Here we explore how scaling time using estimated gradients affect
% theoretical trajectories (obtained from solving odes)

close all;
clear variables

whichlandscape = 'RM';

%% load data
pathname = './forPlotting_AnalyzedData/';
fn = strcat('AllSimData_',whichlandscape);
% fn = strcat('AllSimData_',whichlandscape,'_v3');
load(strcat(pathname,fn),'N_allsets','mub0_allsets','maxT_allsets',...
    'tvec_cell','Fvec_cell','FgradInit_CI_allsets');

%% Extract sets with a given mu
muOI = 1e-5;
muIndsOI = find(mub0_allsets == muOI);
N_rel = N_allsets(muIndsOI);
[N_sorted,N_I] = sort(N_rel);
setInds = muIndsOI(N_I);
numsets = length(muIndsOI);
N_rel = N_sorted;
tvec_rel = tvec_cell(setInds);
Fvec_rel = Fvec_cell(setInds);
Fgrad_rel = FgradInit_CI_allsets(setInds);

namecell = cell(1,numsets);
for jj = 1:numsets
    N = N_rel(jj);
    namecell{jj} = strcat('N = ',num2str(N));
end


%% obtain empirical initial F gradient
% based on when trajectories first hit some target fitness
Fhit_scan = [1.01,1.05,1.1,1.2,1.3,1.4];
initFgradMat = zeros(numsets,length(Fhit_scan));
for jj = 1:numsets
    tvec = tvec_rel{jj};
    Fvec = Fvec_rel{jj};
    for kk = 1:length(Fhit_scan)
        Fhit = Fhit_scan(kk);
        tIndx = find(Fvec>=Fhit,1);
        initgrad = (Fvec(tIndx)-1)/tvec(tIndx);
        initFgradMat(jj,kk) = initgrad;
    end
end

%% plot how estimate of the gradient varies with Fhit
colormat = rand(numsets,3);
mksize = 10;
figure;
for jj = 1:numsets
    scatter(Fhit_scan,initFgradMat(jj,:)./initFgradMat(jj,1),mksize,colormat(jj,:),...
        'markerfacecolor',colormat(jj,:));
%     scatter(Fhit_scan,initFgradMat(jj,:)./Fgrad_rel(jj),mksize,colormat(jj,:),...
%         'markerfacecolor',colormat(jj,:));
    hold on
end
xlabel('F_c');
ylabel('slope');
legend(namecell,'location','best');
    
%% plot scaled trajectories
if numsets > 1
    numrows = ceil(sqrt(length(Fhit_scan)));
    numcols = ceil(length(Fhit_scan)/numrows);
    figure;
    for kk = 1:length(Fhit_scan)
        subplot(numrows,numcols,kk);
        for jj = 1:numsets
            initgrad = initFgradMat(jj,kk)./initFgradMat(1,kk);
            tvec = tvec_rel{jj};
            Fvec = Fvec_rel{jj};
            plot(tvec.*initgrad,Fvec,'color',colormat(jj,:),'LineWidth',1.5);
            hold on
        end        
%         xlim([0 3]);
        xlabel('t scaled');
        ylabel('F');
        title(strcat('Fhit = ',num2str(Fhit_scan(kk))));
        if kk == 1
            legend(namecell,'location','best');
        end        
    end
end

%% save
% fn2save = strcat('Data_varyFc_',whichlandscape,'_v1');
% save(fn2save);
    
    
    
    
    
    
    
