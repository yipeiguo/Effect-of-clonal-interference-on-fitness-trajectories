% This is a function that defines the objective function when we carry out
% fitting of simulation data to numerical solution of ode (to infer 
% parameters of underlying fitness landscape.

% Input parameters:
% - modelparams that should already have the following fields:
%     - modelparams.ifCI = true;
%     - modelparams.k % 1 for moran, 2 for WF, 2*log(D)/(D-1) for dilution protocol
%     - modelparams.soversmean_max

function [LSE, tvec_cell, Fvec_cell] = ObjFunc(xparams,exptdata,modelparams)

    % Extract experimental data
    Nvec = exptdata.Nvec;
    tCoarse_cell = exptdata.tscan_cell;
    gavTrajCoarse_cell = exptdata.gavTraj_cell;
    numsets = length(Nvec);

    % Extract landscape parameters to optimize over
    alpha0 = xparams(1);
    alpha_exponent = xparams(2);
    mub0 = xparams(3);
    mub_exponent = xparams(4);
    modelparams.alphafunc = @(x) alpha0.*x.^alpha_exponent;
    modelparams.mubfunc = @(x) mub0.*exp(-mub_exponent.*(x-1));

    % obtain fitness trajectories from numerically solving ode
    tvec_cell = cell(1,numsets);
    Fvec_cell = cell(1,numsets);
    numdatapts = 0;
    LSE = 0;
    for setIndx = 1:numsets
        N = Nvec(setIndx);
        A = N*log(N)*modelparams.k;
        modelparams.A = A;

        Fcommonfactorfunc = @(x) N*x^2*modelparams.mubfunc(x)/modelparams.alphafunc(x)^2;
        % function representing degree of clonal interference
        lambdafunc = @(x,stilde) (modelparams.A*modelparams.mubfunc(x)).*(1+1./stilde).*exp(-stilde);
        Fintegfunc_CI = @(x,stilde) stilde.^2.*exp(-stilde).*exp(-lambdafunc(x,stilde)); 
        
        modelparams.Fcommonfactorfunc = Fcommonfactorfunc;
        modelparams.Fintegfunc_CI = Fintegfunc_CI;
        
        % solve trajectory
        tCoarseVec = tCoarse_cell{setIndx}; 
        Tmax = max(tCoarseVec);
        gmax = max(gavTrajCoarse_cell{setIndx})*1e6;
        opts = odeset('Events',@(t,y) eventfun(t,y,gmax));
        [tvec,Fvec,te,~,~] = ode45(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], 1, opts);
%         [tvec,Fvec,te,Ffinal,ie] = ode45(@(t,y) dFdtfunc_General(t,y,modelparams), [0 Tmax], 1, opts);
        Fvec_cell{setIndx} = Fvec;
        tvec_cell{setIndx} = tvec;
        
        % Compare with data
        gCoarseVec = gavTrajCoarse_cell{setIndx};
        if ~isempty(te)
            tCoarseVec = tCoarseVec(tCoarseVec <= te);
            gCoarseVec = gCoarseVec(tCoarseVec <= te);
        end
        Fvec_coarse = interp1(tvec,Fvec,tCoarseVec,'linear','extrap');
        LSE = LSE + sum((gCoarseVec - Fvec_coarse).^2);
        numdatapts = numdatapts + length(tCoarseVec);
    end
    LSE = LSE/numdatapts;
    
end

