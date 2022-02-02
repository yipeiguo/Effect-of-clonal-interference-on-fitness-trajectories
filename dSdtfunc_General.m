% This function specifies the functional form of both dF/dt and dS/dt, 
% which we will then use odesolver to solve.
% The equation for dS/dt depends on F, so to solve for S(t), we need to 
% simultaneously solve for F(t).
% INPUTS:
% here y is a vector with y(1) = F and y(2) = S.

function dydt = dSdtfunc_General(t,y,modelparams)

% Retrieve parameters:
if isfield(modelparams,'ifCI')
    ifCI = modelparams.ifCI;
else
    ifCI = true;
end

Fcommonfactorfunc = modelparams.Fcommonfactorfunc;
Scommonfactorfunc = modelparams.Scommonfactorfunc;
if ifCI == true
    Fintegfunc_CI = modelparams.Fintegfunc_CI;
    CIintegral_F = integral(@(stilde) Fintegfunc_CI(y(1),stilde),...
        0,modelparams.soversmean_max);
    Sintegfunc_CI = modelparams.Sintegfunc_CI;
    CIintegral_S = integral(@(stilde) Sintegfunc_CI(y(1),stilde),...
        0,modelparams.soversmean_max);
    
    dydt = [Fcommonfactorfunc(y(1))*CIintegral_F; ...
            Scommonfactorfunc(y(1))*CIintegral_S];
else
    dydt = [Fcommonfactorfunc(y(1))*2; ...
            Scommonfactorfunc(y(1))];
end




end