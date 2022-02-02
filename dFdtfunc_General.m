% This function specifies the functional form of dF/dt, which we will then
% use odesolver to solve.

function dydt = dFdtfunc_General(t,y,modelparams)

% Retrieve parameters:
if isfield(modelparams,'ifCI')
    ifCI = modelparams.ifCI;
else
    ifCI = true;
end
Fcommonfactorfunc = modelparams.Fcommonfactorfunc;
if ifCI == true
    Fintegfunc_CI = modelparams.Fintegfunc_CI;
    CIintegral_F = integral(@(stilde) Fintegfunc_CI(y,stilde),...
        0,modelparams.soversmean_max);
    dydt = Fcommonfactorfunc(y)*CIintegral_F;
else
    dydt = Fcommonfactorfunc(y)*2;
end



end