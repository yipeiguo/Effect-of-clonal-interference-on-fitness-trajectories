% This is a function that carries out evolution in a Moran process using
% Gillespie simulation.

% In this v2, we speed up the code by using u = rand*s rather than dividing
% all elements in the rate vector by the sum of rates

function [gavTraj,nummutavTraj,numgrpsTraj,allpopsCell,finalpop] = ...
    EvolveMoran_Gillispie_v2(initialpop,modelparams,otherparams,namestr)

% Retrieve model parameters
alphafunc = modelparams.alphafunc; % mean beneficial fitness effect
mubfunc = modelparams.mubfunc;
N = modelparams.N;

% Retrieve other parameters (e.g. for simulation, storage, etc.)
tstoreVec = otherparams.tstoreVec;
maxT = max(tstoreVec);
% if isfield(otherparams,'namestr')
%     namestr = otherparams.namestr;
% end

% storage arrays
gavTraj = zeros(1,length(tstoreVec));
nummutavTraj = zeros(1,length(tstoreVec));
numgrpsTraj = zeros(1,length(tstoreVec));
allpopsCell = cell(1,length(tstoreVec));

% Retrieve current state of population
gvec = initialpop.gvec;             % growth rates of the different subgroups
Nvec = initialpop.Nvec;            % initial groups sizes of the different groups
NumMutvec = initialpop.NumMutvec;
smat = initialpop.smat;
sEtimemat = initialpop.sEtimemat;   % time(day) at which mutations emerge
assert(N == sum(Nvec),'inconsistent total number of cells!');

% initialize variables
currNvec = Nvec;
currgvec = gvec;
currNumMutvec = NumMutvec;
currsmat = smat;
currsEtimemat = sEtimemat;
currpop = initialpop;

tcurr = 0;
currNumGrps = length(currNvec);
while tcurr < maxT
    
    % rate of acquiring mutations for each grp
    mubvec = mubfunc(currgvec);  
    rmultiplyVec = currNvec.*currgvec.*(1-mubvec);
    rmutVec = currNvec.*currgvec.*mubvec;
    
    Rmat = (currNvec'*[rmultiplyVec,rmutVec])./N; % slow! (from profiler)
    Rmat(1:currNumGrps+1:currNumGrps^2) = 0;
    rates_all = Rmat(:);
%     rates_tot = sum(rates_all); % a bit slow
%     probvec = rates_all./rates_tot; % a bit slow
%     pcumsum = cumsum(probvec); % slow! (from profiler)
    pcumsum = cumsum(rates_all);
    rates_tot = pcumsum(end);
    
    % draw random event
    u = rand*rates_tot;
    % u = rand;
    whichevent = find(pcumsum>u,1,'first');
%     whichevent = pfsample(rates_all);
    
    % Calculate waiting time 
    dt = -log(rand) / rates_tot;
    % check if it's time to store data
    tstoreInds = find(tstoreVec>tcurr & tstoreVec<=(tcurr+dt));
    if ~isempty(tstoreInds) % store data if relevant
        fprintf('current time =%d, number of grps=%d \n',tcurr,currNumGrps)
        gavTraj(tstoreInds) = (currNvec*currgvec')./N;
        nummutavTraj(tstoreInds) = (currNvec*currNumMutvec')./N;
        numgrpsTraj(tstoreInds) = currNumGrps;
        allpopsCell(tstoreInds) = {currpop};
%         for kk = 1:length(tstoreInds)
%             allpopsCell{tstoreInds(kk)} = currpop;
%         end
        
        assignin('base','currpop',currpop);
        if exist('namestr','var')
            save(namestr);
        end
    end
    
    
    if whichevent <= currNumGrps^2
        % there is no mutation, just replacement of cells (i.e. selection)
        toreduce = mod(whichevent,currNumGrps);
        if toreduce == 0
            toreduce = currNumGrps;
        end
        toincrease = ceil(whichevent/currNumGrps);
        %[toreduce,toincrease] = ind2sub([currNumGrps,currNumGrps],whichevent);
        currNvec(toreduce) = currNvec(toreduce) - 1;
        currNvec(toincrease) = currNvec(toincrease) + 1;
        assert(min(currNvec)>=0,'number of cells cannot be negative!'); 
        
        if currNvec(toreduce) == 0
            currNvec(toreduce) = [];
            currgvec(toreduce) = [];
            currNumMutvec(toreduce) = [];
            currsmat(toreduce,:) = [];
            currsEtimemat(toreduce,:) = [];
            if size(currsmat,2) > max(currNumMutvec)
                currsmat(:,max(currNumMutvec)+1:end) = [];
                currsEtimemat(:,max(currNumMutvec)+1:end) = [];
            end
        end
    else % mutation occurs during division        
        
        newIndx = whichevent - currNumGrps^2;
        
        % find which subgroup mutated
        toreduce = mod(newIndx,currNumGrps);
        if toreduce == 0
            toreduce = currNumGrps;
        end
        whichgrp = ceil(newIndx/currNumGrps);
        %[toreduce,whichgrp] = ind2sub([currNumGrps,currNumGrps],newIndx);
        
        % add new mutant to list of groups
        meansb = 1./alphafunc(currgvec(whichgrp));
        smut = exprnd(meansb);
        gmut = currgvec(whichgrp)*(1+smut);
        nummut_mut = currNumMutvec(whichgrp)+1;
        % mub_mut = mubfunc(gmut);
        
        % store mutant properties
        if currNvec(toreduce) > 1 % subgroup in which death occured still remains
            nextNvec = [currNvec,1];
            nextNvec(toreduce) = nextNvec(toreduce) - 1;
            nextgvec = [currgvec,gmut];
            nextNumMutvec = [currNumMutvec,nummut_mut];
            if nummut_mut > max(currNumMutvec)
                if max(currNumMutvec) > 0
                    nextsmat = [currsmat,zeros(currNumGrps,1);...
                        currsmat(whichgrp,:), smut];
                    nextsEtimemat = [currsEtimemat,zeros(currNumGrps,1);...
                        currsEtimemat(whichgrp,:), tcurr];
                else
                    nextsmat = [zeros(currNumGrps,1);smut];
                    nextsEtimemat = [zeros(currNumGrps,1);tcurr];
                end
                    
            else
%                 size(currsmat)
%                 max(currNumMutvec)
                nextsmat = [currsmat;...
                    currsmat(whichgrp,1:currNumMutvec(whichgrp)),smut,...
                    zeros(1,max(currNumMutvec)-nummut_mut)];
                nextsEtimemat = [currsEtimemat;...
                    currsEtimemat(whichgrp,1:currNumMutvec(whichgrp)),tcurr,...
                    zeros(1,max(currNumMutvec)-nummut_mut)];
            end
            % update variables
            currNvec = nextNvec;
            currgvec = nextgvec;
            currNumMutvec = nextNumMutvec;
            currsmat = nextsmat;
            currsEtimemat = nextsEtimemat;            
        
        else % subgroup in which death occured dies since it only had a single cell
            % replace properties with the recently removed group by the new mutant
            % currNvec(toreduce) remains at 1
            if nummut_mut > max(currNumMutvec)
                nextsmat = [currsmat,zeros(currNumGrps,1)];
                nextsmat(toreduce,:) = [currsmat(whichgrp,:),smut];
                nextsEtimemat = [currsEtimemat,zeros(currNumGrps,1)];
                nextsEtimemat(toreduce,:) = [currsEtimemat(whichgrp,:),tcurr];
                currsmat = nextsmat;
                currsEtimemat = nextsEtimemat;
            else
                currsmat(toreduce,:) = [currsmat(whichgrp,1:currNumMutvec(whichgrp)),...
                    smut,zeros(1,max(currNumMutvec)-nummut_mut)];
                currsEtimemat(toreduce,:) = [currsEtimemat(whichgrp,1:currNumMutvec(whichgrp)),...
                    tcurr,zeros(1,max(currNumMutvec)-nummut_mut)];
            end
            currgvec(toreduce) = gmut;
            currNumMutvec(toreduce) = nummut_mut;
            if max(currNumMutvec) < size(currsmat,2)
                currsmat(:,max(currNumMutvec)+1:end) = [];
                currsEtimemat(:,max(currNumMutvec)+1:end) = [];
            end
            
        end
        
    end
    % update time
    tcurr = tcurr + dt;
    
    currNumGrps = length(currNvec);
    currpop.gvec = currgvec;
    currpop.Nvec = currNvec;
    currpop.smat = currsmat;
    currpop.sEtimemat = currsEtimemat;
    currpop.NumMutvec = currNumMutvec;    
    
end

finalpop = currpop;


end



% This is our own cumsum function, which is much slower than MATLAB's own
% cumsum even though this is much faster in julia!!
%{
function ii = pfsample(wvec)
    s = sum(wvec);
    t = rand() * s;
    ii = 1;
    cw = wvec(1);
    while cw < t && ii < length(wvec)
        ii = ii + 1;
        cw = cw + wvec(ii);
    end
end
%}




