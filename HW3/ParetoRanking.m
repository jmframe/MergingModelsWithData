function [NmbOfFront] = ParetoRanking(ObjVals);
% This function performs Pareto ranking

% First compute size of individual arrays
[nmbOfIndivs nmbOfObjs] = size(ObjVals);
% Pareto-optimal fronts
Front = {[]};
% number of Pareto-optimal front for each individual; 2nd highest priority sorting key
NmbOfFront = zeros( nmbOfIndivs, 1);
% set of individuals a particular individual dominates
Dominated = cell( nmbOfIndivs, 1);
% number of individuals by which a particular individual is dominated
NmbOfDominating = zeros( nmbOfIndivs, 1);

for p = 1:nmbOfIndivs,
    % First replicate current point
    Ptemp = ObjVals(p,1:nmbOfObjs); Pobj = Ptemp(ones(nmbOfIndivs,1),:);
    % Then find set of Dominated points
    [idx] = find((sum(Pobj <= ObjVals,[2]) == nmbOfObjs) & (sum(Pobj < ObjVals,[2]) > 0)); Nidx = length(idx);
    if Nidx > 0,
        Dominated{ p} = idx';
    end
    % Now find set of Nondominated points
    [idx] = find((sum(ObjVals <= Pobj,[2]) == nmbOfObjs) & (sum(ObjVals < Pobj,[2]) > 0)); Nidx = length(idx);
    if Nidx > 0,
        NmbOfDominating( p) = NmbOfDominating( p) + Nidx;
    end;

    if NmbOfDominating( p) == 0
        NmbOfFront( p) = 1;
        Front{ 1}(end + 1) = p;
    end
end
i = 1;

while ~isempty( Front{ i})
    NextFront = [];
    for k = 1:length( Front{ i})
        p = Front{ i}( k);
        % Instead of loop try direct implementation
        q = Dominated{p}; NmbOfDominating( q) = NmbOfDominating( q) - 1;
        % Now do second line
        idx = (NmbOfDominating( q)==0); NmbOfFront(q(idx)) = i + 1;
        % Updat the front
        NextFront = [NextFront q(idx)];
    end
    i = i + 1;
    Front{ end + 1} = NextFront;
end