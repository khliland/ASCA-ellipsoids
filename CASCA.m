function object = CASCA(responses, design, interactions, param, center, rescale_unbalance)
%% ASCA with various parameterizations
% object = ASCAp(response, factor, interaction, param, center);
%
% INPUTS:
% - responses, (nobj x nresp) matrix of continuous response vectors
% - design,    (nobj x nfac)  matrix of discrete design vectors
%                            (non-dummy factor levels, without intercept)
% - interactions, cell(s) of integers coding for pair-wise interaction, e.g. {[1,2]}
% - param,        parameterization ('contrast', 'sumtozero' or 'treatment')
% - center,       column-wise centering of response
% - rescale_unbalance Rescale objects when encountering unbalance
%
% OOUTPUT:
% - object, structure containing main elements and results of CASCA
if nargin < 6
    rescale_unbalance = false;
end
if nargin < 5
    center = true;
end
if nargin < 4 || isempty(param)
    param = 'sumtozero';
end
if nargin < 3 || isempty(interactions)
    interactions = [];
    nint = 0;
end
if nargin < 2
    error('Too few arguments supplied')
end

% Initialization
tol = 10^-12; % Tolerance for retaining a component (tol < singular value)
[nobj, nfac] = size(design);
if iscell(interactions)
    nint = length(interactions);
elseif length(interactions) > 1
    nint = 1;
    interactions = {interactions};
end

% Centering
responses_raw = responses;
if center
    response_mean = mean(responses);
    responses = bsxfun(@minus, responses, response_mean);
end
SSQ_responses = sum(responses(:).^2);

% Main effects
X = ones(nobj,1);
facPlace = zeros(2,nfac);
SSQ_main = zeros(1,nfac); SSQ_interaction = zeros(1,nint);
j = 1;
for i=1:nfac
    D = dummy(design(:,i), nobj, param);
    X = [X, D]; %#ok<AGROW>
    facPlace(:,i) = [(j+1); (j+size(D,2))];
    j = j+size(D,2);
end

% Interactions
if nint > 0
    for i=1:nint
        D = [];
        for k=1:diff(facPlace(:,interactions{i}(2)))+1
            D = [D, bsxfun(@times,X(:,facPlace(1,interactions{i}(1)):facPlace(2,interactions{i}(1))), X(:,facPlace(1,interactions{i}(2))+k-1))]; %#ok<AGROW>
        end
        X = [X, D]; %#ok<AGROW>
        intWidth = (diff(facPlace(:,interactions{i}(1)))+1)*(diff(facPlace(:,interactions{i}(2)))+1);
        facPlace(:,i+nfac) = [(j+1); (j+intWidth)];
        j = j+intWidth;
    end
end

% Check for unbalance
[~,~,u] = unique(X,'rows');
m = max(u);
us = zeros(1,m);
for i=1:m
    us(1,i) = sum(u==i);
end
if ~all(us==us(1)) && rescale_unbalance % Rescale if unbalanced
    c = round(median(us)); % Typical size of groups
    warning(['Unbalance: Assuming number of replicates = ' num2str(c)])
    for i=1:m
        X(u==i,:) = X(u==i,:).*sqrt(c/us(i));
        responses(u==i,:) = responses(u==i,:).*sqrt(c/us(i));
    end
end

% ANOVA
beta = X\responses;

% Main effects
effects    = cell(nfac,1);
scores     = effects; loadings = effects;
projected  = effects; singular = effects;
explained  = effects;
ieffects   = cell(nint,1);
iscores    = ieffects; iloadings = ieffects;
isingular  = ieffects; iexplained = ieffects;
iprojected = ieffects;
residuals  = responses;
for i = 1:nfac
    effects{i} = X(:,facPlace(1,i):facPlace(2,i)) * beta(facPlace(1,i):facPlace(2,i),:);
    residuals  = residuals - effects{i};
    SSQ_main(i) = sum(effects{i}(:).^2);
end
% Interaction effects
if nint > 0
    for i = 1:nint
        ieffects{i} = X(:,facPlace(1,nfac+i):facPlace(2,nfac+i)) * beta(facPlace(1,nfac+i):facPlace(2,nfac+i),:);
        residuals  = residuals - ieffects{i};
        SSQ_interaction(i) = sum(ieffects{i}(:).^2);
    end
end
SSQ_residual = sum(residuals(:).^2);
SSQ = [SSQ_main SSQ_interaction SSQ_residual];
expl = SSQ./SSQ_responses*100;

% SCA
for i = 1:nfac % Main effects
    [scores{i}, loadings{i}, singular{i}, explained{i}] = sca(effects{i}, diff(facPlace(:,i))+1, tol);
    projected{i} = residuals*loadings{i}; % Residuals after all effects have been subtracted, projected onto the loadings
end
factors = struct('scores',{scores}, 'loadings',{loadings}, ...
    'projected', {projected}, 'singular', {singular}, ...
    'explained',{explained});

if nint > 0 % Interactions
    for i = 1:nint
        [iscores{i}, iloadings{i}, isingular{i}, iexplained{i}] = sca(ieffects{i}, diff(facPlace(:,i+nfac))+1, tol);
        iprojected{i} = residuals*iloadings{i}; % Residuals after all main effects have been subtracted, projected onto the loadings
    end
    ifactors = struct('scores',{iscores}, 'loadings',{iloadings}, ...
        'projected', {iprojected}, 'singular', {isingular}, 'explained',{iexplained});
else
    ifactors = [];
end

% Return values
object = struct('effects',{effects}, 'ieffects',{ieffects}, ...
    'factors',factors, 'ifactors',ifactors, 'residuals',residuals, ...
    'responses',responses_raw, 'design',design, ...
    'beta',beta, 'dummy',X, 'explained',expl, 'SSQ', SSQ, ...
    'facID', facPlace, 'interactions', {interactions});


%% Dummy-coding from (possibly messy factor)
function D = dummy(X, n, param)
[~,~,u] = unique(X);
D = zeros(n, max(u));
for i=1:n
    D(i,u(i)) = 1;
end
switch param
    case 'contrast'
        D = D - D(:,end);
        D = D(:,1:(end-1));
    case 'sumtozero'
        D = bsxfun(@minus,D,mean(D)); % FIXME: *n?
        D = D(:,1:(end-1));
    case 'treatment'
        D = D(:,2:end);
    case 'centeredtreatment'
        D = center(D(:,1:(end-1)));
    otherwise
        error('Unsupported parameterization. Use ''contrast'', ''sumtozero'' or ''treatment''')
end


%% SCA
function [scores, loadings, singular, explained] = sca(X, ncomp, tol)
[u,s,v]   = svds(bsxfun(@minus,X,mean(X)), ncomp);
ds        = diag(s);
ncomp     = sum(ds > tol);
scores    = u(:,1:ncomp)*s(1:ncomp,1:ncomp);
loadings  = v(:,1:ncomp);
singular  = ds(1:ncomp);
explained = ds(1:ncomp).^2./sum(ds.^2).*100;
