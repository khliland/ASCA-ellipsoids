function [area, mult] = plotCASCA(casca, factor, ellipsoids, dims, dim1, dim2, dim3, points)
%% Plotting of scores from ASCA
% plotCASCA(casca, 1,2) % Plot model ellipsoids for factor 1
%
% INPUTS
% casca      - ASCA object
% factor     - which factor to plot (negative for interaction number)
% ellipsoids - data ellipsoids (1) or model ellipsoids (2)
% dims       - number of dimensions
% dim1       - first dimension
% dim2       - second dimension
% dim3       - third dimension
% points     - plot points (true or false)
%
% OUTPUTS
% area - area of ellispoids
% mult - multiple testing results
%
if nargin < 2
    error('Too few arguments supplied')
end
if nargin < 3
    ellipsoids = 1;
end
if nargin < 4
    dims = 2;
end
if nargin < 5
    dim1 = 1;
end
if nargin < 6
    dim2 = 2;
end
if nargin < 7
    dim3 = 3;
end
if nargin < 8
    points = true;
end

% Initialization
[nobj,nfac] = size(casca.design);
if factor > nfac
    error('Selected factor not available')
end
if factor > 0 % Factor number
    if size(casca.factors.scores{factor},2) == 1 || dims == 1
        dim1 = 1; dim2 = 1; dim3 = 1;
        dims = 1;
    elseif size(casca.factors.scores{factor},2) == 2
        dim1 = 1; dim2 = 2; dim3 = 2;
        dims = min(dims,2);
    end
    % Extract factors from casca/ASCA object
    explained = casca.factors.explained{factor};
    loadings  = casca.factors.loadings{factor}(:,[dim1,dim2,dim3]);
    projected = casca.factors.scores{factor}(:,[dim1,dim2,dim3]) + ...
        casca.factors.projected{factor}(:,[dim1,dim2,dim3]);
    [~,levs,where] = unique(casca.design(:,factor));
else % Interaction
    factor = -factor;
    if size(casca.ifactors.scores{factor},2) == 1 || dims == 1
        dim1 = 1; dim2 = 1; dim3 = 1;
        dims = 1;
    elseif size(casca.ifactors.scores{factor},2) == 2
        dim1 = 1; dim2 = 2; dim3 = 2;
        dims = min(dims,2);
    end
    explained = casca.ifactors.explained{factor};
    loadings  = casca.ifactors.loadings{factor}(:,[dim1,dim2,dim3]);
    projected = casca.ifactors.scores{factor}(:,[dim1,dim2,dim3]) + ...
        casca.ifactors.projected{factor}(:,[dim1,dim2,dim3]);
    if isstring(casca.design)
        design = casca.design(:,casca.interactions{factor}(1));
    else
        design = num2str(casca.design(:,casca.interactions{factor}(1)));
    end
    for i = 2:length(casca.interactions{factor})
        if isstring(casca.design)
            design = [design casca.design(:,casca.interactions{factor}(i))]; %#ok<AGROW>
        else
            design = [design num2str(casca.design(:,casca.interactions{factor}(i)))]; %#ok<AGROW>
        end
    end
    [~,levs,where] = unique(design,'rows');
    factor = -factor;
end
% Plot symbols and colors
syms = {'o','*','x','d','s','^','+','.','v','p','h','<','>','o','*','x','d','s','^','+','.','v','p','h','<','>', ...
    'o','*','x','d','s','^','+','.','v','p','h','<','>','o','*','x','d','s','^','+','.','v','p','h','<','>', ...
    'o','*','x','d','s','^','+','.','v','p','h','<','>','o','*','x','d','s','^','+','.','v','p','h','<','>', ...
    'o','*','x','d','s','^','+','.','v','p','h','<','>','o','*','x','d','s','^','+','.','v','p','h','<','>'};
lins = {'-','--','-.',':'};
cols = get(groot,'DefaultAxesColorOrder'); cols = [cols;cols;(cols+ones(7,3).*2)./3;cols;cols;cols;cols;cols;cols;cols;cols;cols;cols;cols;cols;cols;cols];
area = zeros(length(levs),dims+1);
mult = [];
if ellipsoids == 2
    colsg = repmat(0.7,size(cols,1),3);
else
    colsg = cols;
end

% Plot scores
hold on
if dims == 1
    ylim([0,length(levs)+1])
    for i=1:length(levs)
        plot([min(projected(:,1)),max(projected(:,1))],[i,i],'-','Color',[0.4,0.4,0.4])
        if points
            plot(projected(where==i,1), i, syms{i},'Color',cols(i,:))
        end
    end
    xlabel('Observations')
    ylabel('Groups')
    title('ASCA data confidence')
elseif dims == 2
    if points
        for i=1:length(levs)
            plot(projected(where==i,1), projected(where==i,2), syms{i},'Color',colsg(i,:))
        end
    end
    xlabel(['Component ' num2str(dim1) ' - ' num2str(round(explained(dim1)*100)/100) '%'])
    ylabel(['Component ' num2str(dim2) ' - ' num2str(round(explained(dim2)*100)/100) '%'])
else
    if points
        for i=1:length(levs)
            plot3(projected(where==i,1), projected(where==i,2), projected(where==i,3), syms{i},'Color',colsg(i,:))
        end
    end
    for i=1:length(levs)
        plot3(mean(projected(where==i,1)), mean(projected(where==i,2)), mean(projected(where==i,3)), 'o','Color',cols(i,:),'MarkerSize',7,'MarkerFaceColor',[0,0,0],'LineWidth',1.5)
    end
    xlabel(['Component ' num2str(dim1) ' - ' num2str(round(explained(dim1)*100)/100) '%'])
    ylabel(['Component ' num2str(dim2) ' - ' num2str(round(explained(dim2)*100)/100) '%'])
    zlabel(['Component ' num2str(dim3) ' - ' num2str(round(explained(dim3)*100)/100) '%'])
end

% Plot data ellipsoids
if ellipsoids == 1 && dims == 2
    if factor < 0
        title(['Interaction ' num2str(-factor) ' (data ellipsoids)'])
    else
        title(['Factor ' num2str(factor) ' (data ellipsoids)'])
    end
    C = normrnd(0,1,[10,2]); C = bsxfun(@times,C,1./sqrt(sum(C.^2,2)));
    
    for i=1:length(levs)
        Sigma = cov(projected(where==i,1:2));
        S     = chol(Sigma,'upper');
        nobji = sum(where==i);

        s1 = bsxfun(@plus, (C*S).*sqrt(2*(nobji-1)/(nobji-2)*finv(0.40,2,nobji-2)), mean(projected(where==i,1:2)));
        s2 = bsxfun(@plus, (C*S).*sqrt(2*(nobji-1)/(nobji-2)*finv(0.68,2,nobji-2)), mean(projected(where==i,1:2)));
        s3 = bsxfun(@plus, (C*S).*sqrt(2*(nobji-1)/(nobji-2)*finv(0.95,2,nobji-2)), mean(projected(where==i,1:2)));
        ellipse_t1 = fit_ellipse( s1(:,1),s1(:,2));
        ellipse_t2 = fit_ellipse( s2(:,1),s2(:,2));
        ellipse_t3 = fit_ellipse( s3(:,1),s3(:,2));
        plot3(ellipse_t1.ellipse(1,:),ellipse_t1.ellipse(2,:),-ones(100,1),'-','Color',cols(i,:))
        plot3(ellipse_t2.ellipse(1,:),ellipse_t2.ellipse(2,:),-ones(100,1),'-','Color',cols(i,:))
        plot3(ellipse_t3.ellipse(1,:),ellipse_t3.ellipse(2,:),-ones(100,1),'-','Color',cols(i,:))

        area(i,1) = ellipse_t3.long_axis*ellipse_t3.short_axis*pi/4;
        area(i,2) = ellipse_t3.short_axis;
        area(i,3) = ellipse_t3.long_axis;
    end
end

% Plot confidence intervals
if ellipsoids == 2 && dims == 1
    UHat      = casca.residuals;
    SigmaHat1 = (UHat'*UHat)./nobj;
    d         = 1;
    if factor < 0
        q         = casca.facID(2,nfac-factor)-casca.facID(1,nfac-factor)+1;
    else
        q         = casca.facID(2,factor)-casca.facID(1,factor)+1;
    end
    scaling   = q/nobj;
    S         = sqrt(loadings(:,1)'*SigmaHat1*loadings(:,1)*scaling*nobj/(nobj-q));
    t1        = sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.40,d,nobj-q-d+1));
    t2        = sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.68,d,nobj-q-d+1));
    t3        = sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.95,d,nobj-q-d+1));
    for i=1:length(levs)
        plot3([1,1]*(mean(projected(where==i,1)) - S*t1), i+0.2*[-1,1],-[1,1],'--k')
        plot3([1,1]*(mean(projected(where==i,1)) - S*t2), i+0.2*[-1,1],-[1,1],':k')
        plot3([1,1]*(mean(projected(where==i,1)) - S*t3), i+0.2*[-1,1],-[1,1],'-.k')
        plot3([1,1]*(mean(projected(where==i,1)) + S*t1), i+0.2*[-1,1],-[1,1],'--k')
        plot3([1,1]*(mean(projected(where==i,1)) + S*t2), i+0.2*[-1,1],-[1,1],':k')
        plot3([1,1]*(mean(projected(where==i,1)) + S*t3), i+0.2*[-1,1],-[1,1],'-.k')
        area(i,1) = S*t3*2;
    end
end

% Plot model ellipsoids
if ellipsoids == 2 && dims == 2
    if factor < 0
        title(['Interaction ' num2str(-factor) ' (model ellipsoids)'])
    else
        title(['Factor ' num2str(factor) ' (model ellipsoids)'])
    end
    UHat      = casca.residuals;
    SigmaHat1 = (UHat'*UHat)./nobj;
    d         = 2;
    projected = projected(:,1:2);
    if factor < 0
        q         = casca.facID(2,nfac-factor)-casca.facID(1,nfac-factor)+1;
    else
        q         = casca.facID(2,factor)-casca.facID(1,factor)+1;
    end
    scaling   = q/nobj;
    ASA       = loadings(:,1:2)'*SigmaHat1*loadings(:,1:2).*(scaling)./(nobj-q).*(nobj);%.*((nobj/2+1)./(nobj));%.*sqrt((k-1)./nobj);%.*(nobj)./(nobj-k);
    ASA       = (ASA+ASA').*0.5; % Force symmetry
    S         = chol(ASA,'upper');
    C = normrnd(0,1,[10,2]); C = bsxfun(@times,C,1./sqrt(sum(C.^2,2)));
    if nobj > q+1
        mult = multComp(projected, where, ASA, nobj, q, d, length(levs));
    else
        mult = [];
    end
    
    for i=1:length(levs)
        s1 = bsxfun(@plus, (C*S).*sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.40,d,nobj-q-d+1)), mean(projected(where==i,1:2)));
        s2 = bsxfun(@plus, (C*S).*sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.68,d,nobj-q-d+1)), mean(projected(where==i,1:2)));
        s3 = bsxfun(@plus, (C*S).*sqrt((nobj-q)*d./(nobj-q-d+1).*finv(0.95,d,nobj-q-d+1)), mean(projected(where==i,1:2)));
        ellipse_t1 = fit_ellipse( s1(:,1), s1(:,2) );
        ellipse_t2 = fit_ellipse( s2(:,1), s2(:,2) );
        ellipse_t3 = fit_ellipse( s3(:,1), s3(:,2) );
        plot3(ellipse_t1.ellipse(1,:), ellipse_t1.ellipse(2,:), ones(100,1),lins{floor((i-1)/7)+1},'Color',cols(i,:))
        plot3(ellipse_t2.ellipse(1,:), ellipse_t2.ellipse(2,:), ones(100,1),lins{floor((i-1)/7)+1}, 'Color',cols(i,:))
        plot3(ellipse_t3.ellipse(1,:), ellipse_t3.ellipse(2,:), ones(100,1),lins{floor((i-1)/7)+1},'Color',cols(i,:))

        area(i,1) = ellipse_t3.long_axis*ellipse_t3.short_axis*pi/4;
        area(i,2) = ellipse_t3.short_axis;
        area(i,3) = ellipse_t3.long_axis;
    end
end

% 3D data ellipsoids
if ellipsoids == 1 && dims == 3
    if factor < 0
        title(['Interaction ' num2str(-factor) ' (data ellipsoids)'])
    else
        title(['Factor ' num2str(factor) ' (data ellipsoids)'])
    end
    
    colormap(cols)
    for i=1:length(levs)
        Sigma = cov(projected(where==i,:));
        
        nobji = sum(where==i);
        C = Sigma.*(3*(nobji-1)/(nobji-3)*finv(0.95,3,nobji-3));
        M = mean(projected(where==i,:));
        [U,L] = eig(C);
        
        % For N standard deviations spread of data, the radii of the eliipsoid will
        % be given by N*SQRT(eigenvalues).
        
        N = 1; % choose your own N
        radii = N*sqrt(diag(L));
        
        % Generate data for "unrotated" ellipsoid
        [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
        
        % Rotate data with orientation matrix U and center M
        a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
        data = a+b+c; n = size(data,2);
        x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
        
        % Now plot the rotated ellipse
        surf(x,y,z,'FaceColor',cols(i,:),'FaceAlpha',0.1,'EdgeAlpha',0.5);
        area(i,1) = prod(radii)*4/3*pi;
        area(i,2:end) = radii;
    end
end

% Plot model ellipsoids in 3D
if ellipsoids == 2 && dims == 3
    if factor < 0
        title(['Interaction ' num2str(-factor) ' (model ellipsoids)'])
    else
        title(['Factor ' num2str(factor) ' (model ellipsoids)'])
    end
    UHat      = casca.residuals;
    SigmaHat1 = (UHat'*UHat)./nobj;
    if factor < 0
        q         = casca.facID(2,nfac-factor)-casca.facID(1,nfac-factor)+1;
    else
        q         = casca.facID(2,factor)-casca.facID(1,factor)+1;
    end
    scaling = q/nobj;
    d         = 3;
    ASA       = loadings(:,1:3)'*SigmaHat1*loadings(:,1:3).*(scaling)./(nobj-q).*(nobj);%.*((nobj/2+1)./(nobj));%.*sqrt((k-1)./nobj);%.*(nobj)./(nobj-k);
    
    for i=1:length(levs)
        Sigma = ASA;
        C = Sigma.*((nobj-q)*d./(nobj-q-d+1).*finv(0.95,d,nobj-q-d+1));
        M = mean(projected(where==i,:));
        [U,L] = eig(C);
        
        % For N standard deviations spread of data, the radii of the eliipsoid will
        % be given by N*SQRT(eigenvalues).
        
        N = 1; % Choose your own N
        radii = N*sqrt(diag(L));
        
        % Generate data for "unrotated" ellipsoid
        [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
        
        % Rotate data with orientation matrix U and center M
        a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
        data = a+b+c; n = size(data,2);
        x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
        
        % Now plot the rotated ellipse
        surf(x,y,z,'FaceColor',cols(i,:),'FaceAlpha',0.1,'EdgeAlpha',0.5);
        area(i,1) = prod(radii)*4/3*pi;
        area(i,2:4) = radii;
    end
end

%% Add factor level centres
if dims == 2
    for i=1:length(levs)
        plot3(mean(projected(where==i,1)), mean(projected(where==i,2)),1, 'o','Color',cols(i,:),'MarkerSize',7,'MarkerFaceColor',[0,0,0],'LineWidth',1.5)
    end
end


%% Multiple comparisons
function mult = multComp(projected, where, ASA, nobj, q, d, levs)
K = meanBy(projected(:,1:2),where);
dists = pdist(K,'mahalanobis',ASA./sqrt((nobj-q)*d./(nobj-q-d+1)));
fc = fcdf(dists,d,nobj-q-d+1)';
A = zeros(levs,levs);
names  = repmat('  -  ', length(fc),1);
col1   = repmat(' ',length(fc),1);
col2   = repmat(' ',length(fc),1);
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
k = 0;
for i=1:levs-1
    for j=i+1:levs
        k = k+1; A(j,i) = fc(k);
        names(k,1) = letters(j);
        names(k,5) = letters(i);
        col1(k,1) = letters(j);
        col2(k,1) = letters(i);
    end
end
A = A+A';
mult = [];
mult.means  = K;
mult.names  = names;
mult.col1   = col1;
mult.col2   = col2;
mult.dists  = dists;
mult.pvalue = 1-fc;
mult.pmat   = 1-A;
