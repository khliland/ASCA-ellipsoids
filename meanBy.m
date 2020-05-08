function [means,orig] = meanBy(X,y, prc)
%% Calculate means of X by groups in y
% [means,orig] = meanBy(X,y)

[~,~,y] = unique(y,'rows');
n = max(y);
p = size(X,2);
means = zeros(n,p);

for i=1:n
    if sum(y==i) > 1
        if nargin == 3
            means(i,:) = trimmean(X(y==i,:),prc);
        else
            means(i,:) = nanmean(X(y==i,:));
        end
    else
        means(i,:) = X(y==i,:);
    end
end

if nargout == 2
    orig = zeros(n,p);
    for i=1:n
        orig(y==i,:) = repmat(means(i,:),sum(y==i),1);
    end
end