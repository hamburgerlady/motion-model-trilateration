function [T,inliers] = ransac_mini_estsim(xx,yy,dd2,bnd,iters,linefix)
% function [T,inliers] = ransac_mini_estsim(xx,yy,dd,bnd,iters)
%
% input:    xx 2xnn local platform coordinates
%           yy 2xnn global receiver coordinates
%           dd 1xnn squared distances
%           bnd inlier bound
%           iters number of ransac iterations
%
% output:   T best transformation [sa -sb tx;sb sa ty;0 0 1]
%           inliers 1xnn logical inlier set
if nargin<6
    linefix = 0;
end
if nargin<5
    iters = 1000;
end
if nargin<4
    bnd = 0.05;
end
lille = 1e-9;
bestins = 0;
nn = size(dd2,2);

for iter = 1:iters
    if linefix>0
        ids = randi(linefix,1,4) + [0 linefix 2*linefix 3*linefix];
    else
    ids = randperm(nn,4);
    end
    data = [xx(:,ids);yy(:,ids);dd2(:,ids)];
    sols = solver_mini_estsim(data(:));
    sols(:,sum(abs(imag(sols)))>lille)=[];
    nsol = size(sols,2);
    for iii = 1:nsol
        a = sols(1,iii);
        b = sols(2,iii);
        tx = sols(3,iii);
        ty = sols(4,iii);
        res = abs(sqrt((a*xx(1,:)-b*xx(2,:)+tx-yy(1,:)).^2+(b*xx(1,:)+a*xx(2,:)+ty-yy(2,:)).^2)-sqrt(dd2));
        resins = sum(res<bnd);
        if resins>bestins
            bestins = resins;
            T = [a -b tx;b a ty;0 0 1];
        end
    end
end

a = T(1,1);
b = T(2,1);
tx = T(1,3);
ty = T(2,3);

inliers = abs(sqrt((a*xx(1,:)-b*xx(2,:)+tx-yy(1,:)).^2+(b*xx(1,:)+a*xx(2,:)+ty-yy(2,:)).^2)-sqrt(dd2))<bnd;

