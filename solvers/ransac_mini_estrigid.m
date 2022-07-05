function [T,inliers,inliers2] = ransac_mini_estrigid(xx,yy,dd2,bnd,iters)
% function [T,inliers] = ransac_mini_estrigid(xx,yy,dd,bnd,iters)
%
% input:    xx 2xnn local platform coordinates
%           yy 2xnn global receiver coordinates
%           dd 1xnn squared distances
%           bnd inlier bound
%           iters number of ransac iterations
%
% output:   T best transformation [a -b tx;b a ty;0 0 1]
%           inliers 1xnn logical inlier set

if nargin<5
    iters = 1000;
end
if nargin<4
    bnd = 0.05;
end

fac = 10;
bestins = 0;
bestnorm = inf;
nn = size(dd2,2);

for iter = 1:iters
    ids = randperm(nn,3);
    [aa,bb,ttx,tty] = fullsolver_mini_estrigid(xx(:,ids),yy(:,ids),dd2(:,ids));
    nsol = size(aa,2);
    for iii = 1:nsol
        a = aa(iii);
        b = bb(iii);
        tx = ttx(iii);
        ty = tty(iii);
        res = abs(sqrt((a*xx(1,:)-b*xx(2,:)+tx-yy(1,:)).^2+(b*xx(1,:)+a*xx(2,:)+ty-yy(2,:)).^2)-sqrt(dd2));
        resins = sum(res<bnd);
        resnorm = norm(res(res<bnd));
        if resins>bestins
            bestins = resins;
            bestnorm = resnorm;
            T = [a -b tx;b a ty;0 0 1];
        end
        if resins==bestins && resnorm<bestnorm
            bestins = resins;
            bestnorm = resnorm;
            T = [a -b tx;b a ty;0 0 1];
        end
            
    end
end

a = T(1,1);
b = T(2,1);
tx = T(1,3);
ty = T(2,3);

inliers = abs(sqrt((a*xx(1,:)-b*xx(2,:)+tx-yy(1,:)).^2+(b*xx(1,:)+a*xx(2,:)+ty-yy(2,:)).^2)-sqrt(dd2))<bnd;
inliers2 = abs(sqrt((a*xx(1,:)-b*xx(2,:)+tx-yy(1,:)).^2+(b*xx(1,:)+a*xx(2,:)+ty-yy(2,:)).^2)-sqrt(dd2))<(fac*bnd);

%disp(bestnorm)
