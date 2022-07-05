function [aa,bb,ttx,tty] = fullsolver_mini_estrigid(xi,yi,di2)
% function [a,b,tx,ty] = fullsolver_mini_estrigid(xi,yi,di2)
%
% input xi: 2x3 points to be transformed
%       yi: 2x3 points in world
%       di2: squared distances between points
% output transformation parameters 1xnsol 
% T = [a -b tx;b a ty;0 0 1]

lille = 1e-9;
lille2 = 1e-5;
cc = [xi;yi;di2];
[data,cc2] = getdata_mini_estrigid(cc(:));
sol = solver_mini_estrigid_v2(data);
sol(:,sum(abs(imag(sol)))>lille)=[];
sol(:,abs(sol(1,:).^2+sol(2,:).^2-1)>lille2)=[];
nsol = size(sol,2);
aa = sol(1,:);
bb = sol(2,:);
ttx = aa;
tty = aa;
for iii = 1:nsol
    a = aa(iii);
    b = bb(iii);
    ttx(iii) = (cc2(8)*cc2(18) - cc2(9)*cc2(17) - cc2(8)*cc2(27) + cc2(9)*cc2(26) + cc2(17)*cc2(27) - cc2(18)*cc2(26) + a*cc2(2)*cc2(18) - a*cc2(3)*cc2(17) + a*cc2(8)*cc2(12) - a*cc2(9)*cc2(11) - a*cc2(2)*cc2(27) + a*cc2(3)*cc2(26) - a*cc2(8)*cc2(21) + a*cc2(9)*cc2(20) + a*cc2(11)*cc2(27) - a*cc2(12)*cc2(26) + a*cc2(17)*cc2(21) - a*cc2(18)*cc2(20) + b*cc2(5)*cc2(18) - b*cc2(6)*cc2(17) + b*cc2(8)*cc2(15) - b*cc2(9)*cc2(14) - b*cc2(5)*cc2(27) + b*cc2(6)*cc2(26) - b*cc2(8)*cc2(24) + b*cc2(9)*cc2(23) + b*cc2(14)*cc2(27) - b*cc2(15)*cc2(26) + b*cc2(17)*cc2(24) - b*cc2(18)*cc2(23) + a^2*cc2(2)*cc2(12) - a^2*cc2(3)*cc2(11) - a^2*cc2(2)*cc2(21) + a^2*cc2(3)*cc2(20) + a^2*cc2(11)*cc2(21) - a^2*cc2(12)*cc2(20) + b^2*cc2(5)*cc2(15) - b^2*cc2(6)*cc2(14) - b^2*cc2(5)*cc2(24) + b^2*cc2(6)*cc2(23) + b^2*cc2(14)*cc2(24) - b^2*cc2(15)*cc2(23) + a*b*cc2(2)*cc2(15) - a*b*cc2(3)*cc2(14) + a*b*cc2(5)*cc2(12) - a*b*cc2(6)*cc2(11) - a*b*cc2(2)*cc2(24) + a*b*cc2(3)*cc2(23) - a*b*cc2(5)*cc2(21) + a*b*cc2(6)*cc2(20) + a*b*cc2(11)*cc2(24) - a*b*cc2(12)*cc2(23) + a*b*cc2(14)*cc2(21) - a*b*cc2(15)*cc2(20))/(cc2(7)*cc2(17) - cc2(8)*cc2(16) - cc2(7)*cc2(26) + cc2(8)*cc2(25) + cc2(16)*cc2(26) - cc2(17)*cc2(25) + a*cc2(1)*cc2(17) - a*cc2(2)*cc2(16) + a*cc2(7)*cc2(11) - a*cc2(8)*cc2(10) - a*cc2(1)*cc2(26) + a*cc2(2)*cc2(25) - a*cc2(7)*cc2(20) + a*cc2(8)*cc2(19) + a*cc2(10)*cc2(26) - a*cc2(11)*cc2(25) + a*cc2(16)*cc2(20) - a*cc2(17)*cc2(19) + b*cc2(4)*cc2(17) - b*cc2(5)*cc2(16) + b*cc2(7)*cc2(14) - b*cc2(8)*cc2(13) - b*cc2(4)*cc2(26) + b*cc2(5)*cc2(25) - b*cc2(7)*cc2(23) + b*cc2(8)*cc2(22) + b*cc2(13)*cc2(26) - b*cc2(14)*cc2(25) + b*cc2(16)*cc2(23) - b*cc2(17)*cc2(22) + a^2*cc2(1)*cc2(11) - a^2*cc2(2)*cc2(10) - a^2*cc2(1)*cc2(20) + a^2*cc2(2)*cc2(19) + a^2*cc2(10)*cc2(20) - a^2*cc2(11)*cc2(19) + b^2*cc2(4)*cc2(14) - b^2*cc2(5)*cc2(13) - b^2*cc2(4)*cc2(23) + b^2*cc2(5)*cc2(22) + b^2*cc2(13)*cc2(23) - b^2*cc2(14)*cc2(22) + a*b*cc2(1)*cc2(14) - a*b*cc2(2)*cc2(13) + a*b*cc2(4)*cc2(11) - a*b*cc2(5)*cc2(10) - a*b*cc2(1)*cc2(23) + a*b*cc2(2)*cc2(22) - a*b*cc2(4)*cc2(20) + a*b*cc2(5)*cc2(19) + a*b*cc2(10)*cc2(23) - a*b*cc2(11)*cc2(22) + a*b*cc2(13)*cc2(20) - a*b*cc2(14)*cc2(19));
    tty(iii) = -(cc2(7)*cc2(18) - cc2(9)*cc2(16) - cc2(7)*cc2(27) + cc2(9)*cc2(25) + cc2(16)*cc2(27) - cc2(18)*cc2(25) + a*cc2(1)*cc2(18) - a*cc2(3)*cc2(16) + a*cc2(7)*cc2(12) - a*cc2(9)*cc2(10) - a*cc2(1)*cc2(27) + a*cc2(3)*cc2(25) - a*cc2(7)*cc2(21) + a*cc2(9)*cc2(19) + a*cc2(10)*cc2(27) - a*cc2(12)*cc2(25) + a*cc2(16)*cc2(21) - a*cc2(18)*cc2(19) + b*cc2(4)*cc2(18) - b*cc2(6)*cc2(16) + b*cc2(7)*cc2(15) - b*cc2(9)*cc2(13) - b*cc2(4)*cc2(27) + b*cc2(6)*cc2(25) - b*cc2(7)*cc2(24) + b*cc2(9)*cc2(22) + b*cc2(13)*cc2(27) - b*cc2(15)*cc2(25) + b*cc2(16)*cc2(24) - b*cc2(18)*cc2(22) + a^2*cc2(1)*cc2(12) - a^2*cc2(3)*cc2(10) - a^2*cc2(1)*cc2(21) + a^2*cc2(3)*cc2(19) + a^2*cc2(10)*cc2(21) - a^2*cc2(12)*cc2(19) + b^2*cc2(4)*cc2(15) - b^2*cc2(6)*cc2(13) - b^2*cc2(4)*cc2(24) + b^2*cc2(6)*cc2(22) + b^2*cc2(13)*cc2(24) - b^2*cc2(15)*cc2(22) + a*b*cc2(1)*cc2(15) - a*b*cc2(3)*cc2(13) + a*b*cc2(4)*cc2(12) - a*b*cc2(6)*cc2(10) - a*b*cc2(1)*cc2(24) + a*b*cc2(3)*cc2(22) - a*b*cc2(4)*cc2(21) + a*b*cc2(6)*cc2(19) + a*b*cc2(10)*cc2(24) - a*b*cc2(12)*cc2(22) + a*b*cc2(13)*cc2(21) - a*b*cc2(15)*cc2(19))/(cc2(7)*cc2(17) - cc2(8)*cc2(16) - cc2(7)*cc2(26) + cc2(8)*cc2(25) + cc2(16)*cc2(26) - cc2(17)*cc2(25) + a*cc2(1)*cc2(17) - a*cc2(2)*cc2(16) + a*cc2(7)*cc2(11) - a*cc2(8)*cc2(10) - a*cc2(1)*cc2(26) + a*cc2(2)*cc2(25) - a*cc2(7)*cc2(20) + a*cc2(8)*cc2(19) + a*cc2(10)*cc2(26) - a*cc2(11)*cc2(25) + a*cc2(16)*cc2(20) - a*cc2(17)*cc2(19) + b*cc2(4)*cc2(17) - b*cc2(5)*cc2(16) + b*cc2(7)*cc2(14) - b*cc2(8)*cc2(13) - b*cc2(4)*cc2(26) + b*cc2(5)*cc2(25) - b*cc2(7)*cc2(23) + b*cc2(8)*cc2(22) + b*cc2(13)*cc2(26) - b*cc2(14)*cc2(25) + b*cc2(16)*cc2(23) - b*cc2(17)*cc2(22) + a^2*cc2(1)*cc2(11) - a^2*cc2(2)*cc2(10) - a^2*cc2(1)*cc2(20) + a^2*cc2(2)*cc2(19) + a^2*cc2(10)*cc2(20) - a^2*cc2(11)*cc2(19) + b^2*cc2(4)*cc2(14) - b^2*cc2(5)*cc2(13) - b^2*cc2(4)*cc2(23) + b^2*cc2(5)*cc2(22) + b^2*cc2(13)*cc2(23) - b^2*cc2(14)*cc2(22) + a*b*cc2(1)*cc2(14) - a*b*cc2(2)*cc2(13) + a*b*cc2(4)*cc2(11) - a*b*cc2(5)*cc2(10) - a*b*cc2(1)*cc2(23) + a*b*cc2(2)*cc2(22) - a*b*cc2(4)*cc2(20) + a*b*cc2(5)*cc2(19) + a*b*cc2(10)*cc2(23) - a*b*cc2(11)*cc2(22) + a*b*cc2(13)*cc2(20) - a*b*cc2(14)*cc2(19));
end


