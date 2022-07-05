% this script sets up a simulated problem and compares the proposed
% method to dead reckoning and optimal trilateration

%% simulation parameters


rtt_std = 0.1; % meters
step_relerror = 0.1;
mean_steps = 5;
mean_stepstd = 2;
step_std = mean_steps*step_relerror; % no unit
step_length = 0.75;
heading_std = 0.005; % radians
nr_senders = 20;
room_sz = 50;
nr_segments = 10;
nr_receivers = 5;
nr_tot  = nr_segments*nr_receivers;
outliers = 0.1;
nr_meas = 3;

mean_dd2 = 1000;



%% setup simulation

yy = rand(2,nr_senders)*room_sz-room_sz/2;
m_ids = zeros(nr_meas,nr_tot);
for iii = 1:nr_tot
    m_ids(:,iii) = randperm(nr_senders,nr_meas);
end


heads = cumsum(randn(1,nr_tot)/2);

steps = round(mean_steps+rand(1,nr_tot)*mean_stepstd*2-mean_stepstd);
xx = rand(2,nr_tot)*room_sz-room_sz/2;
for iii = 2:nr_tot
    xx(:,iii) = xx(:,iii-1)+steps(iii-1)*step_length*[cos(heads(iii-1));sin(heads(iii-1))];
end

dd2 = zeros(nr_meas,nr_tot);
for iii = 1:nr_tot
    dd2(:,iii) = (sum((xx(:,iii)-yy(:,m_ids(:,iii))).^2))';
end

dd2n = (sqrt(dd2)+randn(nr_meas,nr_tot)*rtt_std).^2;

dn = nr_meas*nr_tot;
nr_out = round(dn*outliers);
dd2n(randperm(dn,nr_out)) = rand(1,nr_out)*mean_dd2;

stepsn = round(steps+randn(1,nr_tot)*step_std);
headsn = heads+randn(1,nr_tot)*heading_std;

%% test opttrilat

xxo = zeros(2,nr_tot);
for iii = 1:nr_tot
    r0 = yy(:,m_ids(:,iii));
    d = sqrt(dd2n(:,iii))';
    xo = solver_opttrilat(r0,d,1./d.^2);
    xxo(:,iii) = xo;
end

%% test dead reckoning
xxn = xx;
for iii = 2:nr_tot
    xxn(:,iii) = xxn(:,iii-1)+stepsn(iii-1)*step_length*[cos(headsn(iii-1));sin(headsn(iii-1))];
end

%% test proposed rigid

bnd = 0.5;
iters = 1000;
xxer = xx;

for iii = 1:nr_segments
    xxl = zeros(2,nr_receivers*nr_meas);
    
    for jjj = 2:nr_receivers
        id = nr_receivers*(iii-1)+(jjj-1);
        xxl(:,(1:nr_meas)+(jjj-1)*nr_meas) = repmat(xxl(:,(jjj-1)*nr_meas)+steps(id)*step_length*[cos(heads(id));sin(heads(id))],1,nr_meas);
    end
    yyl = zeros(2,nr_receivers*nr_meas);
    for jjj = 1:nr_receivers
        id = nr_receivers*(iii-1)+jjj;
        yyl(:,(1:nr_meas)+(jjj-1)*nr_meas) = yy(:,m_ids(:,id));
    end
    
    dd2l = zeros(1,nr_receivers*nr_meas);
    for jjj = 1:nr_receivers
        id = nr_receivers*(iii-1)+jjj;
        dd2l(1,(1:nr_meas)+(jjj-1)*nr_meas) = dd2n(:,id)';
    end
    
    
    [T,inliers] = ransac_mini_estrigid(xxl,yyl,dd2l,bnd,iters);
    
    xe = T(1:2,1:2)*xxl(:,1:nr_meas:end)+T(1:2,3);
    xxer(:,(1:nr_receivers)+(iii-1)*nr_receivers)=xe;
    
end

%% test proposed sim

bnd = 0.5;
iters = 1000;
xxes = xx;

for iii = 1:nr_segments
    xxl = zeros(2,nr_receivers*nr_meas);
    
    for jjj = 2:nr_receivers
        id = nr_receivers*(iii-1)+(jjj-1);
        xxl(:,(1:nr_meas)+(jjj-1)*nr_meas) = repmat(xxl(:,(jjj-1)*nr_meas)+steps(id)*[cos(heads(id));sin(heads(id))],1,nr_meas);
    end
    yyl = zeros(2,nr_receivers*nr_meas);
    for jjj = 1:nr_receivers
        id = nr_receivers*(iii-1)+jjj;
        yyl(:,(1:nr_meas)+(jjj-1)*nr_meas) = yy(:,m_ids(:,id));
    end
    
    dd2l = zeros(1,nr_receivers*nr_meas);
    for jjj = 1:nr_receivers
        id = nr_receivers*(iii-1)+jjj;
        dd2l(1,(1:nr_meas)+(jjj-1)*nr_meas) = dd2n(:,id)';
    end
    
    
    [T,inliers] = ransac_mini_estsim(xxl,yyl,dd2l,bnd,iters);
    
    xe = T(1:2,1:2)*xxl(:,1:nr_meas:end)+T(1:2,3);
    xxes(:,(1:nr_receivers)+(iii-1)*nr_receivers)=xe;
    
end

%% errors 
e_o = rms(xx(:)-xxo(:));
e_n = rms(xx(:)-xxn(:));
e_er = rms(xx(:)-xxer(:));
e_es = rms(xx(:)-xxes(:));

disp('Distance to ground truth (rms)')
disp(['Optimal trilateration: ' num2str(e_o) ' m']);
disp(['Dead reckoning: ' num2str(e_n) ' m']);
disp(['Proposed Euclidean: ' num2str(e_er) ' m']);
disp(['Proposed similarity: ' num2str(e_es) ' m']);


%% plot results


ms = 10;
figure(1);
clf
ll = plot(yy(1,:),yy(2,:),'x');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
hold on
ll = plot(xx(1,:),xx(2,:),'*');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxer(1,:),xxer(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxn(1,:),xxn(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxo(1,:),xxo(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);

for iii = 1:nr_segments
    ll = plot(xx(1,(1:nr_receivers)+nr_receivers*(iii-1)),xx(2,(1:nr_receivers)+nr_receivers*(iii-1)),'-');
    set(ll,'LineWidth',2);
    set(ll,'color','c');
end
legend({'Senders', 'GT receivers','Estimated receivers (proposed Euclidean)',...
    'Estimated receivers (dead reckoning)','Estimated receivers (optimal trilateration)'},'Location','northwest');
title('Proposed Euclidean solver comparison');

axis equal


figure(2);
clf

ids = (nr_receivers*(nr_segments-3)+1):(nr_receivers*nr_segments);
ll = plot(xx(1,ids),xx(2,ids),'*');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
hold on
ll = plot(xxer(1,ids),xxer(2,ids),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxn(1,ids),xxn(2,ids),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);

for iii = (nr_segments-2):nr_segments
    ll = plot(xx(1,(1:nr_receivers)+nr_receivers*(iii-1)),xx(2,(1:nr_receivers)+nr_receivers*(iii-1)),'-');
    set(ll,'LineWidth',2);
    set(ll,'color','c');
end
legend({'GT receivers','Estimated receivers (proposed Euclidean)',...
    'Estimated receivers (dead reckoning)'},'Location','northwest');
axis equal
title('Close up of last three segments (Euclidean case)')


figure(3);
clf
ll = plot(yy(1,:),yy(2,:),'x');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
hold on
ll = plot(xx(1,:),xx(2,:),'*');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxes(1,:),xxes(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxn(1,:),xxn(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxo(1,:),xxo(2,:),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);

for iii = 1:nr_segments
    ll = plot(xx(1,(1:nr_receivers)+nr_receivers*(iii-1)),xx(2,(1:nr_receivers)+nr_receivers*(iii-1)),'-');
    set(ll,'LineWidth',2);
    set(ll,'color','c');
end
legend({'Senders', 'GT receivers','Estimated receivers (proposed similarity)',...
    'Estimated receivers (dead reckoning)','Estimated receivers (optimal trilateration)'},'Location','northwest');
axis equal
title('Proposed similarity solver comparison');

figure(4);
clf

ids = (nr_receivers*(nr_segments-3)+1):(nr_receivers*nr_segments);
ll = plot(xx(1,ids),xx(2,ids),'*');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
hold on
ll = plot(xxes(1,ids),xxes(2,ids),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);
ll = plot(xxn(1,ids),xxn(2,ids),'o');
set(ll,'MarkerSize',ms);
set(ll,'LineWidth',2);

for iii = (nr_segments-2):nr_segments
    ll = plot(xx(1,(1:nr_receivers)+nr_receivers*(iii-1)),xx(2,(1:nr_receivers)+nr_receivers*(iii-1)),'-');
    set(ll,'LineWidth',2);
    set(ll,'color','c');
end
legend({'GT receivers','Estimated receivers (proposed similarity)',...
    'Estimated receivers (dead reckoning)'},'Location','northwest');
axis equal
title('Close up of last three segments (similarity case)')

