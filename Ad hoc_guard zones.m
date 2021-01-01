%Simulation window parameters
r=10; %radius of disk
 
xx0=0; yy0=0; %centre of disk
areaTotal=pi*r^2; %area of disk
%Point process parameters
lambda=1; %intensity (ie mean density) of the Poisson process
%Simulate Poisson point process
numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
theta=2*pi*(rand(numbPoints,1)); %angular coordinates
rho=r*sqrt(rand(numbPoints,1)); %radial coordinates
%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
%Shift centre of disk to (xx0,yy0)
xx=xx+xx0;
yy=yy+yy0;
tx = [xx,yy];
 
dmax = 1/2.5;
theta2=2*pi*(rand(numbPoints,1));
dMax = dmax*ones(numbPoints,1);
[xxr,yyr]=pol2cart(theta2,dMax);
xx1 = xx+xxr;
yy1 = yy+yyr;
rx = [xx1,yy1];
 
%No point outside the radius
for i=1:1:numbPoints;
    if norm([xx1(i),yy1(i)]-[xx0,yy0])>r;
        while 1
            x = xx1(i)-xx(i);
            y = yy1(i)-yy(i);
            [thet,dis]=cart2pol(x,y);
            thet = 2*pi*rand;
            [x1,y1] = pol2cart(thet,dis);
            xx1(i) = x1+xx(i);
            yy1(i) = y1+yy(i);
            if norm([xx1(i),yy1(i)]-[xx0,yy0])<r;
                break
            end
        end
    end
end
 
All = [xx, yy, xx1, yy1];
%Plotting No Guard:
subplot(2,3,2)
F=scatter(xx,yy,10,'filled');
title(['No Guard Zones, Nodes = ', num2str(numbPoints)]);
xlabel('x');ylabel('y');
axis square;
viscircles([xx0,yy0],r);
hold on
scatter(xx1, yy1, 10, 'x');
hold on
for i=1:1:numbPoints;
    hold all
    plot([xx(i),xx1(i)],[yy(i),yy1(i)],'r')
end
 
 
%With Guards
 
iter = 1;
 
for k =1:iter;
    delta = 0.5;
    D = (1+delta)*dmax;
    % Guard implementation
    
    Tx=[xx(1),yy(1)];
    %Txy=[yy(1)];
    Rx=[xx1(1),yy1(1)];
    %Rxy=[yy1(1)];
    
    for i = 2:1:numbPoints;
        for j = 1:1:size(Tx,1);
            if norm([xx(i),yy(i)]-[Rx(j,1),Rx(j,2)])<D;
                break
            end
            if norm([xx1(i),yy1(i)]-[Tx(j,1),Tx(j,2)])<D;
                break
            end
            if j==(size(Tx,1));
                Tx = [Tx;[xx(i),yy(i)]];
                Rx = [Rx;[xx1(i),yy1(i)]];
            end
        end
    end
 
    %Plotting With Guard:
    txx = Tx(:,1);
    txy = Tx(:,2);
    rxx = Rx(:,1);
    rxy = Rx(:,2);
    subplot(2,3,5)
    scatter(txx,txy,10,'filled');
    % %Virtual rings as guard zones around receivers
    %     for m = 1:1:size(Tx,1)
    %         viscircles([rxx(m),rxy(m)],D,'LineWidth',1);
    %     end
    
    title(['Delta = ', num2str(delta),', Nodes = ', num2str(size(Tx,1))]);
    xlabel('x');ylabel('y');
    axis square;
    viscircles([xx0,yy0],r);
    hold on
    scatter(rxx, rxy, 10, 'x');
    hold on
    for i=1:1:size(Tx,1);
        hold all
        plot([txx(i),rxx(i)],[txy(i),rxy(i)],'r')
    end
    axis([-r r -r r])
end
 
 
%Delta vs No. of nodes & Okumura Hata Model
 
iter = 50;
 
aa = [];
bb = [];
 
dd = [];
f = 960*10^6;
hb = 1.5;
hm = 1.5;
 
Ch = 0.8+(1.1*log(f)-0.7)*hm-1.56*log(f);
p_sum = 0;
ag_sum = 0;
ee = [];
ff = [];
P = 1;
alpha = 4;
OkHM = [];
cc = [];
OkHM_Sb =[];
OkHM_R = [];
 
for k =0:iter;
    OHM = 0;
    delta = 0.1*k-1;
    D = (1+delta)*dmax;
    % Guard implementation
    S = 0;
    SINR = 0;
    Tx=[xx(1),yy(1)];
    %Txy=[yy(1)];
    Rx=[xx1(1),yy1(1)];
    %Rxy=[yy1(1)];
 
    
    for i = 2:1:numbPoints;
        for j = 1:1:size(Tx,1);
            if norm([xx(i),yy(i)]-[Rx(j,1),Rx(j,2)])<D;
                break
            end
            if norm([xx1(i),yy1(i)]-[Tx(j,1),Tx(j,2)])<D;
                break
            end
            if j==(size(Tx,1));
                Tx = [Tx;[xx(i),yy(i)]];
                Rx = [Rx;[xx1(i),yy1(i)]];
            end
        end
    end
    aa = [aa,delta];
    bb = [bb,length(Tx)];
    
    % Okumura Hata Model, Aggregate Interference and SINR
    ee = [];
    for m = 1:length(Tx);
        
        ag_sum = 0;
        p_sum = 0;
        for n = 1:length(Tx);
            d = norm([Rx(m,1),Rx(m,2)]-[Tx(n,1),Tx(n,2)]);
            
            aggregate = P/(d^alpha);
            
            if m ==1;
                Lu = 690.55 + 26.16*log(f) - 13.82*log(hb) - Ch + (44.9-6.55*log(hb))*log(d);
                OHM = OHM+Lu;
            end
            
            if n==m;
                continue
            end
            ag_sum = ag_sum+aggregate;
        end
        
        ee = [ee;[ag_sum]];
        
    end
    
    
    
    avg_ag = mean(ee);
    
    %SINR
    S = P/(dmax^alpha);
    K = 1.3807*10^-23;
    T = 273;
    B = f;
    N = K*B*T;
    SINR = S/(N+avg_ag);
    
    cc = [cc;[SINR]];
        
    if avg_ag>1000;
        avg_ag = 1000;
    end
    ff = [ff;[avg_ag]];
    OHM_Sb = OHM-4.78*log(f)^2-18.33*log(f)-40.94;
    OHM_R = OHM-2*(log(f/28))^2-5.4;
    OkHM_Sb = [OkHM_Sb;[OHM_Sb]];
    d = dmax;
    OkHM = [OkHM;[OHM]];
    OkHM_R = [OkHM_R;[OHM_R]];
    
    
end
    
 
    % Lu = 69.55 + 26.16*log(f) - 13.82*log(hb) - Ch + (44.9-6.55*log(hb))*log(d);
    %
    % own_pl = Lu*ones(length(dd),1);
    %
    % ratio = dd./own_pl;
    
    
    subplot(2,3,1)
    plot(aa,OkHM);
    title('Okumura-Hata Model');
    hold on 
    plot(aa,OkHM_Sb);
    plot(aa,OkHM_R);
    hold off
    ylabel('Path Loss via Okumura Hata Model');
    xlabel('Delta')
    
    subplot(2,3,4)
    plot(aa,ff)
    hold on
    plot(aa,ff,'x');
    hold off
    title('Aggregate Interference');
    xlabel('Delta');
    ylabel('Aggregate Interference');
    %axis([1 4 0 5])
    
    subplot(2,3,6)
    plot(aa, cc);
    hold on
    stem(aa,cc);
    hold off
    title('Signal to Noise plus Interference Ratio')
    xlabel('delta');
    ylabel('SINR');
    
    
    t_t=round(numbPoints*2/3);
    dist    = abs(t_t - bb);
    minDist = min(dist);
    Idx     = find(dist == minDist);
    idx = Idx(1);
    subplot(2,3,3)
    plot(aa,bb)
    title('No. of nodes w.r.t delta')
    xlabel('delta')
    ylabel('No. of nodes')
    hold on
    axis([min(aa) max(aa) 0 max(bb)+5]);
    refline(0,bb(idx))
    hold on
    line([aa(idx) aa(idx)], [0 max(bb)+5]);