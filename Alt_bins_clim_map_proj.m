% bin and plot Atlantic storms vs el nino and NAO winters.


%Bin Storm tracks
clear 
close all
addpath(genpath('/Users/johnlodise/Desktop/ERA5/m_map'));

Files=dir('Alt_stms_v2_*.mat'); 
filenames = {Files.name};

Nino_data = readtable('meiv2.txt');

Nino_idx = table2array(Nino_data(1:42,1:13));
Nino_1d = reshape(Nino_idx(:,2:13)',[],1);

ONI = readtable('ONI_Nino_index1979_2020.txt'); %Oceanic Niño Index
ino_idx = table2array(ONI(1:42,1:13));
ino_1d = reshape(ino_idx(:,2:13)',[],1);

PDO = readtable('ersst.v5.pdo.txt');
PDO_t = table2array(PDO(126:167,1:13));
PDO_1d = reshape(PDO_t(:,2:13)',[],1);


NAO = readtable('NAO_index1950_2021.txt');
NAO_t = table2array(NAO(30:end-1,1:13));
NAO_1d = reshape(NAO_t(:,2:13)',[],1);

figure(10)
hold on
plot(ino_1d)
plot(Nino_1d)
plot(PDO_1d)
plot(NAO_1d)
% for m = 1:42
%    Nyr(m)  =  Nino_idx(m,1) ;
%    
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read in monthy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat_all = [];
lon_all = [];
u_all = [];
v_all = [];
Lap_press_all = [];
press_all = [];
time_all = [];
time_all_RoC = [];
num_storms = 0;
num_trops = 0;
min_press = [];
dur_stm = [];

u2 = [];
v2 = [];

press_RoC = [];
lat_RoC = [];
lon_RoC = [];

storms_PM = NaN(1,length(Files));

for ii=1:length(Files) %where N is the number of images
 
    fileName = Files(ii).name;
    load(fileName)
    
    for jj = 1:length(A_stms)
        if A_stms(jj).trop_idx == 1 %if storm originiates in the tropics do not include 
            num_trops = num_trops + 1;
            continue
        end
%         
        if ~any(A_stms(jj).mslp <= 980,'all') %if storm doesnt go below 960 exclude  
            %DC = DC + 1;
            continue
        end
        
        num_storms = num_storms+1;      
        lat_all = cat(2,lat_all,A_stms(jj).lat');
        lon_all = cat(2,lon_all,A_stms(jj).lon');
        u_all = cat(2,u_all,A_stms(jj).u');
        v_all = cat(2,v_all,A_stms(jj).v');
        press_all = cat(2,press_all,A_stms(jj).mslp');
        Lap_press_all = cat(2,Lap_press_all,A_stms(jj).mslp_Lap');
        time_all = cat(2,time_all,A_stms(jj).time'); 
         
        press = A_stms(jj).mslp';
        min_press = cat(2,min_press,min(press));
       
        lond =  A_stms(jj).lon';
        latd = A_stms(jj).lat';
        dur_stm = cat(2,dur_stm,length(lond));
        
        press_RoC = cat(2,press_RoC,diff(press));
        lat_RoC = cat(2,lat_RoC,latd(1:end-1));
        lon_RoC = cat(2,lon_RoC,lond(1:end-1));
        time2 = A_stms(jj).time';
        time_all_RoC = cat(2,time_all_RoC,time2(1:end-1));
       
        
%         stm_time = A_stms(jj).time;
%         time_str = datestr(stm_time); 
%        
%         N_idx = cat(2,N_idx)     
    end
    
    storms_PM(ii) = jj;  %note number of storms per month  
end


figure(1)
hold on
plot(Nino_1d,storms_PM,'b.','Markersize',15)%,'Linestyle','none')
%plot(ino_1d,storms_PM,'r.','Markersize',15)%,'Linestyle','none')

figure(11)
hold on
plot(PDO_1d,storms_PM,'b.','Markersize',15)%,'Linestyle','none') 

% lon_all(lat_all < 25) = [];
% u_all(lat_all < 25) = [];
% v_all(lat_all < 25) = [];
% press_all(lat_all < 25) = [];
% time_all(lat_all < 25) = [];
% lat_all(lat_all < 25) = [];
u_all = u_all.*3.6;
v_all = v_all.*3.6;
Vmag_all = sqrt(u_all.^2 + v_all.^2);

ybins =   15:2:90;
xbins = -100:2:40;

[lon_bins,lat_bins] = meshgrid(xbins,ybins);
% 
% ddistyy = NaN(size(lon_bins));
% ddistxx = NaN(size(lat_bins));
% a=6378.137;
% b=6356.7523;   %Earth radii at equator and poles.
% 
% for i = 1:length(ybins)-1
%     lat0 = lat_bins(i,:)*(pi/180);    %approx lat of the data-set.  U can change it to the drifter lat.
%     R=sqrt(((  (a^2).*cos(lat0)).^2+((b^2).*sin(lat0)).^2)./   ((a.*cos(lat0)).^2 + (b.*sin(lat0)).^2)); 
%     
%     long = lon_bins(i,:);
%     d1=2*pi.*R./360; %length of 1 degree (km) meridional
%   %  d1=d1*1.e3;    % convert to meters
% 
%     ddistxx(i,1:end-1)=diff(long).*d1(1:end-1).*cos(lat0(1:end-1));
%     ddistxx(i,end)=ddistxx(i,end-1); 
% end
% 
% %ddistxx(1,:) = 1;
% %ddistxx(end,:) = ddistxx(end-1,:);
% 
% 
% lat00 = lat_bins(:,1)*(pi/180);
% R=sqrt((((a^2).*cos(lat00)).^2+((b^2).*sin(lat00)).^2)./   ((a.*cos(lat00)).^2 + (b.*sin(lat00)).^2)); 
% d1=2*pi.*R./360; %length of 1 degree (km) meridional
% %d1=d1*1.e3; 
% ddistyy(1:end-1,:) = diff(lat_bins,1,1).*d1(1:end-1); 
% ddistyy(end,:) = ddistyy(end-1,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Area_bin = ddistxx .* ddistyy;

[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all,lat_all,xbins,ybins); 
% N = N./Area_bin';
% N = N;
%N(N==0) = NaN;
%N = N/42;
load coastlines
figure(88)
hold on
set(gcf,'color','w')
title('1979 - 2020')
set(gca,'FontSize',30)
h = histogram2(lon_all,lat_all,xbins,ybins,'DisplayStyle','tile','FaceColor','flat'); 
colorbar
caxis([0 600])
box on
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])

[X,Y] = meshgrid(mean([Xedges(1:end-1);Xedges(2:end)],1),...
    mean([Yedges(1:end-1);Yedges(2:end,1)],1));


figure(99)
levs = 0:.005:.1;
hold on
set(gcf,'color','w')
title('1979 - 2020')
set(gca,'FontSize',30)
contourf(X,Y,N',levs,'Linestyle','none')
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
colorbar
caxis([0 .1])
box on


time_str = datestr(time_all,'mm-dd-yyyy HH:MM:SS');
idt = time_str(:,1:2);
idt_yr = time_str(:,7:10);
int_mm = NaN(1,length(idt(1,:)));
int_yr= NaN(1,length(idt(1,:)));
N_idx= NaN(1,length(idt(1,:)));
PDO_idx= NaN(1,length(idt(1,:)));
NAO_idx= NaN(1,length(idt(1,:)));

for i = 1:length(idt)
    int_mm(i) = str2double(idt(i,:));
    int_yr(i) = str2double(idt_yr(i,:));
    N_idx(i) = Nino_idx( Nino_idx(:,1) == int_yr(i), int_mm(i)+1);           
    PDO_idx(i) = PDO_t(PDO_t(:,1) == int_yr(i),int_mm(i)+1);
    NAO_idx(i) = NAO_t(NAO_t(:,1) == int_yr(i),int_mm(i)+1);
end

time_str2 = datestr(time_all_RoC,'mm-dd-yyyy HH:MM:SS');
idt2 = time_str2(:,1:2);
idt_yr2 = time_str2(:,7:10);
int_mmRoC = NaN(1,length(idt2(1,:)));
int_yrRoC = NaN(1,length(idt2(1,:)));

for i = 1:length(idt2)
    int_mmRoC(i) = str2double(idt2(i,:));
    int_yrRoC(i) = str2double(idt_yr2(i,:));
end


mo_array = repmat(1:1:12,[42 1]);
mo_array_1d = reshape(mo_array',[],1);

mo_win = find(mo_array_1d < 3 | mo_array_1d > 11);
Nino_1d_w = Nino_1d(mo_win);

num_nino = length(Nino_1d(Nino_1d_w >= 0.5) );
num_nina = length(Nino_1d(Nino_1d_w <= 0.5) );

%season_str = ['Winter (DJF)'; 'Spring (MAM)'; 'Summer (JJA)'; 'Fall (SON) ' ];

season_num = [12 1 2 3 4 5 6 7 8 9 10 11];

tt = 1; 
mm = find(int_mm == season_num(tt) | int_mm ==  season_num(tt+1) | int_mm == season_num(tt+2) );

nno = find(N_idx >= 0.5);
nna = find(N_idx <= 0.5);

nao_p = find(NAO_idx > 0);
nao_n = find(NAO_idx < 0);

nno_w = intersect(mm,nno);
nna_w = intersect(mm,nna);

nao_p_w = intersect(mm,nao_p);
nao_n_w = intersect(mm,nao_n);

% 
% 
% figure(2)
% set(gca,'FontSize',30)
% [N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nno_w),lat_all(nno_w),xbins,ybins);
% hold on
% set(gcf,'color','w')
% title('El Nino winter total' )
% hold on
% set(gcf,'color','w')
% set(gca,'FontSize',30)
% contourf(X,Y,N'./num_nino,'Linestyle','none')
% plot(coastlon,coastlat,'k','Linewidth',2)
% axis([-110 60 10 90 ])
% colorbar
% caxis([0 2.5])
% box on

[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nno_w),lat_all(nno_w),xbins,ybins);
figure(2)
subplot(2,1,1)
hold on
set(gcf,'color','w')
title('El Nino winter total' )
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
m_proj('lambert','long',[-100 40],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(X,Y,N'./num_nino)
m_coast('patch',[1 1 1])
%m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
caxis([0 2.5])
a = colorbar;
a.Label.String = 'Avg storm density / month';
colormap(flipud(hot));

subplot(2,1,2)
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
histogram(press_all(nno_w),40,'normalization', 'probability')
    

press_nno_w = press_all(nno_w);
press_bins_nno_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_nno_w(l,m) =  mean(press_nno_w(aa),'omitnan');
        end

    end
end

figure(3)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('El Nino Winter Ave Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_nno_w',levs,'Linestyle','none');
h = colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ])

u_bins_nno_w = NaN(length(xbins),length(ybins));
v_bins_nno_w =  NaN(length(xbins),length(ybins));
vmag_bins_nno_w =  NaN(length(xbins),length(ybins));

u_nno_w = u_all(nno_w);
v_nno_w = v_all(nno_w);
Vmag_nno_w = Vmag_all(nno_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_nno_w(l,m) =  mean(u_nno_w(aa),'omitnan');
            v_bins_nno_w(l,m) =  mean(v_nno_w(aa),'omitnan');
            vmag_bins_nno_w(l,m) =  mean(Vmag_nno_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 2;
figure(33)
hold on
set(gcf,'color','w')
title('El Nino winter Avg Speed  ')
set(gca,'FontSize',30)
Levs = 0:.4:20;
contourf(xx,yy,vmag_bins_nno_w',Levs,'Linestyle','none')
hold on
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
h =colorbar;
% colormap(redblue)
set(get(h,'title'),'string','ms^-^1');
caxis([0 20])
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_nno_w(1:skip:end,1:skip:end)',v_bins_nno_w(1:skip:end,1:skip:end)','k','Linewidth',2)



%%%%%%%%%%%%%%%%%%%% La Nina Winters Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(4)
% %levs = .0001:.1:.55;
% set(gca,'FontSize',30)
% [N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nna_w),lat_all(nna_w),xbins,ybins);
% hold on
% set(gcf,'color','w')
% title('La Nina winter total')
% hold on
% set(gcf,'color','w')
% set(gca,'FontSize',30)
% contourf(X,Y,N'./num_nina,'Linestyle','none')
% plot(coastlon,coastlat,'k','Linewidth',2)
% axis([-99.5 39.50 15 87 ])
% colorbar
% caxis([0 2.5])
% box on


[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nna_w),lat_all(nna_w),xbins,ybins);
figure(4)
subplot(2,1,1)
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
title('La Nina winter total')
m_proj('lambert','long',[-100 40],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(X,Y,N'./num_nina)
m_coast('patch',[1 1 1])
%m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
caxis([0 2.5])
a = colorbar;
a.Label.String = 'Avg storm density / month';
colormap(flipud(hot));

subplot(2,1,2)
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
histogram(press_all(nna_w),40,'normalization', 'probability')


press_nna_w = press_all(nna_w);
press_bins_nna_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_nna_w(l,m) =  mean(press_nna_w(aa),'omitnan');
        end

    end
end

figure(5)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('La Nina Winter Ave Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_nna_w',levs,'Linestyle','none');
h= colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ])


u_bins_nna_w = NaN(length(xbins),length(ybins));
v_bins_nna_w =  NaN(length(xbins),length(ybins));
vmag_bins_nna_w =  NaN(length(xbins),length(ybins));

u_nna_w = u_all(nna_w);
v_nna_w = v_all(nna_w);
Vmag_nna_w = Vmag_all(nna_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_nna_w(l,m) =  mean(u_nna_w(aa),'omitnan');
            v_bins_nna_w(l,m) =  mean(v_nna_w(aa),'omitnan');
            vmag_bins_nna_w(l,m) =  mean(Vmag_nna_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 2;
figure(55)
hold on
set(gcf,'color','w')
title('La Nina winter Avg Speed  ')
set(gca,'FontSize',30)
Levs = 0:.4:20;
contourf(xx,yy,vmag_bins_nna_w',Levs,'Linestyle','none')
hold on
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
h =colorbar;
% colormap(redblue)
set(get(h,'title'),'string','ms^-^1');
caxis([0 20])
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_nna_w(1:skip:end,1:skip:end)',v_bins_nna_w(1:skip:end,1:skip:end)','k','Linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nino_press_diff_plot%%%%%%%%%%%%%%%%%%%%%%

% figure(555)
% levs =-40:1:40;
% hold on
% set(gcf,'color','w')
% title('La Nina Ave Pressures - El Nino Ave Pressures ')
% set(gca,'FontSize',30)
% contourf(lon_bins,lat_bins,(press_bins_nna_w'- press_bins_nno_w'),levs,'Linestyle','none');
% h= colorbar;
% set(get(h,'title'),'string','m s^-^1');
% %caxis([950 1020])
% plot(coastlon,coastlat,'k','Linewidth',2)
% colormap(redblue)
% box on
% axis([-110 60 10 90 ])
% caxis([-40 40])

figure(555)
hold on
set(gcf,'color','w')
% title('Pos NAO Ave Pressures - NegNAO Ave Pressures ')
set(gca,'FontSize',30)
m_proj('lambert','long',[-100 40],'lat',[24 87]);
m_pcolor(lon_bins,lat_bins,(press_bins_nno_w'- press_bins_nna_w'))
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_coast('patch',[1 1 1])
m_grid('box','fancy','tickdir','in');

colormap(redblue)
h= colorbar;
set(get(h,'title'),'string','hPa');
box on
caxis([-35 35])





%%%%%%%%%%%%%%%%% NAO Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mo_win = find(mo_array_1d < 3 | mo_array_1d > 11);
NAO_win_1d = NAO_1d(mo_win);

num_NAO_p = length(NAO_win_1d(NAO_win_1d > 0) );
num_NAO_n = length(NAO_win_1d(NAO_win_1d < 0) );


% figure(6)
% set(gca,'FontSize',30)
% [N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nao_p_w),lat_all(nao_p_w),xbins,ybins);
% hold on
% set(gcf,'color','w')
% title('Positive NAO winter total' )
% hold on
% set(gcf,'color','w')
% set(gca,'FontSize',30)
% contourf(X,Y,N'./num_NAO_p,'Linestyle','none')
% plot(coastlon,coastlat,'k','Linewidth',2)
% axis([-110 60 10 90 ])
% colorbar
% caxis([0 .55])
% box on

[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nao_p_w),lat_all(nao_p_w),xbins,ybins);
figure(6)
%subplot(2,1,1)
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
m_proj('lambert','long',[-90 30],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(X,Y,N'./num_NAO_p)
m_coast('patch',[1 1 1])
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
caxis([0 7])
a = colorbar;
a.Label.String = 'Avg storm density / month';
colormap(flipud(hot));

% 
% subplot(2,1,2)
% hold on
% set(gcf,'color','w')
% %title('Yearly Averaged Storm Density 1979 - 2020')
% set(gca,'FontSize',30)
% histogram(press_all(nao_p_w),40,'normalization', 'probability')

press_nao_p_w = press_all(nao_p_w);
press_bins_nao_p_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_nao_p_w(l,m) =  mean(press_nao_p_w(aa),'omitnan');
        end

    end
end

figure(7)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('Positive NAO Ave Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_nao_p_w',levs,'Linestyle','none');
h = colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ])


u_bins_nao_p_w = NaN(length(xbins),length(ybins));
v_bins_nao_p_w =  NaN(length(xbins),length(ybins));
vmag_bins_nao_p_w =  NaN(length(xbins),length(ybins));

u_nao_p_w = u_all(nao_p_w);
v_nao_p_w = v_all(nao_p_w);
Vmag_nao_p_w = Vmag_all(nao_p_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_nao_p_w(l,m) =  mean(u_nao_p_w(aa),'omitnan');
            v_bins_nao_p_w(l,m) =  mean(v_nao_p_w(aa),'omitnan');
            vmag_bins_nao_p_w(l,m) =  mean(Vmag_nao_p_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 1;
figure(77)
hold on
set(gcf,'color','w')
%title('Positive NAO ')
% set(gca,'FontSize',30)
% Levs = 0:.4:16;
% contourf(xx,yy,vmag_bins_nao_p_w',Levs,'Linestyle','none')
% hold on
% plot(coastlon,coastlat,'k','Linewidth',2)
% axis([-110 60 10 90 ])
% h =colorbar;
% % colormap(redblue)
% set(get(h,'title'),'string','ms^-^1');
% caxis([0 16])
% quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_nao_p_w(1:skip:end,1:skip:end)',v_bins_nao_p_w(1:skip:end,1:skip:end)','k','Linewidth',2)

hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
m_proj('lambert','long',[-90 30],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(xx,yy,vmag_bins_nao_p_w')
%m_plot(coastlon,coastlat,'k','Linewidth',2)
m_quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_nao_p_w(1:skip:end,1:skip:end)',v_bins_nao_p_w(1:skip:end,1:skip:end)','k','Linewidth',2)
m_coast('patch',[1 1 1])
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
a = colorbar;
a.Label.String = 'km hr^-^1';
caxis([0 70])


%%%%%%%%%%%%%%%%%%%% NAO Negative Winters Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% figure(8)
% set(gca,'FontSize',30)
% [N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nao_n_w),lat_all(nao_n_w),xbins,ybins);
% hold on
% set(gcf,'color','w')
% title('Negative NAO winter total')
% hold on
% set(gcf,'color','w')
% set(gca,'FontSize',30)
% contourf(X,Y,N'./num_NAO_n,'Linestyle','none')
% plot(coastlon,coastlat,'k','Linewidth',2)
% axis([-110 60 10 90 ])
% colorbar
% caxis([0 .55])
% box on



[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(nao_n_w),lat_all(nao_n_w),xbins,ybins);
figure(8)
% subplot(2,1,1)
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
m_proj('lambert','long',[-90 30],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(X,Y,N'./num_NAO_n)
m_coast('patch',[1 1 1])
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
caxis([0 7])
a = colorbar;
a.Label.String = 'Avg storm density / month';
colormap(flipud(hot));

% 
% subplot(2,1,2)
% hold on
% set(gcf,'color','w')
% %title('Yearly Averaged Storm Density 1979 - 2020')
% set(gca,'FontSize',30)
% histogram(press_all(nao_n_w),40,'normalization', 'probability')
% 


press_nao_n_w = press_all(nao_n_w);
press_bins_nao_n_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_nao_n_w(l,m) =  mean(press_nao_n_w(aa),'omitnan');
        end

    end
end


figure(9)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('Negative NAO Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_nao_n_w',levs,'Linestyle','none');
h= colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ])

figure(999)
hold on
set(gcf,'color','w')
% title('Pos NAO Ave Pressures - NegNAO Ave Pressures ')
set(gca,'FontSize',30)
m_proj('lambert','long',[-100 40],'lat',[24 87]);
m_pcolor(lon_bins,lat_bins,(press_bins_nao_p_w'- press_bins_nao_n_w'))
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');

colormap(redblue)
h= colorbar;
set(get(h,'title'),'string','hPa');
box on
caxis([-30 30])



u_bins_nao_n_w = NaN(length(xbins),length(ybins));
v_bins_nao_n_w =  NaN(length(xbins),length(ybins));
vmag_bins_nao_n_w =  NaN(length(xbins),length(ybins));

u_nao_n_w = u_all(nao_n_w);
v_nao_n_w = v_all(nao_n_w);
Vmag_nao_n_w = Vmag_all(nao_n_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_nao_n_w(l,m) =  mean(u_nao_n_w(aa),'omitnan');
            v_bins_nao_n_w(l,m) =  mean(v_nao_n_w(aa),'omitnan');
            vmag_bins_nao_n_w(l,m) =  mean(Vmag_nao_n_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 1;
figure(99)
hold on
set(gcf,'color','w')
%title('Negative NAO ')
hold on
set(gcf,'color','w')
%title('Yearly Averaged Storm Density 1979 - 2020')
set(gca,'FontSize',30)
m_proj('lambert','long',[-90 30],'lat',[24 87]);
%Levs = 0:.5:20;
%m_contourf(X,Y,N',Levs,'Linestyle','none');
m_pcolor(xx,yy,vmag_bins_nao_n_w')
% m_plot(coastlon,coastlat,'k','Linewidth',2)
m_quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_nao_n_w(1:skip:end,1:skip:end)',v_bins_nao_n_w(1:skip:end,1:skip:end)','k','Linewidth',2)
m_coast('patch',[1 1 1])
m_plot(coastlon,coastlat,'k','Linewidth',2)
m_grid('box','fancy','tickdir','in');
colorbar;
a = colorbar;
a.Label.String = 'km hr^-^1';
caxis([0 70])



%%%%%%%%%%%%%%% PDO Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pdo_p = find(PDO_idx > 0);
pdo_n = find(PDO_idx < 0);

pdo_p_w = intersect(mm,pdo_p);
pdo_n_w = intersect(mm,pdo_n);

PDO_1d_w = PDO_1d(mo_win);

num_PDO_p = length(PDO_1d(PDO_1d_w > 0) );
num_PDO_n = length(PDO_1d(PDO_1d_w < 0) );

figure(11)
set(gca,'FontSize',30)
[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(pdo_p_w),lat_all(pdo_p_w),xbins,ybins);
hold on
set(gcf,'color','w')
title('Positive PDO winter total' )
hold on
set(gcf,'color','w')
set(gca,'FontSize',30)
contourf(X,Y,N'./num_PDO_p,'Linestyle','none')
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
colorbar
caxis([0 2])
box on
    

press_pdo_p_w = press_all(pdo_p_w);
press_bins_pdo_p_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_pdo_p_w(l,m) =  mean(press_pdo_p_w(aa),'omitnan');
        end

    end
end

figure(12)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('Positive PDO Ave Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_pdo_p_w',levs,'Linestyle','none');
h = colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ]) 



u_bins_pdo_p_w = NaN(length(xbins),length(ybins));
v_bins_pdo_p_w =  NaN(length(xbins),length(ybins));
vmag_bins_pdo_p_w =  NaN(length(xbins),length(ybins));

u_pdo_p_w = u_all(pdo_p_w);
v_pdo_p_w = v_all(pdo_p_w);
Vmag_pdo_p_w = Vmag_all(pdo_p_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_pdo_p_w(l,m) =  mean(u_pdo_p_w(aa),'omitnan');
            v_bins_pdo_p_w(l,m) =  mean(v_pdo_p_w(aa),'omitnan');
            vmag_bins_pdo_p_w(l,m) =  mean(Vmag_pdo_p_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 2;
figure(1212)
hold on
set(gcf,'color','w')
title('Positive PDO winter Avg Speed  ')
set(gca,'FontSize',30)
Levs = 0:.4:16;
contourf(xx,yy,vmag_bins_pdo_p_w',Levs,'Linestyle','none')
hold on
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
h =colorbar;
% colormap(redblue)
set(get(h,'title'),'string','ms^-^1');
caxis([0 16])
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_pdo_p_w(1:skip:end,1:skip:end)',v_bins_pdo_p_w(1:skip:end,1:skip:end)','k','Linewidth',2)


%%%%%%%%%%%%%%%%%%%% PDO Negative Winters Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure(13)
set(gca,'FontSize',30)
[N,Xedges,Yedges, binX, binY] = histcounts2(lon_all(pdo_n_w),lat_all(pdo_n_w),xbins,ybins);
hold on
set(gcf,'color','w')
title('Negative PDO winter total')
hold on
set(gcf,'color','w')
set(gca,'FontSize',30)
contourf(X,Y,N'./num_PDO_n,'Linestyle','none')
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
colorbar
%caxis([0 .55])
box on


press_pdo_n_w = press_all(pdo_n_w);
press_bins_pdo_n_w = NaN(length(xbins),length(ybins));
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        %[tf,loc] = ismember([binX',binY'],binToFind,'row');
        aa = find(  binX == binToFind(1) & binY == binToFind(2));        
        if ~isempty(aa)
            press_bins_pdo_n_w(l,m) =  mean(press_pdo_n_w(aa),'omitnan');
        end

    end
end

figure(14)
levs = 950:2:1020;
hold on
set(gcf,'color','w')
title('Negative PDO Winter Pressure at Centers ')
set(gca,'FontSize',30)
contourf(lon_bins,lat_bins,press_bins_pdo_n_w',levs,'Linestyle','none');
h= colorbar;
set(get(h,'title'),'string','m s^-^1');
caxis([950 1020])
plot(coastlon,coastlat,'k','Linewidth',2)
box on
axis([-110 60 10 90 ])



u_bins_pdo_n_w = NaN(length(xbins),length(ybins));
v_bins_pdo_n_w =  NaN(length(xbins),length(ybins));
vmag_bins_pdo_n_w =  NaN(length(xbins),length(ybins));

u_pdo_n_w = u_all(pdo_n_w);
v_pdo_n_w = v_all(pdo_n_w);
Vmag_pdo_n_w = Vmag_all(pdo_n_w);
    
for l = 1:length(xbins)
    for m = 1:length(ybins)
        binToFind = [l m];
        aa = find(  binX == binToFind(1) & binY == binToFind(2));
  
        if ~isempty(aa)
            u_bins_pdo_n_w(l,m) =  mean(u_pdo_n_w(aa),'omitnan');
            v_bins_pdo_n_w(l,m) =  mean(v_pdo_n_w(aa),'omitnan');
            vmag_bins_pdo_n_w(l,m) =  mean(Vmag_pdo_n_w(aa),'omitnan');
            
        end

    end
end

[xx,yy] = meshgrid(xbins,ybins);
skip = 2;
figure(1414)
hold on
set(gcf,'color','w')
title('Negative PDO winter Avg Speed  ')
set(gca,'FontSize',30)
Levs = 0:.4:16;
contourf(xx,yy,vmag_bins_pdo_n_w',Levs,'Linestyle','none')
hold on
plot(coastlon,coastlat,'k','Linewidth',2)
axis([-110 60 10 90 ])
h =colorbar;
% colormap(redblue)
set(get(h,'title'),'string','ms^-^1');
caxis([0 16])
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),u_bins_pdo_n_w(1:skip:end,1:skip:end)',v_bins_pdo_n_w(1:skip:end,1:skip:end)','k','Linewidth',2)










