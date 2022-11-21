%%%%%%%%%%%%%%%Extratropical Cyclone Tracker-JGR Oceans %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This Extratropical Cyclone Tracker is that used in:
%
% Lodise, J., 
% Merrifield, S. T., Collins, C., Rogowski, P., Behrens, & J., Terrill,E, 
% (In Review). Global Climatology of Extratropical Cyclones From a New 
% Tracking Approach and Associated Wave Heights from Satellite Radar 
% Altimeter. Journal of Geophysical Research: Oceans. https://doi.org/10.1029/2022JC018925
%
%
%All Input data was downloaded from the Copernicus Climate Change Service
%(C3S) Climate Date Store > https://cds.climate.copernicus.eu/#!/search?text=ERA5&type=dataset
% from the dataset titled, "ERA5 hourly data on single levels from 1959 to
% present" in the form of monthly data files from 1979 up to and including 2020
%
%
%This code is only to serve as an example of the tracking described in the full publication.
%This code demostrations only the tracking perfomred for one month in the the North
%Atlantic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Correspondence: jlodise@ucsd.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all

format long
ncdisp('ERA5_atm_2018_06.nc') % monthly file containing mslp downloaded from C3S

ncFile = ('ERA5_atm_2018_06.nc'); %Must download data from C3S to run storm tracker

lon = ncread(ncFile,'longitude');
lat = ncread(ncFile,'latitude');


mslp =   ncread(ncFile,'msl')./100 ; %mean sea level pressure
mslp_A = mslp(:,lat >  0,:);

%Read in constant fields form ERA5 including lake cover and land mask
%provided by ECMWF 
load('ERA5_masks.mat'); 

lat_N = lat(lat > 0);

% %% %%%%%%%%%%%%%% Defining Atlantic Tracking Region %%%%%%%%%%%%%%%%%%%%%%%%

x_pol_A = [-87.5  -115  -115 40   40    40.5    9.5    -68    -87.5  ]; % North Atlantic region 
y_pol_A = [  15    20    87  87   68    35.95    5       5     15   ];


[lon_E5_A,lat_E5_A] = meshgrid(lon,lat_N); 
in_Atl = inpolygon(lon_E5_A,lat_E5_A,x_pol_A,y_pol_A);
in_A = in_Atl';

%%% Calculating distances between grid points of lat/lon to calucle Laplace
%%% of MSLP
ddistyy = NaN(size(lon_E5_A));
ddistxx = NaN(size(lon_E5_A));
a=6378.137;
b=6356.7523;   %Earth radii at equator and poles.

for i = 1:length(lat_N)
    lat0 = lat_E5_A(i,:)*(pi/180);    %approx lat of the data-set. 
    R=sqrt(((  (a^2).*cos(lat0)).^2+((b^2).*sin(lat0)).^2)./   ((a.*cos(lat0)).^2 + (b.*sin(lat0)).^2)); 
    
    long = lon_E5_A(i,:);
    
    d1=2*pi.*R./360; %length of 1 degree (km) meridional

    ddistxx(i,1:end-1)=diff(long).*d1(1:end-1).*cos(lat0(1:end-1));
    ddistxx(i,end)=ddistxx(i,end-1); 
end

ddistxx(1,:) = 1;
ddistxx(end,:) = ddistxx(end-1,:);


lat00 = lat_E5_A(:,1)*(pi/180);
R=sqrt((((a^2).*cos(lat00)).^2+((b^2).*sin(lat00)).^2)./   ((a.*cos(lat00)).^2 + (b.*sin(lat00)).^2)); 
d1=2*pi.*R./360; %length of 1 degree (km) meridional 
ddistyy(1:end-1,:) = diff(lat_E5_A,1,1).*d1(1:end-1); 
ddistyy(end,:) = ddistyy(end-1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Calculate the discrete Laplace of MSLP %%%%%%%%%%%%%%%%%%%

mslp_Lap = NaN(size(mslp_A));
du2_dx = NaN(size(mslp_A));
du2_dy = NaN(size(mslp_A));

mslpt = mslp_A.*100;

du2_dy(:,2:end-1,:) = (mslpt(:,1:end-2,:) - 2.*mslpt(:,2:end-1,:) + mslpt(:,3:end,:));  
du2_dx(2:end-1,:,:) = (mslpt(1:end-2,:,:) - 2.*mslpt(2:end-1,:,:) + mslpt(3:end,:,:));

du2_dx(1,:,:) = du2_dx(2,:,:);
du2_dx(end,:,:) = du2_dx(end-1,:,:);
du2_dy(:,1,:) = du2_dy(:,2,:);
du2_dy(:,end,:) = du2_dy(:,end-1,:);

mslp_Lap_CL = du2_dx./(ddistxx'.^2)  + du2_dy./(ddistyy'.^2);
mslp_Lap_CL(:,1:2,:) = NaN;
mslp_Lap_CL(:,end-1:end,:) = NaN;

%%%%%%%%%%%%% Apply Gaussian filter to Laplace(MSLP) field)%%%%%%%%%%%%%
for ll = 1: length(mslp(1,1,:))
    mslp_Lap(:,:,ll) = imgaussfilt(mslp_Lap_CL(:,:,ll),5);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% define time from .nc file%%%%%%%%%%%%%%%%%%%%
time_nc = double(ncread(ncFile,'time'))./24.0  +  datenum( '1900-01-01 00:00:00.0') ;
time = datestr(time_nc,'mm-dd-yyyy HH:MM:SS');  % time in UTC!!!
%time_new = cat(1,time2con,time);

month = str2double(time(1,1:2));
year = str2double(time(1,7:10));

  
% %%%%%%%%%%%%%Track storms in the Atlantic%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This for loop turns feilds into grayscale images and defines lattidue and longitudes of potetnial storm
%centers. 
%Storm centers are found by grouping pixels using bwconncomp and then
%assinging the center as the minimum pressure within each group of pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_lons = NaN(length(time),400);
min_lats = NaN(length(time),400); 
se = strel('disk',2,4); %Morphological structuring element  
for k = 1:length(time)
    test_mslp_A = mslp_A(:,:,k);
    test_mslp_A(in_A==0) = NaN;
    
    %convert feilds of MSLP to inverted Grayscale Image%%%%%%%%%%%%%%%%
    test_mslp_I = -(test_mslp_A);
    I=mat2gray(test_mslp_I,[-1030 max(max(test_mslp_I,[],'omitnan'))])';

    %convert feilds of Laplace(MSLP) to Grayscale Image%%%%%%%%%%%%%%%%
    mslp_L_gf = mslp_Lap(:,:,k);
    mslp_L_gf(in_A==0) = NaN;
    I_lap = mat2gray(mslp_L_gf,[.02 .038])'; %this is in Pa/ km ^2

    %%%%%%%%%%%Erode, Reconstruct, and find regional Max in both images%%%%%%%%%%%%
    Ie0 = imerode(I,se);
    Ire0 = imreconstruct(Ie0,I);
    fgm0 = imregionalmax(Ire0);

    Ie = imerode(I_lap,se);
    Ire = imreconstruct(Ie,I_lap);
    fgm = imregionalmax(Ire);
  
    %%%%%%%%%%%%%%Combine maximums in laplace and mslp fields in one image %%%%%%%%%%%%%%
    fgm(fgm0 == 1) = 1;
         
    %set values outside region of interest to 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fgm(in_A' == 0) = 0;
    
    
    %%%%%Group points within the same identified pressure minimum regions %%%%%%%
    % bwconncomp converts areas identified in grayscale images to pixel
    % groups
    CC = bwconncomp(fgm,8);
    stats = regionprops('struct',CC,'all');

    pixels_2D_all = cat(1,stats.PixelList);
    pixels_2D = regionprops('struct',CC,'PixelList');
    centroid_2D = regionprops('struct',CC,'Centroid');
    if ~isempty(pixels_2D_all) 
        stm_pix = NaN(length(pixels_2D),length(pixels_2D_all(:,1)),length(pixels_2D_all(1,:)));
        stm_mslp = NaN(length(pixels_2D),length(pixels_2D_all(:,1)));
        min_pixels = NaN(length(pixels_2D),2);
        
        %%%%%%%%%%%%%Assign MSLP values to pixel locations%%%%%%%%%
        
        for i = 1:length(pixels_2D)
            pixs = pixels_2D(i).PixelList;
            stm_pix(i,1:length(pixs(:,1)),1:length(pixs(1,:))) =  pixs;

            for j = 1:length(pixs(:,1))
                stm_mslp(i,j) = test_mslp_A(stm_pix(i,j,1),stm_pix(i,j,2) );
            end

            %%%%%%%%% Find minimum pressure location%%%%%%%%%% 
            min_p = min(stm_mslp(i,:),[],'omitnan');
            aa = find(stm_mslp(i,:) == min_p);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Setting center to minimum pressure location
            %If more than one minimum is found within group, the pixel
            %closer to the centroid of the grouping is chosen. 
            %If they are equally distant, one is chosen at random 
            
            if length(aa) > 1
                troid = centroid_2D(i).Centroid;
                dist = sqrt((squeeze(stm_pix(i,aa,1)) - troid(1)).^2   +  (squeeze(stm_pix(i,aa,2)) - troid(2)).^2);
                min_2_troid = find(dist == min(dist));
                aa = aa(min_2_troid);
                aa = aa(1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            min_pixels(i,1:2) = [stm_pix(i,aa,1) stm_pix(i,aa,2)];
        end

        %%%%%%%%%%%%%Record minimum pressure location in lat and lon into arrays%%%%%%%%%%%%%%%%%%%
        
        for h = 1 : length(min_pixels(:,1))
            min_lons(k,h) = lon_E5_A(1,min_pixels(h,1));
            min_lats(k,h) = lat_E5_A(min_pixels(h,2),1); 

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This next section links suspected storm cetners in space and time,
%followed by exlcuding non-cyclone related pressure pertubations from the dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set up arrays to by filled with Storm Variables 

All_Storms_lon_A = NaN(1000,length(time)); %number of storms by time matrix to store storm tracks
All_Storms_lat_A = NaN(1000,length(time));
All_Storms_mslp_A = NaN(1000,length(time));
All_Storms_mslp_Lap_A = NaN(1000,length(time));
All_Storms_LSM_A = NaN(1000,length(time));
All_Storms_LC_A = NaN(1000,length(time));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1; %Start at time = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign first time step of identified pressure centers to new arrays 
%%%% Along with other variables assigned to each pressure center

All_Storms_lon_A(1:sum(~isnan(min_lons(t,:))),t) = min_lons(t,~isnan(min_lons(t,:)));
All_Storms_lat_A(1:sum(~isnan(min_lats(t,:))),t) = min_lats(t,~isnan(min_lats(t,:)));

s = sum(~isnan(All_Storms_lat_A(:,1)));  %initiate storm counter 


for p = 1:sum(~isnan(min_lons(t,:)))
    All_Storms_mslp_A(p,t) = mslp_A(lon == All_Storms_lon_A(p,t),lat_N == All_Storms_lat_A(p,t),t);
    All_Storms_mslp_Lap_A(p,t) = mslp_Lap(lon == All_Storms_lon_A(p,t),lat_N == All_Storms_lat_A(p,t),t);
    All_Storms_LSM_A(p,t) = lsm_A(lon == All_Storms_lon_A(p,t),lat_N == All_Storms_lat_A(p,t));
    All_Storms_LC_A(p,t) = lc_A(lon == All_Storms_lon_A(p,t),lat_N == All_Storms_lat_A(p,t));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Start at Second time step and begin linking centers in time based
%%%%%%%%% on criteria that each congruoys 
for t = 2:length(time) %start loop on second timestep
    clear  stm_ctrs2   stm_ctrs1 
    stm_ctrs2 = [min_lats(t,~isnan(min_lats(t,:))); min_lons(t,~isnan(min_lons(t,:)))]';  %Storm centers at time t
    stm_ctrs1 = NaN(size(stm_ctrs2));
    stm_ctrs1(1:sum(~isnan(min_lats(t-1,:))),:) = [min_lats(t-1,~isnan(min_lats(t-1,:))); min_lons(t-1,~isnan(min_lons(t-1,:)))]'; %Storm centers at time t-1

    for j = 1:length(stm_ctrs2(:,1))

        stm_lon_new = stm_ctrs2(j,2);
        stm_lat_new = stm_ctrs2(j,1);

        Hav_dist =  NaN(length(stm_ctrs1(:,1)),1);
        % Calculate distances between each storm center at
        % previous time step
        for k = 1:length(stm_ctrs1(:,1))  
            latlon2 = stm_ctrs2(j,:);
            latlon1 = stm_ctrs1(k,:);
            [Hav_dist(k),d2km] = lldistkm(latlon1,latlon2);
            if Hav_dist(k) <= 160.00 % If distance between storm cetner of interest and storm cetner at previous time step is less than 160 km, continue 
                storm0_lat = stm_ctrs1(k,1);
                storm0_lon = stm_ctrs1(k,2);

                old_storm_lons = All_Storms_lon_A(:,t-1);
                old_storm_lats = All_Storms_lat_A(:,t-1);

                bb = find(old_storm_lons == storm0_lon & old_storm_lats == storm0_lat);
                if  length(bb) > 1  % if more than one storm center is within 160km at given time we take the storm track with longer tail thus far. Implemented for storms that merge early in lifetime
                    LT = sum(~isnan(All_Storms_lon_A(bb,:)),2) ; 
                    bb = bb(LT == max(LT));
                    bb = bb(1);
                end


                if ~isnan(All_Storms_lon_A(bb,t))  %If storm track has already been assigned a point at this time step, test new point to see if it has a lower mslp.. if so replace with current storm center. 
                    if mslp_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t) < All_Storms_mslp_A(bb,t)
                        All_Storms_lon_A(bb,t) = stm_ctrs2(j,2);
                        All_Storms_lat_A(bb,t) = stm_ctrs2(j,1);
                        All_Storms_mslp_A(bb,t) = mslp_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t);
                        All_Storms_mslp_Lap_A(bb,t) = mslp_Lap(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t); 
                        All_Storms_LSM_A(bb,t) = lsm_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
                        All_Storms_LC_A(bb,t) = lc_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
                    end

                end
                
                if  isnan(All_Storms_lon_A(bb,t)) % If point has not been assigned to storm track yet, assign it 
                    All_Storms_lon_A(bb,t) = stm_ctrs2(j,2);
                    All_Storms_lat_A(bb,t) = stm_ctrs2(j,1);
                    All_Storms_mslp_A(bb,t) = mslp_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t);
                    All_Storms_mslp_Lap_A(bb,t) = mslp_Lap(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t); 
                    All_Storms_LSM_A(bb,t) = lsm_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
                    All_Storms_LC_A(bb,t) = lc_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
                end



            end

        end
        
        if ~any(Hav_dist <= 160.00 ,'all') % if none of the points fall with 160km of previous time steps identified centers, start a new storm track
            s = s+1;
            All_Storms_lon_A(s,t) = stm_ctrs2(j,2);
            All_Storms_lat_A(s,t) = stm_ctrs2(j,1);
            All_Storms_mslp_A(s,t) = mslp_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t);
            All_Storms_mslp_Lap_A(s,t) = mslp_Lap(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1),t);
            All_Storms_LSM_A(s,t) = lsm_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
            All_Storms_LC_A(s,t) = lc_A(lon == stm_ctrs2(j,2),lat_N == stm_ctrs2(j,1));
        end
        if s > (length(All_Storms_lon_A(:,1))-100)  % pad end of array with extra NaN values to fill 
            ExtraNaNs = NaN(1000,length(time));
            All_Storms_lon_A = cat(1,All_Storms_lon_A,ExtraNaNs);
            All_Storms_lat_A = cat(1,All_Storms_lat_A,ExtraNaNs);
            All_Storms_mslp_A = cat(1,All_Storms_mslp_A,ExtraNaNs);
            All_Storms_mslp_Lap_A = cat(1,All_Storms_mslp_Lap_A,ExtraNaNs);
            All_Storms_LSM_A = cat(1,All_Storms_LSM_A,ExtraNaNs);
            All_Storms_LC_A =  cat(1,All_Storms_LC_A,ExtraNaNs);
        end


    end
end   


%%%%% Eliminate storms lasting more than 14 days  %%%%%%%%%%%%%%%%%%%%%%

All_Storms_lon_A( (sum(~isnan(All_Storms_lon_A),2)) > (24*14),:) = []; 
All_Storms_lat_A( (sum(~isnan(All_Storms_lat_A),2)) > (24*14),:) = [];
All_Storms_mslp_A( (sum(~isnan(All_Storms_mslp_A),2)) > (24*14),:) = []; 
All_Storms_mslp_Lap_A( (sum(~isnan(All_Storms_mslp_Lap_A),2)) > (24*14),:) = []; 
All_Storms_LSM_A( (sum(~isnan(All_Storms_LSM_A),2)) > (24*14),:) = []; 
All_Storms_LC_A( (sum(~isnan(All_Storms_LC_A),2)) > (24*14),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Eliminate storms lasting less than 1 day %%%%%%%%%%%%%%%%%%%%%%%%%

All_Storms_lon_A( (sum(~isnan(All_Storms_lon_A),2)) < 24,:) = [];
All_Storms_lat_A( (sum(~isnan(All_Storms_lat_A),2)) < 24,:) = [];
All_Storms_mslp_A( (sum(~isnan(All_Storms_mslp_A),2)) < 24,:) = []; 
All_Storms_mslp_Lap_A( (sum(~isnan(All_Storms_mslp_Lap_A),2)) < 24,:) = []; 
All_Storms_LSM_A( (sum(~isnan(All_Storms_LSM_A),2)) < 24,:) = []; 
All_Storms_LC_A( (sum(~isnan(All_Storms_LC_A),2)) < 24,:) = []; 


%%%%%%%%%% Eliminate storms that don't cross into region of interest %%%%%%

del_idx = NaN(1,length(All_Storms_lat_A(:,1)));
for k = 1:length(All_Storms_lon_A(:,1))


    ss = find(All_Storms_lat_A(k,:) >= 25);
    nn = find(All_Storms_lat_A(k,:) <= 85);
    ww = find(All_Storms_lon_A(k,:) >= -80);
    ee = find(All_Storms_lon_A(k,:) <= 10);


    if isempty(ee) || isempty(ss) || isempty(nn)|| isempty(ww) 

        del_idx(k) = 1;

    end
end


All_Storms_lon_A(~isnan(del_idx),:) = [];
All_Storms_lat_A(~isnan(del_idx),:) = [];
All_Storms_mslp_A(~isnan(del_idx),:) = []; 
All_Storms_mslp_Lap_A(~isnan(del_idx),:) = [];
All_Storms_LSM_A(~isnan(del_idx),:) = []; 
All_Storms_LC_A(~isnan(del_idx),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Eliminate storms that don't travel 1000 km from its starting point %%%%

trop_idx = zeros(length(All_Storms_lon_A(:,1)),1);
dist_tot = NaN(length(All_Storms_lon_A(:,1)),1);
for hh = 1:length(All_Storms_lon_A(:,1))
    ss_lon = All_Storms_lon_A(hh,~isnan(All_Storms_lon_A(hh,:)));
    ss_lat = All_Storms_lat_A(hh,~isnan(All_Storms_lat_A(hh,:)));

    dist = NaN(1,length(ss_lon));
    for jj = 2:length(ss_lon)
        [Hav,d2km] = lldistkm( [ss_lat(jj) ss_lon(jj) ],[ss_lat(1) ss_lon(1)]  ) ;
        dist(jj) = Hav;
    end 

    if any(dist >= 1000,'all')
        dist_tot(hh) = 1;
    end

    %%%Label cyclones that form south of 25 N%%%%%%%%%%%%%%%%%%%
    if ss_lat(1) < 25
        trop_idx(hh) = 1; % set trop index to indicate tropical cyclone formation 
    end
end

All_Storms_lon_A(isnan(dist_tot),:) = [];
All_Storms_lat_A(isnan(dist_tot),:) = [];
All_Storms_mslp_A(isnan(dist_tot),:) = []; 
All_Storms_mslp_Lap_A(isnan(dist_tot),:) = []; 
All_Storms_LSM_A(isnan(dist_tot),:) = []; 
All_Storms_LC_A(isnan(dist_tot),:) = []; 
trop_idx(isnan(dist_tot),:) = []; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Apply Savitzky-Golay filtering to storm tracks %%%%%%%%%%%%%%%%%%%%%

All_Storms_lat_filt_A = NaN(size(All_Storms_lat_A));
All_Storms_lon_filt_A = NaN(size(All_Storms_lat_A));
for i = 1:length(All_Storms_lat_A(:,1))
    All_Storms_lat_filt_A(i,~isnan(All_Storms_lat_A(i,:)))  = sgolayfilt(All_Storms_lat_A(i,~isnan(All_Storms_lat_A(i,:))),3,15);
    All_Storms_lon_filt_A(i,~isnan(All_Storms_lon_A(i,:)))  = sgolayfilt(All_Storms_lon_A(i,~isnan(All_Storms_lon_A(i,:))),3,15);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Plot storm tracks and place variables into matlab struct %%%%%%%%%

A_stms = struct;

figure(1)
load coastlines
hold on
set(gcf,'color','w')
set(gca,'FontSize',30)
plot(coastlon,coastlat,'k') 
box on
axis([-150 80 15 90])
title('Comparison of Storm Tracks')
xlabel('Longitude')
ylabel('Latitude')
hold on

for ss = 1:length(All_Storms_lon_filt_A(:,1))

    A_stms(ss).time = datenum(time(~isnan(All_Storms_lon_filt_A(ss,:)),:)); 
    A_stms(ss).lat = All_Storms_lat_filt_A(ss,~isnan(All_Storms_lat_filt_A(ss,:)))';
    A_stms(ss).lon = All_Storms_lon_filt_A(ss,~isnan(All_Storms_lon_filt_A(ss,:)))';
    A_stms(ss).lat_raw = All_Storms_lat_A(ss,~isnan(All_Storms_lat_filt_A(ss,:)))';
    A_stms(ss).lon_raw = All_Storms_lon_A(ss,~isnan(All_Storms_lon_filt_A(ss,:)))';
    A_stms(ss).name = sprintf('%i',ss);
    A_stms(ss).mslp =  All_Storms_mslp_A(ss,~isnan(All_Storms_lon_filt_A(ss,:)))';
    A_stms(ss).mslp_Lap = All_Storms_mslp_Lap_A(ss,~isnan(All_Storms_lon_filt_A(ss,:)))';
    A_stms(ss).trop_idx = trop_idx(ss);

    plot(A_stms(ss).lon,A_stms(ss).lat,'b')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Save data into a monthly file %%%%%%%%%%%%%%%%%
fname = sprintf('Atlantic_ECs_github_ex%0.4i_%0.2i.mat',year,month);

save(fname, 'A_stms')

