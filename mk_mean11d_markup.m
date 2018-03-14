warning off

%% Initialize variables.
filename = '/media/derek/data/TESS/TDA-3 data/Data_Batch_TDA3_2.txt';
delimiter = ',';
startRow = 14;

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
star_name = dataArray{:, 1};
Tmag = dataArray{:, 2};
cadence = dataArray{:, 3};
tobs = dataArray{:, 4};
eclat = dataArray{:, 5};
eclon = dataArray{:, 6};
teff = dataArray{:, 7};
teff_err = dataArray{:, 8};
logg = dataArray{:, 9};
logg_err = dataArray{:, 10};
star_type = dataArray{:, 11};


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


% we are only interested in files that are sampled the same way...
%calculate distances between points

dist = zeros(numel(star_name),2);

ts_files = dir('Star*.noisy');
num_files = numel({ts_files.name});
star = zeros(num_files,3);
star = struct('name',star_name,'time',[],'flux',[],'fcorr',[]);
for file_list = 1:num_files
    fname = strcat(string(star_name(file_list)),'.noisy');
    test = dlmread(fname,'',5,0);
    temp = find(test(:,3)); %find nonzero flag entries
    test(temp,:) = []; %remove nonzero flag entries
    star(file_list).star = star_name(file_list); %populate structure
    star(file_list).time = test(:,1);
    star(file_list).flux = test(:,2);
    fmean(file_list) = mean(star(file_list).flux);
    fstd(file_list) = std(diff(diff(star(file_list).flux)));
    range(file_list) = prctile(star(file_list).flux,95)-prctile(star(file_list).flux,5);
    range(file_list) = range(file_list)/fmean(file_list);
    range(file_list) = abs(range(file_list));
    drange(file_list) = std(diff(star(file_list).flux));
    drange(file_list) = drange(file_list)/fmean(file_list);
    drange(file_list) = abs(drange(file_list));
    
            sflag = strfind(star_type(file_list),'Eclipse','ForceCellOutput',1);
            sflag = ~cellfun(@isempty,sflag);
            if (sflag==0)
                sflag = strfind(star_type(file_list),'Transit','ForceCellOutput',1);
                sflag = ~cellfun(@isempty,sflag);
            end
            
            if (sflag==0)
                sflag = strfind(star_type(file_list),'LPV','ForceCellOutput',1);
                sflag = ~cellfun(@isempty,sflag);
            end
            
            ttflag(file_list) = sflag;
            if(cadence(file_list)==20)
                ttflag(file_list) = 0;
            end
    
            %this last is only really needed as a rough way to deal with
            %the negative weirdnesses present in some of the current set of light curves
            if min(star(file_list).flux)<0
                ttflag(file_list) = 0;
            end
            
end


%Now go through each star, find comparison targets and detrend

for file_list=1:num_files
    for star_no = 1:num_files
        tdist = sphere_dist(eclat(star_no),eclon(star_no),eclat(file_list),eclon(file_list));
        dist(star_no,1) = star_no;
        dist(star_no,2) = tdist; %note distances are in radians!
    end
   

dist(file_list,2) = 10*pi; %artifically max distance for central target so that it moves to end of list
%check number of nearby stars to determine search radius
sort_dist = sortrows(dist,2);

search_radius = sort_dist(20,2); %20 works well for 20s cadence...more for longer?

time_start = min(star(file_list).time)
time_end = max(star(file_list).time)
full_time = [0 365];
min_range = -2.0;
min_range0 = min_range;
flag = 1;

while 1   

num_star = 0;

% go through the first time and see which has min, max time

full_time = [];
full_flux = [];
full_flag = [];
full_weight = [];
tflux = [];
  
       comp_list = [];
     
        %Need to loop through all the other stars here to make the ensemble
        for test_star = 1:num_files
            if (dist(test_star,2)<search_radius & ttflag(test_star) == 0 & log10(range(test_star)) < min_range & drange(test_star) < 10*drange(file_list)  & min(star(test_star).flux)>0 ) 
                
            %%need to fix to ensure minimum number of stars and that times overlap!
            num_star = num_star+1;
            comp_list = [comp_list;test_star];
            clear test
            test(:,1) = star(test_star).time;
            test(:,2) = star(test_star).flux;  
            test(:,2) = test(:,2)/fmean(test_star); %OK, but need to weight these...
            weight = ones(numel(test(:,1)),1);
            weight = weight*fmean(test_star)/fstd(test_star); %weight based on whitened noise levels
            full_time = [full_time;test(:,1)];
            full_flux = [full_flux;test(:,2).*weight];
      
            full_weight = [full_weight;weight];
            tflux = [tflux;test(:,2)];
            end
        end
     

  %probably need to make the ensemble in pieces as well...could save some time by only looking at pieces that overlap the target star time series...
  
  %%need to figure out a way to check if there's more than N stars (3-4?) in each time segment
  %%we want to detrend and if not, increase the search radius and try
  %%again...
  %how about step through each time step for the star we're working on and
  %check how many points there are? if not >3, break and redo...
 
  gx = [time_start:0.5:time_end]; %0.2 was default, but is this too small?
  [n,edges, bin] = histcounts(full_time,gx);
  %disp(min(n))
  
  if min(n) < 2000 %10000 for 20s cadence?

    min_range = min_range+0.3; %values in this section (0.3, 0.1, 1.2) are a bit rough-and-ready
 
        if min_range > log10(max(range))
            if (search_radius < 0.5)
                search_radius = search_radius+0.1;
            else
                search_radius = search_radius*1.2;
            end
            min_range = min_range0;
        end
     
  
        if search_radius > pi/4 %probably could shrink this for real data
            break
        end
  else
      break
  end
  
end

        
       

full_time(isnan(full_flux)) = [];
full_weight(isnan(full_flux)) = [];
tflux(isnan(full_flux)) = [];
full_flux(isnan(full_flux)) = [];

[full_time,idx] = sort(full_time);
full_flux = full_flux(idx);
full_weight = full_weight(idx);

%Cycle through shrinking binning window, using median binning

%set up temporary files
temp_time = full_time;
temp_flux = full_flux;
temp_weight = full_weight;

btime0 = temp_time;
bflux0 = temp_flux;

temp_time(full_time<time_start | full_time>time_end) = [];
temp_flux(full_time<time_start | full_time>time_end) = [];
temp_weight(full_time<time_start | full_time>time_end) = [];

btime = temp_time;

%Find breaks in ensemble time series

break_locs = int32(find(diff(full_time)>0.1));
     if numel(break_locs)>0
        break_locs = break_locs+1;
        if (break_locs(end) < numel(full_time))
            break_locs = [break_locs;numel(full_time)];
        end
        [bincts2 edges2 bidx2] = histcounts(full_time,full_time(break_locs)); %binning for ensemble
        bidx2 = bidx2+1; 
        num_segs = numel(break_locs);
     else
        [bincts2 edges2 bidx2] = histcounts(full_time,[full_time(1);full_time(end)]); %binning for ensemble
        num_segs = 1;
         break_locs = [1 numel(full_time)];
     end 

    for iseg = 1:num_segs
        influx = full_flux(bidx2==iseg);
        inweight = full_weight(bidx2==iseg);
        intime = full_time(bidx2==iseg);
        
        bin_size = 4.0; %chosen because it seems to work, but somewhat arbitrary

        for ib = 1:4
            clear ttflux ttweight ttime
            bin_size = bin_size/2.0; %should move to end of loop to avoid redundancy
            gx = [min(intime)-0.02:bin_size:max(intime)+0.02]; %sizes of 0.02 given by data minimum cadence
            if numel(gx) < 3
                gx = [min(intime):0.5*(max(intime)-min(intime)):max(intime)];
            end
            [b1 n1 s1] = bindata(intime,influx,gx);
            [bw nw sw] = bindata(intime,inweight,gx);
            if numel(intime) < 2*numel(gx)
                break
            end
            counter = numel(intime);
            while counter > 0
                pp_ensemble(iseg) = pchip(gx,b1./bw); %spline(gx,b);
                diff1 = (influx./inweight)-ppval(pp_ensemble(iseg),intime);
                sdiff = 3*std(diff1);
                counter = numel(diff1(abs(diff1)>sdiff));
                intime(abs(diff1)>sdiff) = [];
                influx(abs(diff1)>sdiff) = [];
                inweight(abs(diff1)>sdiff) = [];
            end
        end
    end
   
    [nn ee bb] = histcounts(star(file_list).time,full_time(break_locs));
    bb = bb+1;
    for iseg = 1: num_segs
        scale = 1.0;
        star(file_list).fcorr(bb==iseg) = scale*star(file_list).flux(bb==iseg)./ppval(pp_ensemble(iseg),star(file_list).time(bb==iseg))   ;
    end
    
btime2 = temp_time; 
 
bin_size = 4.0;
for ib = 1:10
    clear ttflux ttweight ttime
    bin_size = bin_size/2;
    gx = [time_start:bin_size:time_end];
    [n,edges, bin] = histcounts(temp_time,gx);
    if (min(n) < 10 && ib>2) %fix so that if we're on the first pass....
        break
    end
    for ix = 1: numel(n)
        ttweight(ix) = mean(temp_weight(bin==ix));
        ttime(ix) = mean(temp_time(bin==ix));    
        ttflux(ix) = median(temp_flux(bin==ix)./temp_weight(bin==ix));
    end
    ottime = ttime;
    otflux = ttflux;
    pp = spline(ttime,ttflux);
    
    break_locs = find(diff(star(file_list).time)>0.1);

    if numel(break_locs)>0
        break_locs = break_locs+1;
        if (break_locs(end) < numel(star(file_list).time))
            break_locs = [break_locs;numel(star(file_list).time)];
        end
        [bincts edges bidx] = histcounts(star(file_list).time,star(file_list).time(break_locs)); %binning for star
        bidx = bidx+1; %start counting from one
        [bincts2, edges2 bidx2] = histcounts(temp_time,star(file_list).time(break_locs)); %binning for ensemble
        bidx2 = bidx2+1; 
        num_segs = numel(break_locs);
    else
        [bincts edges bidx] = histcounts(star(file_list).time,[star(file_list).time(1);star(file_list).time(end)]); %binning for star
        [bincts2 edges2 bidx2] = histcounts(temp_time,[star(file_list).time(1);star(file_list).time(end)]); %binning for ensemble
        num_segs = 1;
    end
    
    for iseg = 1:num_segs
        %scale(iseg) = mean(temp_flux(bidx==iseg)./temp_weight(bidx==iseg))
        influx = temp_flux(bidx2==iseg);
        inweight = temp_weight(bidx2==iseg);
        intime = temp_time(bidx2==iseg);
        
        influx = star(file_list).flux(bidx==iseg);
        intime = star(file_list).time(bidx==iseg);
        fun = @(x) sum(((influx./median(influx))-x*ppval(pp,intime)).^2);
        scale(iseg) = fminbnd(fun,0.9,1.5);
        
    end
end


%do I want to rescale the segments? maybe only for ones that don't vary
%much?

tscale = scale(bidx);
tscale = tscale(:);
cflux = star(file_list).flux./(tscale.*ppval(pp,star(file_list).time)); 
cflux2 = cflux;

cflux_mean = mean(cflux);

tflag = strfind(star_type(file_list),'LPV','ForceCellOutput',1);
tflag = ~cellfun(@isempty,tflag);
tflag = 0;

frange = prctile(cflux,95)-prctile(cflux,5);
frange = frange/fmean(file_list);
frange = abs(range(file_list));

seg_flag = 0;
for iseg=1:num_segs
    seg_mean(iseg) = mean(cflux(bidx==iseg));
    cflux(bidx==iseg) = cflux(bidx==iseg)*cflux_mean/seg_mean(iseg); 
  %  ok, this
  %  is not useful when there's a very low frequency component....what if
  %  we detrend first? probably not the best way, so need to figure out how
  %  to handle this...if LPV! (or if max freq is low compared to length of
  %  typical segment?
  
  %maybe just check if it's better before or after? just compare standard
  %deviation?
  
  %alternatively, use find to see if there's big break in the average...abs
  %greater than 0.1, maybe? Making this too complicated!! Just examine the
  %means of all the segments and look for big changes. If those exist, then
  %rescale...yes, I think that will work. Try 10% as a change to start?
  %Want to see if one of the segment means changes by more than 10% from
  %the previous one
  if iseg>1
      if abs((seg_mean(iseg)-seg_mean(iseg-1))/cflux_mean)>0.1
          seg_flag = 1;
      end
  end
end


if seg_flag>0
    %disp('set flag')
    cflux = cflux2;
end

star(file_list).fcorr = cflux;
outfile = strcat('toutput/',string(star_name(file_list)),'.noisy_detrend')
out_data = [star(file_list).time,star(file_list).flux,star(file_list).fcorr];
if min(fmean(file_list))>0 & min(star(file_list).flux)>0
    dlmwrite(char(outfile),out_data,'delimiter','\t','precision','%10.8f')
end

seg_flag = 0;
end
    
