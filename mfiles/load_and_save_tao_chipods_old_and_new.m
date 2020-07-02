% load_and_save_tao_chipods_old_and_new.m
%
% Chipods that were processed with the old software have different naming
% conventions from chipods that are processed with the chipod_gust software.
% 
% This code with load in both old and new and combine into single datasets
%
% **************************************************************
% **** make sure you create the info file:          ************
% **** e.g. ./chipod_info/chipod_info_at_0_104W.m   ************
% **************************************************************
%
% written by Sally Warner    
    
clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;

%% use IC
% decide whether or not to include IC estimates in the final averaged
% product.
% DON'T TURN THIS OFF. IT IS NOW STANDARD TO INCLUDE IC AND VC IN ALL CASES
% WHERE THEY ARE INCLUDED IN TURB.MAT

ynIC = 1;

%% load in chipodinfo
% You need to create an m-file that defines specifics of 
% the chipod deployments at each location:
% e.g. ./chipod_info/chipod_info_at_0_140W.m

run(['./chipod_info/chipod_info_at_' latstrshortnop '_' lonstrnop '.m']);


%% make dummy variables (10 min averages at all depths)

% define the structures
depths = chipodinfo.depths;
for ii = 1:length(depths)
    depthnames{ii} = ['chi' num2str(depths(ii))];
end

% define a continuous timeseries
for ii = 1:length(depths)
    c.(depthnames{ii}).time = time;
end

% define variables that will go in new structure
if ynIC
    newcases = {'mi11';'mi22';'mm1';'mm2';'pi11';'pi22';'pm1';'pm2';...
        'mi11_ic';'mi22_ic';'mm1_ic';'mm2_ic';'pi11_ic';'pi22_ic';'pm1_ic';'pm2_ic'};
else
    newcases = {'mi11';'mi22';'mm1';'mm2';'pi11';'pi22';'pm1';'pm2'};
end
newvars =    {'dTdz';  'T';  'chi';  'eps';  'Kt';  'Jq';  'N2'};
%%% note newvars must correspond to vars1 and vars2 below!! %%%

% predefine variables with NaNs
for ii = 1:length(depths)
    for vv = 1:length(newvars)
        c.(depthnames{ii}).(newvars{vv}) = NaN*time;
    end
    for jj = 1:length(newcases) 
        c.(depthnames{ii}).(newcases{jj}).time = time;
        for vv = 1:length(newvars)                    
            c.(depthnames{ii}).(newcases{jj}).(newvars{vv}) = NaN*time;
        end
    end
end
        

%% load in old data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% old processing internal and mooring stratifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if chipodinfo.ynoldchipods == 1

for kk = 1:2
    clear cpd
    if kk == 1 & ~strcmp(chipodinfo.filenames.oldproc_i,'none')
        % load 10 min chipod data from 2005-2012 (chipod uses internal stratification)
        filename = chipodinfo.filenames.oldproc_i;
        disp(filename)
        load(filename)
        vars1 =     {'dT1dz'; 'T1'; 'chi1'; 'eps1'; 'Kt1'; 'Jq1'};
        vars2 =     {'dT2dz'; 'T2'; 'chi2'; 'eps2'; 'Kt2'; 'Jq2'};
        foundolddata = 1;
    elseif kk == 2 & ~strcmp(chipodinfo.filenames.oldproc_m,'none')
        % load 10 min chipod data from 2005-2012 (chipod uses mooring stratification)
        filename = chipodinfo.filenames.oldproc_m;
        disp(filename)
        load(filename)
        vars1 =     {'dTdz'; 'T1'; 'chi1'; 'eps1'; 'Kt1'; 'Jq1'; 'N2'};
        vars2 =     {'dTdz'; 'T2'; 'chi2'; 'eps2'; 'Kt2'; 'Jq2'; 'N2'}; 
        foundolddata = 1;
    else
        foundolddata = 0;
    end

    if foundolddata
        
        % change sign of Jq
        for ii = 1:length(cpd)
            cpd(ii).Jq1 = -cpd(ii).Jq1;
            cpd(ii).Jq2 = -cpd(ii).Jq2;   
        end

        % index of each depth in cpd structure for each deployment
        inddeps = chipodinfo.inddeps;

        for dd = 1:length(depths)
            for ii = 1:length(cpd)
                % define the index of newdepths(dd) for cpd(ii)
                clear j
                j = inddeps(dd,ii);

                % make sure an index has been defined, if not, no data at that depth
                if j > 0
                    for vv = 1:length(vars1) % note: length(vars1) is longer for mm than mi
                        clear dummy1 dummy2
                        dummy1 = interp1(cpd(ii).time,cpd(ii).(vars1{vv})(j,:),time,'nearest');
                        dummy2 = interp1(cpd(ii).time,cpd(ii).(vars2{vv})(j,:),time,'nearest');
                        if kk == 1
                            c.(depthnames{dd}).mi11.(newvars{vv})(~isnan(dummy1)) ...
                                = dummy1(~isnan(dummy1));
                            c.(depthnames{dd}).mi22.(newvars{vv})(~isnan(dummy2)) ...
                                = dummy2(~isnan(dummy2));
                        elseif kk == 2
                            c.(depthnames{dd}).mm1.(newvars{vv})(~isnan(dummy1)) ...
                                = dummy1(~isnan(dummy1));  
                            c.(depthnames{dd}).mm2.(newvars{vv})(~isnan(dummy2)) ...
                                = dummy2(~isnan(dummy2)); 
                        end
                    end  
                    if kk == 1
                        % add in N2 variable
                        alpha1 = nanmean(sw_alpha(c.(depthnames{dd}).mi11.T*0+35,...
                            c.(depthnames{dd}).mi11.T,...
                            str2num(depthnames{1}(4:length(depthnames{1})))));
                        alpha2 = nanmean(sw_alpha(c.(depthnames{dd}).mi22.T*0+35,...
                            c.(depthnames{dd}).mi22.T,...
                            str2num(depthnames{1}(4:length(depthnames{1})))));
                        c.(depthnames{dd}).mi11.N2 = 9.81*alpha1*c.(depthnames{dd}).mi11.dTdz;
                        c.(depthnames{dd}).mi22.N2 = 9.81*alpha2*c.(depthnames{dd}).mi22.dTdz;
                    end
                end
            end
        end
    end
end

end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% new processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpl = chipodinfo.dpl;
newdepths = chipodinfo.newdepths;
cpds = chipodinfo.cpds;
               
for dd = 1:length(dpl)
    for nn = 1:length(newdepths)
        unit = cpds(dd,nn);
        depth = newdepths(nn);
        if unit ~= 0
            depthname = ['chi' num2str(depth)];

            % load the correct chipod
            load([chipodinfo.basedir.(dpl{dd}) num2str(unit) '/proc/Turb.mat'])
            disp([chipodinfo.basedir.(dpl{dd}) num2str(unit) '/proc/Turb.mat'])

            for ii = 1:length(newcases)
                newcase = 'none';
                if isfield(Turb,[newcases{ii}])
                    newcase = newcases{ii};
                elseif isfield(Turb,['chi_' newcases{ii}])
                    newcase = ['chi_' newcases{ii}];
                end
                if ~strcmp(newcase,'none')
                    disp(newcase)
                    for vv = 1:length(newvars)
                        indgood = find(~isnan(Turb.(newcase).time));
                        clear dummy                  
                        dummy = interp1(Turb.(newcase).time(indgood),...
                            Turb.(newcase).(newvars{vv})(indgood),time);
                        c.(depthname).(newcases{ii}).(newvars{vv})(~isnan(dummy)) ...
                            = dummy(~isnan(dummy));
                    end
                end
            end
        end
       
    end
end

%% average all of the processing cases

disp('averaging onto hourly, daily and monthly time grids...')

for dd = 1:length(depthnames)
    for vv = 1:length(newvars)
        clear dummy
        for ii = 1:length(newcases)
            dummy(ii,:) = c.(depthnames{dd}).(newcases{ii}).(newvars{vv});
        end
        c.(depthnames{dd}).(newvars{vv}) = nanmean(dummy,1);
%         c.(depthnames{dd}).([newvars{vv} '_std']) = nanstd(dummy,1);
    end
end


%% combine into hourly and daily averages
% dt = 6 for hourly averages of 10 minute data
% dt = 6*24 = 144 for daily averages of 10 minute data

for dd = 1:length(depthnames)
    chr.(depthnames{dd}) = bin_average(c.(depthnames{dd}),6);
    cdy.(depthnames{dd}) = bin_average(c.(depthnames{dd}),144);
    for ii = 1:length(newcases)
        chr.(depthnames{dd}).(newcases{ii}) = ...
            bin_average(c.(depthnames{dd}).(newcases{ii}),6);
        cdy.(depthnames{dd}).(newcases{ii}) = ...
            bin_average(c.(depthnames{dd}).(newcases{ii}),144);
    end
end



%% calculate monthly averages

for dd = 1:length(depthnames)
    cmo.(depthnames{dd}) = monthly_average(cdy.(depthnames{dd}));
    for ii = 1:length(newcases)
        disp([depthnames{dd} '.' newcases{ii}])
        cmo.(depthnames{dd}).(newcases{ii}) = ...
            monthly_average(cdy.(depthnames{dd}).(newcases{ii}));
    end
end


%% create structure with all of the data wrt depth

% combine all 10 min, hourly, daily, and monthly data in 4 strucutres (c,
% chr, cdy, and cmo)

% define time vector 
c.time = c.(depthnames{1}).time;
chr.time = chr.(depthnames{1}).time;
cdy.time = cdy.(depthnames{1}).time;
cmo.time = cmo.(depthnames{1}).time;

% define depth vector
c.depth = depths';
chr.depth = depths';
cdy.depth = depths';
cmo.depth = depths';

% load first depth into strucutre
for vv = 1:length(newvars)
    c.(newvars{vv}) = c.(depthnames{1}).(newvars{vv});
    chr.(newvars{vv}) = chr.(depthnames{1}).(newvars{vv});
    cdy.(newvars{vv}) = cdy.(depthnames{1}).(newvars{vv});
    cmo.(newvars{vv}) = cmo.(depthnames{1}).(newvars{vv});
end

% load subsequent depths into structure (start at index 2)
for dd = 2:length(depthnames) 
    for vv = 1:length(newvars)
        c.(newvars{vv}) = [c.(newvars{vv}); c.(depthnames{dd}).(newvars{vv})];
        chr.(newvars{vv}) = [chr.(newvars{vv}); chr.(depthnames{dd}).(newvars{vv})];
        cdy.(newvars{vv}) = [cdy.(newvars{vv}); cdy.(depthnames{dd}).(newvars{vv})];
        cmo.(newvars{vv}) = [cmo.(newvars{vv}); cmo.(depthnames{dd}).(newvars{vv})];   
    end
end



%% save data

disp('saving...')

save([savedir 'chipods_' latstrshort '_' lonstr '_10min.mat'],'c','-v7.3')
save([savedir 'chipods_' latstrshort '_' lonstr '_hourly.mat'],'chr')
save([savedir 'chipods_' latstrshort '_' lonstr '_daily.mat'],'cdy')
save([savedir 'chipods_' latstrshort '_' lonstr '_monthly.mat'],'cmo')


return

%% make plots

set_location;

dd = dir([savedir 'chipod*']);

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-5:end),'y.mat')
        load([savedir dd(ii).name])
    end
end

vars = {'T';'dTdz';'N2';'chi';'eps';'Kt';'Jq'};


% %%%%%%%%%%%%%% plot all variabiles at each depth %%%%%%%%%%%%%%
% 
% for jj = 1:length(c.depth)
% 
%     figure(265+jj)
%     clf
%         set(gcf,'position',[495          57         766        1048])
% 
% 
%     for ii = 1:length(vars)
% 
%         ax(ii) = subplot(length(vars),1,ii);
%         
%         if strcmp(vars{ii},'chi') | strcmp(vars{ii},'eps') | strcmp(vars{ii},'Kt')
%             
%             semilogy(c.time,c.(vars{ii})(jj,:))
%             hold on
%             semilogy(chr.time,chr.(vars{ii})(jj,:))
%             semilogy(cdy.time,cdy.(vars{ii})(jj,:))
%             semilogy(cmo.time,cmo.(vars{ii})(jj,:))
%         
%         else
% 
%             plot(c.time,c.(vars{ii})(jj,:))
%             hold on
%             plot(chr.time,chr.(vars{ii})(jj,:))
%             plot(cdy.time,cdy.(vars{ii})(jj,:))
%             plot(cmo.time,cmo.(vars{ii})(jj,:))
%         end
% 
%         datetick
% 
%         ylabel(vars{ii})
%         set(gca,'tickdir','out')
% 
%         if ii == 1
%             title(['Chipods | ' num2str(c.depth(jj)) 'm | ' latstrshort ', ' lonstr])
%             legend('10min','hr','daily','monthly','location','northwest')
%             legend boxoff
%         end
%     end
% 
%     export_fig([figdir 'chi_' num2str(c.depth(jj)) 'm.png'],'-r200')
%     
% end




%%%%%%%%%%%%%% plot just epsilon at all depths in one plot %%%%%%%%%%%%%%

figure(674)
clf
set(gcf,'position',[495          57         766        1048])

for ii = 1:length(c.depth)

    ax(ii) = subplot(length(c.depth),1,ii);

    semilogy(cdy.time,cdy.eps(ii,:))
    hold on
    semilogy(cmo.time,cmo.eps(ii,:))

    ylim(10.^[-10 -5])
    set(gca,'ytick',10.^(-10:-4))
    datetick

    ylabel({[num2str(c.depth(ii)) 'm'];'\epsilon [W/kg]'})
    set(gca,'tickdir','out')

    if ii == 1
        title(['Chipods   |   \epsilon   |   all depths   |   ' latstrshort ', ' lonstr])
    end
    if ii == length(c.depth)
        legend('daily','monthly','location','southwest')
        legend boxoff
    end
    
    linkaxes
    grid on
end

    export_fig([figdir 'chi_all_depths_just_eps.png'],'-r200')
    

    
    
%%%%%%%%%%%%%% plot just Jq at all depths in one plot %%%%%%%%%%%%%%

figure(675)
clf
set(gcf,'position',[495          57         766        1048])

for ii = 1:length(c.depth)

    ax(ii) = subplot(length(c.depth),1,ii);

    plot(cdy.time,cdy.Jq(ii,:))
    hold on
    plot(cmo.time,cmo.Jq(ii,:))
    
    ylim([-400 0])
    set(gca,'ytick',-500:100:0)
    datetick

    ylabel({[num2str(c.depth(ii)) 'm'];'Jq [W/m^2]'})
    set(gca,'tickdir','out')

    if ii == 1
        title(['Chipods   |   Jq   |   all depths   |   ' latstrshort ', ' lonstr])
    end
    if ii == length(c.depth)
        legend('daily','monthlly','location','northwest')
        legend boxoff
    end
    
    linkaxes
    grid on
end

    export_fig([figdir 'chi_all_depths_just_Jq.png'],'-r200')
    
    
    
    
%%%%%%%%%%%%%% plot just Kt at all depths in one plot %%%%%%%%%%%%%%

figure(676)
clf
set(gcf,'position',[495          57         766        1048])

for ii = 1:length(c.depth)

    ax(ii) = subplot(length(c.depth),1,ii);

    semilogy(cdy.time,cdy.Kt(ii,:))
    hold on
    semilogy(cmo.time,cmo.Kt(ii,:))

    ylim(10.^[-7 -1])
    set(gca,'ytick',10.^(-10:-1))
    datetick

    ylabel({[num2str(c.depth(ii)) 'm'];'K_T [m^2 s^{-1}]'})
    set(gca,'tickdir','out')

    if ii == 1
        title(['Chipods   |   K_T   |   all depths   |   ' latstrshort ', ' lonstr])
    end
    if ii == length(c.depth)
        legend('daily','monthly','location','northwest')
        legend boxoff
    end
    
    linkaxes
    grid on
end

    export_fig([figdir 'chi_all_depths_just_Kt.png'],'-r200')
