% plot_all_data
%
% This function plots the processed data at a given location in the
% TaoTritonPirataRama array.
%
%

clear
addpath('./utilities/')

%% define location
% this creates correct strings and folders

%%%%%%%%% be sure to set the correct location in set_location.m! %%%%%%%%%%

set_location;

%% load common tao time

time = taotime;

%% If desired set time limits that highlight part of the dataset

% If you want to set special time limits, set yncutx = 1 and set special
% time limits in the variable timelims. Set yncutx = 0 to see only the
% timeseries over the entire length of taotime.
% Note, that when special time limits are defined, both figures showing the
% entire time series and the truncated time series will be saved.

yncutx = 1;
if yncutx
    timelims = [datenum(2014,1,1) datenum(2020,1,1)];
end





%% plot vel, adcp, cur

dd = dir([savedir 'vel*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v';'dudz';'dvdz';'Ssq';'S'};


figure(489)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    pcolor(veldy.time,veldy.depth,veldy.(vars{ii}))
    shading flat
    axis ij
    colormap(redblue3)
    colorbar
    if ii == 1 | ii == 2
        caxis([-1 1]*3)
    elseif ii == 3 | ii == 4
        caxis([-1 1]*0.04)
    elseif ii == 5
        caxis([0 2]*10^-3)
    else
        caxis([0 0.04])
    end
        
    ylabel(vars{ii})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'vel.png'],'-r200')


if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'vel_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end


%%%%%%%%%%%%%%%%%%%% ADCP %%%%%%%%%%%%%%%%%%%%

dd = dir([savedir 'adcp*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v';'dudz';'dvdz';'Ssq';'S'};


figure(589)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    pcolor(adcpdy.time,adcpdy.depth,adcpdy.(vars{ii}))
    shading flat
    axis ij
    colormap(redblue3)
    colorbar
    if ii == 1 | ii == 2
        caxis([-1 1]*3)
    elseif ii == 3 | ii == 4
        caxis([-1 1]*0.04)
    elseif ii == 5
        caxis([0 2]*10^-3)
    else
        caxis([0 0.04])
    end
        
    ylabel(vars{ii})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['ADCP Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'adcp.png'],'-r200')


if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'adcp_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end

%%%%%%%%%%%%%%%%%%%% Current meters ****************


dd = dir([savedir 'cur*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-8:end),'daily.mat')
        disp(dd(ii).name)
        load([savedir dd(ii).name])
    end
end

vars = {'u';'v'};


figure(4689)
clf
set(gcf,'position',[376          57        1157        1048])

for ii = 1:length(curdy.depth)
    
    ax(ii) = subplot(length(curdy.depth),1,ii);
    
    plot(curdy.time,curdy.u(ii,:))
    hold on
    plot(curdy.time,curdy.v(ii,:))
    ylim([-2 2])
        
    ylabel({[num2str(curdy.depth(ii)) ' m'];'u,v [m s^{-1}]'})
    datetick
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Current Meter Velocity at ' latstrshort ', ' lonstr])
    end
end

linkaxes

export_fig([figdir 'cur.png'],'-r200')



if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'cur_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end

%% plot temperature and stratification

dd = dir([savedir 'T_*']);

if ~isempty(dd)

vars = {'T';'dTdz';'NTsq'};

if length(dd) == 4
    toproc = [2 3 4];
    yn3met = 0;
elseif length(dd) == 5
    toproc = [3 4 5];
    yn3met = 1;
end

for jj = toproc
    
    clear t
    
    load([savedir dd(jj).name])
    
    if jj == 1
        lab = '10-minute';
    elseif (jj == 2 & yn3met == 0) | (jj == 3 & yn3met == 1) 
        t = tdy;
        lab = 'daily';
    elseif (jj == 3 & yn3met == 0) | (jj == 4 & yn3met == 1) 
        t = thr;
        lab = 'hourly';
    elseif (jj == 4 & yn3met == 0) | (jj == 5 & yn3met == 1) 
        t = tmo;
        lab = 'monthly';
    end
    
    figure(346)
    clf
    set(gcf,'position',[376          57        1157        1048])

    for ii = 1:length(vars)

        ax(ii) = subplot(length(vars),1,ii);

        pcolor(t.time,t.depth,t.(vars{ii}))
        shading flat
        axis ij
        colormap(redblue3)
        cb = colorbar;
        if ii == 1
            caxis([10 30])
            units = '[°C]';
        elseif ii == 2
            caxis([-1 1]*0.2)
            units = '[°C m^{-1}]';
        elseif ii == 3
            caxis([0 1]*10^-3)
            units = '[s^{-2}]';
        end        
        ylabel({vars{ii};units})
        datetick
        set(gca,'tickdir','out')
        if ii == 1
            title(['Temperature at ' latstrshort ', ' lonstr ', ' lab ' temporal resolution'])
        end
    end

    linkaxes

    export_fig([figdir 'T_pcolor_' lab '.png'],'-r200')
    if jj ~= 1
        % Saving the 10 minute data as a figure is slow and nearly 1.5GB
        savefig([figdir 'T_pcolor_' lab '.fig'])
    end
    
    if yncutx
        for ii = 1:length(vars)
            axes(ax(ii));
            xlim(timelims)
            datetick('x','keeplimits')
        end
        export_fig([figdir 'T_pcolor_' datestr(timelims(1),'yyyymmdd') ...
            '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
    end
end



%% plot

load([savedir 'T_' latstrshort '_' lonstr '_daily.mat'])


%%%%%%%%% plot temperature %%%%%%%%%
clear ax

figure(346)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(311);
pcolor(tdy.time,tdy.depth,tdy.T)
shading flat
axis ij
hold on
% plot(tdy.time,tdy.depth20,'k')
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
title(['Temperature and stratification at ' latstrshort ', ' lonstr])
ylabel({'Temperature [°C]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(312);
pcolor(tdy.time,tdy.depth,tdy.dTdz)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*0.25)
ylabel({'dT/dz [°C m^{-1}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(313);
pcolor(tdy.time,tdy.depth,tdy.NTsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N_T^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Temperature.png'],'-r200')

if yncutx
    for ii = 1:3
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'Temperature_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end


%%%%%%%%% plot salinity %%%%%%%%%
clear ax

figure(347)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(311);
pcolor(tdy.time,tdy.depth,tdy.S)
shading flat
axis ij
hold on
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
title(['Salinity and stratification at ' latstrshort ', ' lonstr])
ylabel({'Salinity [psu]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(312);
pcolor(tdy.time,tdy.depth,tdy.dSdz)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*0.25)
ylabel({'dS/dz [psu m^{-1}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(313);
pcolor(tdy.time,tdy.depth,tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*1e-3)
ylabel({'N_S^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Salinity.png'],'-r200')

if yncutx
    for ii = 1:3
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'Salinity_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end



%%%%%%%%% compare NTsq, NSsq, Nsq %%%%%%%%%
clear ax

cax = [-1 1]*4e-4;

figure(348)
clf
set(gcf,'position',[512          57        1117        1048])


ax(1) = subplot(511);
pcolor(tdy.time,tdy.depth,tdy.NTsq)
shading flat
axis ij
hold on
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis(cax)
title(['Compare buoyancy frequencies (NTsq, NSsq, Nsq) at ' latstrshort ', ' lonstr])
ylabel({'NTsq [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(2) = subplot(512);
pcolor(tdy.time,tdy.depth,tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis(cax)
ylabel({'NSsq [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')


ax(3) = subplot(513);
pcolor(tdy.time,tdy.depth,tdy.Nsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis(cax)
ylabel({'N^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(4) = subplot(514);
pcolor(tdy.time,tdy.depth,tdy.NTsq+tdy.NSsq)
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis(cax)
ylabel({'N_T^2+N_S^2 [s^{-2}]';'depth [m]'})
set(gca,'tickdir','out')

ax(5) = subplot(515);
pcolor(tdy.time,tdy.depth,(tdy.NTsq - tdy.Nsq)./tdy.NTsq*100)
% pcolor(tdy.time,tdy.depth,(tdy.NTsq - tdy.Nsq))
shading flat
axis ij
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
colorbar
caxis([-1 1]*100)
% caxis(cax/10)
ylabel({'error [%]';'(N_T^2-N^2)/N_T^2*100';'depth [m]'})
% ylabel({'error';'N_T^2-N^2';'depth [m]'})
set(gca,'tickdir','out')

linkaxes

export_fig([figdir 'Nsq_compare.png'],'-r200')

if yncutx
    for ii = 1:5
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'Nsq_compare_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end


%%%%%%%%% plot mld and z20 %%%%%%%%%
clear ax

figure(349)
clf
set(gcf,'position',[512          57        1117        1048])

ax(1) = axes;
pcolor(tdy.time,tdy.depth,tdy.T)
shading flat
axis ij
hold on
plot(tdy.time,tdy.depth20,'k')
plot(tdy.time,tdy.mld,'k')
colormap(gca,flipud(lbmap(24,'RedBlue')))
datetick
axc = colorbar;
title(['Temperature, MLD & z_{20} at ' latstrshort ', ' lonstr])
title(axc,'temp [°C]')
ylabel({'temp [°C]';'depth [m]'})
set(gca,'tickdir','out')

export_fig([figdir 'T_mld_z20.png'],'-r200')

if yncutx
    for ii = 1:1
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'T_mld_z20_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end



%%%%%%%%% plot mld at different time-spacing %%%%%%%%%
clear ax

dd = dir([savedir 'mld*']);

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(mld);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(mld.time,mld.(vars{ii}))
    hold on
    plot(mldhr.time,mldhr.(vars{ii}))
    plot(mlddy.time,mlddy.(vars{ii}))
    plot(mldmo.time,mldmo.(vars{ii}))
    
    datetick
    
    if ii == 1
        ylabel({'MLD [m]';'depth at which \DeltaT = 0.015°C'})
    elseif ii == 2
        ylabel({'MLD [m] (for velocity interpolation)';'depth at which \DeltaT = 0.04°C'})
    end
%     ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Mixed Layer depth at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'mld.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'mld_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end




%% plot winds
clear ax


dd = dir([savedir 'wind*']);

if ~isempty(dd)

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

vars = {'u';'v';'speed';'dir';'tau';'taux';'tauy'};


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(wind.time,wind.(vars{ii}))
    hold on
    plot(windhr.time,windhr.(vars{ii}))
    plot(winddy.time,winddy.(vars{ii}))
    plot(windmo.time,windmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Winds at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'wind.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'wind_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end


%% plot Jq0
clear ax

dd = dir([savedir 'Jq0*']);

if ~isempty(dd)

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(Jq0);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(Jq0.time,Jq0.(vars{ii}))
    hold on
    plot(Jq0hr.time,Jq0hr.(vars{ii}))
    plot(Jq0dy.time,Jq0dy.(vars{ii}))
    plot(Jq0mo.time,Jq0mo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Jq0 at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'Jq0.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'Jq0_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end


end

%% plot SST
clear ax

dd = dir([savedir 'sst*']);

if ~isempty(dd)

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(sst);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(sst.time,sst.(vars{ii}))
    hold on
    plot(ssthr.time,ssthr.(vars{ii}))
    plot(sstdy.time,sstdy.(vars{ii}))
    plot(sstmo.time,sstmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['SST at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'sst.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'sst_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end


%% plot EUC depth
clear ax

dd = dir([savedir 'euc*']);

if ~isempty(dd)

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(euc);
vars = varsrm([2 4 5]);

figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(euc.time,euc.(vars{ii}))
    hold on
    plot(euchr.time,euchr.(vars{ii}))
    plot(eucdy.time,eucdy.(vars{ii}))
    plot(eucmo.time,eucmo.(vars{ii}))
    
    datetick
    
    if ii == 1
        ylabel({'speed (max(u))';'[m s^{-1}]'})
        ylim([0 2.5])
    elseif ii == 2
        ylabel({'z_{EUC}: depth of EUC';'(i.e. depth of max(u)) [m]'})
    elseif ii == 3
        ylabel({'w = dz_{EUC}/dt';'[m day^-1]'})
    else
        ylabel(vars{ii})
    end
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['Equatorial undercurrent depth at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'euc.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'euc_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end

end


%% plot TIW index
clear ax

dd = dir([savedir 'tiw*']);

if ~isempty(dd)

for ii = 1:length(dd)
    load([savedir dd(ii).name])
end

varsrm = fields(tiw);
vars = varsrm(2:end-1);


figure(265)
clf
    set(gcf,'position',[495          57         766        1048])


for ii = 1:length(vars)
    
    ax(ii) = subplot(length(vars),1,ii);
    
    plot(tiw.time,tiw.(vars{ii}))
    hold on
    plot(tiwhr.time,tiwhr.(vars{ii}))
    plot(tiwdy.time,tiwdy.(vars{ii}))
    plot(tiwmo.time,tiwmo.(vars{ii}))
    
    datetick
    
    ylabel(vars{ii})
    set(gca,'tickdir','out')
    
    if ii == 1
        title(['TIW KE at ' latstrshort ', ' lonstr])
        legend('10min','hr','daily','monthly','location','northwest')
        legend boxoff
    end
end

export_fig([figdir 'tiw.png'],'-r200')

if yncutx
    for ii = 1:length(vars)
        axes(ax(ii));
        xlim(timelims)
        datetick('x','keeplimits')
    end
    export_fig([figdir 'tiw_' datestr(timelims(1),'yyyymmdd') ...
        '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
end


end



%% plot chipods
clear ax

dd = dir([savedir 'chipod*']);

if ~isempty(dd)

for ii = 1:length(dd)
    if strcmp(dd(ii).name(end-5:end),'ly.mat')
        load([savedir dd(ii).name])
    end
end




%%%%%%%%%%%%%% plot just chi at all depths in one plot %%%%%%%%%%%%%%
%%%%%%%%%%%%%% plot just epsilon at all depths in one plot %%%%%%%%%%%%%%
%%%%%%%%%%%%%% plot just Kt at all depths in one plot %%%%%%%%%%%%%%
%%%%%%%%%%%%%% plot just Jq at all depths in one plot %%%%%%%%%%%%%%
clear ax


vars = {'chi';'eps';'Kt';'Jq'};
ylims = {10.^[-9 -4];10.^[-9 -4];10.^[-7 -1];[-1 0]*500};
yticks = {10.^(-9:-4);10.^(-9:-4);10.^(-10:0);-1000:100:0};
labs = {'\chi_T [K^2 s^{-1}]';'\epsilon [W kg^{-1}]';'K_T [m^2 s^{-1}]';'Jq [W m^{-2}]'};

for vv = 1:length(vars)

    figure(674+vv)
    clf
    set(gcf,'position',[495          57         766        1048])

    for ii = 1:length(cmo.depth)

        ax(ii) = subplot(length(cmo.depth),1,ii);

        if strcmp(vars{vv},'eps') | strcmp(vars{vv},'Kt') | strcmp(vars{vv},'chi')
            semilogy(cdy.time,cdy.(vars{vv})(ii,:))
            hold on
            semilogy(cmo.time,cmo.(vars{vv})(ii,:))
        else
            plot(cdy.time,cdy.(vars{vv})(ii,:))
            hold on
            plot(cmo.time,cmo.(vars{vv})(ii,:))
        end

        ylim(ylims{vv})
        set(gca,'ytick',yticks{vv})
        datetick

        ylabel({[num2str(cmo.depth(ii)) 'm'];labs{vv}})
        set(gca,'tickdir','out')

        if ii == 1
            title(['Chipods   |   ' labs{vv} '   |   all depths   |   ' latstrshort ', ' lonstr])
            legend('daily','monthly','location','northwest')
            legend boxoff
        end

        linkaxes
        grid on
    end

    export_fig([figdir 'chi_all_depths_just_' vars{vv} '.png'],'-r200')

    if yncutx
        for ii = 1:length(cmo.depth)
            axes(ax(ii));
            xlim(timelims)
            datetick('x','keeplimits')
        end
        export_fig([figdir 'chi_all_depths_just_' vars{vv} '_' ...
            datestr(timelims(1),'yyyymmdd') ...
            '_to_' datestr(timelims(2),'yyyymmdd') '.png'],'-r200')
    end   
    
end


end