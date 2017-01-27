
% params
prefix='29_';
beamfile='/home/pracownicy/jnowak/holo/beam/beamParams_29.mat';
output='/home/pracownicy/jnowak/holo/beam/plots';


frm='png';
res=300;

lambda=532e-9;
zR=69e-3;


% load file
load(beamfile)
L=length(cR);
time=(1:L)'/15;
w0=zR*lambda/pi./wR;


% prepare to plot
vars={'xR','yR','cR','wR','w0'};
labels={'x_0 [mm]','y_0 [mm]','C [a.u.]','w_R [mm]','w_0 [\mum]'};
factor=[1e3 1e3 1 1e3 1e6];
legStr=num2str(cutoff'*1e3,'lowpass %3.1f mm');
ifleg=[1 0 0 0 0];


% plot
for i=1:numel(vars)
    
    f=figure('Color','white',...
            'PaperUnits','centimeters',...
            'PaperSize',[21 29.7],...
            'PaperPosition',[2.5 2.5 16 7]);
    ax=axes('Parent',f,'Color','none','Box','off',...
        'Position',[0.08 0.13 0.88 0.83],'FontSize',8);
    hold on
    
    plot(time,eval(vars{i})*factor(i))
    xlabel('Time [s]')
    ylabel(labels{i})
    
    if ifleg(i), legend(legStr,'Location','southeast'); end
    set(ax,'XLim',[time(1) time(end)],'XGrid','on','YGrid','on','GridAlpha',0.5)
    
    print(f,[output,filesep,prefix,vars{i}],['-d',frm],['-r',num2str(res)]);
    
end


f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
scatter(xR(:,1)*1e3,yR(:,1)*1e3,3,time,'filled')
colormap jet
c=colorbar;
c.Label.String='Time [s]';
xlabel('x [mm]')
ylabel('y [mm]')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,...
    'XLim',[min(xR(:,1)) max(xR(:,1))]*1e3,'YLim',[min(yR(:,1)) max(yR(:,1))]*1e3)
print(f,[output,filesep,prefix,'xy'],['-d',frm],['-r',num2str(res)]);

set(ax,'XLim',[-3.25 3.25],'YLim',[-3.25 3.25])
print(f,[output,filesep,prefix,'xyEx'],['-d',frm],['-r',num2str(res)]);
