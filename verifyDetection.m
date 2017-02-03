% Jakub Nowak 2017 02 02

function [pairedInd,pairedData,falseInd,falseData] = verifyDetection (sim,mes,maxDist,plotoutput)

if nargin<4, plotoutput=''; end

%% get measured pos, diam and vel
zR=abs(sim.beam(1)); xR=sim.beam(2); yR=sim.beam(3);

if ~isfield(mes,'pos')
    zpos=zR*mes.zpos./(zR+mes.zpos);
    xpos=(zR-zpos)/zR.*(mes.xpos-xR)+xR;
    ypos=(zR-zpos)/zR.*(mes.ypos-yR)+yR;
    mes.pos=[zpos xpos ypos];
end
   
if ~isfield(mes,'diam')
    mes.diam=(zR-mes.pos(:,1))/zR.*mes.minsiz;
end

Nsim=length(sim.diam);
Nmes=length(mes.diam);

if ~isfield(mes,'vel')
    ang=-mes.orient*pi/180;
    ang(ang<0)=ang(ang<0)+pi;
    v=(mes.majsiz-mes.minsiz)/sim.expT.*(zR-mes.pos(:,1))/zR;
    mes.vel=[zeros(Nmes,1) v.*cos(ang) v.*sin(ang)];
end

if all(sim.pos(:,1)<0)
    sim.pos(:,1)=-sim.pos(:,1);
end


%% check detection
pairedInd=zeros(Nsim,1);
pairedData.pos=nan(Nsim,3);
pairedData.diam=nan(Nsim,1);
pairedData.vel=nan(Nsim,3);


for i=1:Nsim
    distances=sqrt(sum((repmat(sim.pos(i,:),Nmes,1)-mes.pos).^2,2));
    [closeDist,closeInd]=min(distances);
    if closeDist<maxDist
        pairedInd(i)=closeInd;
        pairedData.pos(i,:)=mes.pos(closeInd,:);
        pairedData.diam(i,:)=mes.diam(closeInd);
        pairedData.vel(i,:)=mes.vel(closeInd,:);
    end
end


%% list of false detections
falseInd=setdiff((1:Nmes)',pairedInd);
falseData.pos=mes.pos(falseInd,:);
falseData.diam=mes.diam(falseInd);
falseData.vel=mes.vel(falseInd,:);


%% plot
if ~isempty(plotoutput)
    dt=4*sim.expT;
    
    indSimFound=(pairedInd>0);
    zt=sim.pos(:,1);
    xt=sim.pos(:,2); yt=sim.pos(:,3);
    dx=sim.vel(:,2)*dt; dy=sim.vel(:,3)*dt;
    myscatter(xt,yt,zR-zt,sim.diam,indSimFound,[plotoutput,'_sim'],dx,dy);
    
    %ztc=sim.pos(:,1)*zR./(zR-sim.pos(:,1));
    ztc=zR-sim.pos(:,1);
    xtc=zR*(sim.pos(:,2)-xR)./(zR-zt)+xR; ytc=zR*(sim.pos(:,3)-yR)./(zR-zt)+yR;
    dxc=sim.vel(:,2)*dt*zR./(zR-zt); dyc=sim.vel(:,3)*dt*zR./(zR-zt);
    myscatter(xtc,ytc,ztc,sim.diam,indSimFound,[plotoutput,'_simC'],dxc,dyc);
   
    
    indMesTrue=~ismember((1:Nmes)',falseInd);
    zt=mes.pos(:,1);
    xt=mes.pos(:,2); yt=mes.pos(:,3);
    dx=mes.vel(:,2)*dt; dy=mes.vel(:,3)*dt;
    myscatter(xt,yt,zR-zt,mes.diam,indMesTrue,[plotoutput,'_mes'],dx,dy);
    
    %ztc=mes.pos(:,1)*zR./(zR-mes.pos(:,1));
    ztc=zR-mes.pos(:,1);
    xtc=zR*(mes.pos(:,2)-xR)./(zR-zt)+xR; ytc=zR*(mes.pos(:,3)-yR)./(zR-zt)+yR;
    dxc=mes.vel(:,2)*dt*zR./(zR-zt); dyc=mes.vel(:,3)*dt*zR./(zR-zt);
    myscatter(xtc,ytc,ztc,mes.diam,indMesTrue,[plotoutput,'_mesC'],dxc,dyc);
    
end


end


function f=myscatter(x,y,z,d,ind,output,dx,dy)
frm='png';
res=300;

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none','Box','on',...
    'Position',[0.075 0.09 0.89 0.9],'FontSize',8);
hold on

scatter(x(ind)*1e3,y(ind)*1e3,10*d(ind)*1e6,z(ind)*1e3,'filled')
scatter(x(~ind)*1e3,y(~ind)*1e3,10*d(~ind)*1e6,z(~ind)*1e3)
colormap jet, c=colorbar('FontSize',8); c.Label.String='z [mm]'; c.Label.FontSize=8;

if nargin>6
    q=quiver(x*1e3,y*1e3,dx*1e3,dy*1e3);
    q.AutoScale='off'; q.Color='black';
    %q1.MaxHeadSize=10; q1.AlignVertexCenters='on';
end

xlabel('x [mm]'), ylabel('y [mm]')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,...
    'XLim',[-3.25 3.25],'YLim',[-3.25 3.25],'FontSize',8,...
    'YDir','reverse')
print(f,output,['-d',frm],['-r',num2str(res)]);

end