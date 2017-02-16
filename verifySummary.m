
function [mDet,sD] = verifySummary (sim,zR,plotpath)

frm='png'; res=300;
if nargin<3 || isempty(plotpath), ifplot=false; else ifplot=true; end


%% unwrap xy positions
L=length(sim.diam);
x=sim.pos(:,2);
y=sim.pos(:,3);
z=sim.pos(:,1);

rotMinv=cat(3,[1 0; 0 1],[0 1; -1 0],[-1 0; 0 -1],[0 -1; 1 0]);
xy=sim.pos(:,2:3)'; vxy=sim.vel(:,2:3)';
for i=1:L
    if x(i)>0 && y(i)>=0
        xy(:,i)=rotMinv(:,:,1)*xy(:,i);
        vxy(:,i)=rotMinv(:,:,1)*vxy(:,i);
    elseif x(i)<=0 && y(i)>0
        xy(:,i)=rotMinv(:,:,2)*xy(:,i);
        vxy(:,i)=rotMinv(:,:,2)*vxy(:,i);
    elseif x(i)<0 && y(i)<=0
        xy(:,i)=rotMinv(:,:,3)*xy(:,i);
        vxy(:,i)=rotMinv(:,:,3)*vxy(:,i);
    elseif x(i)>=0 && y(i)<0
        xy(:,i)=rotMinv(:,:,4)*xy(:,i);
        vxy(:,i)=rotMinv(:,:,4)*vxy(:,i);
    else
        error('Dupa dupa.')
    end
end
xy=xy'; vxy=vxy';


%% unique values
xyc=round(xy.*[zR./(zR-z) zR./(zR-z)],5);
z=round(sim.pos(:,1),5);
d=round(sim.diam,9);
vx=vxy(:,1); vy=vxy(:,2);
v=round(sqrt(vx.^2+vy.^2),5);
ang=round(atan(vy./(vx+eps)),5);
ang(round(vx,5)==0)=pi/2;
ang(all([round(vx,5)==0 round(vy,5)==0],2))=0;

xycU=unique(xyc,'rows');
zU=unique(z);
dU=unique(d);
vU=unique(v);
angU=unique(ang);


%% translate list to detection matrix
LL=[length(xycU) length(zU) length(dU) length(vU) length(angU)];
mTrans=nan(LL); mDet=nan(LL);
for cnt1=1:LL(1)
    for cnt2=1:LL(2)
        for cnt3=1:LL(3)
            for cnt4=1:LL(4)
                for cnt5=1:LL(5)
                    ind=find(all([xyc(:,1)==xycU(cnt1,1) xyc(:,2)==xycU(cnt1,2) z==zU(cnt2) d==dU(cnt3) v==vU(cnt4) ang==angU(cnt5)],2));
                    if length(ind)==1
                        mTrans(cnt1,cnt2,cnt3,cnt4,cnt5)=ind;
                        mDet(cnt1,cnt2,cnt3,cnt4,cnt5)=sim.indNubbly(ind);
                    elseif length(ind)>1
                        error('More than 1 index found.')
                    end
                end
            end
        end
    end
end

sD=struct('xyc',xycU,'z',zU,'diam',dU,'v',vU,'phi',angU,'detect',mDet);


%% detection report
disp('XY detection stats')
xyDet=zeros(LL(1),1); xyAll=zeros(LL(1),1);
for i=1:LL(1)
    temp=mDet(i,:,:,:,:);
    xyDet(i)=sum(temp(:)>0);
    xyAll(i)=sum(~isnan(temp(:)));
    fprintf('xc=%4.2f yc=%4.2f  %4d/%4d (%5.2f%%) \n',...
        xycU(i,1)*1e3,xycU(i,2)*1e3,xyDet(i),xyAll(i),xyDet(i)/xyAll(i)*100);
end
sD.xyDet=xyDet;

disp('Z detection stats')
zDet=zeros(LL(2),1); zAll=zeros(LL(2),1);
for i=1:LL(2)
    temp=mDet(:,i,:,:,:);
    zDet(i)=sum(temp(:)>0);
    zAll(i)=sum(~isnan(temp(:)));
    fprintf('z=%4.1f  %4d/%4d (%5.2f%%) \n',...
        zU(i)*1e3,zDet(i),zAll(i),zDet(i)/zAll(i)*100);
end
sD.zDet=zDet;

disp('D detection stats')
dDet=zeros(LL(3),1); dAll=zeros(LL(3),1);
for i=1:LL(3)
    temp=mDet(:,:,i,:,:);
    dDet(i)=sum(temp(:)>0);
    dAll(i)=sum(~isnan(temp(:)));
    fprintf('d=%2.0f  %4d/%4d (%5.2f%%) \n',...
        dU(i)*1e6,dDet(i),dAll(i),dDet(i)/dAll(i)*100);
end
sD.dDet=dDet;

disp('Summed detection stats')
fprintf('all %4d/%4d (%5.2f%%) \n',sum(xyDet),sum(xyAll),sum(xyDet)/sum(xyAll)*100);


%% plot series
if ifplot
    rcU=sqrt(xycU(:,1).^2+xycU(:,2).^2);
    [Z,R,D]=meshgrid((zR-zU)*1e3,rcU*1e3,dU*1e6);
    xcU=unique(xycU(:,1)); LX=length(xcU);
    indX=false(LL(1),LX);
    for cntX=1:LX
        indX(:,cntX)=(xycU(:,1)==xcU(cntX));
    end
   
    for cnt5=1:LL(5)
        for cnt4=1:LL(4)
            
            f=figure('Color','white',...
                'PaperUnits','centimeters',...
                'PaperSize',[21 29.7],...
                'PaperPosition',[2.5 2.5 12 10]);
            ax=axes('Parent',f,'Color','none','Box','on',...
                'Position',[0.1 0.09 0.88 0.88],'FontSize',8);
            hold on
            co=get(ax,'ColorOrder');
            
            for cntX=1:LX
                indDet=(mDet(indX(:,cntX),:,:,cnt4,cnt5)>0); indDet=indDet(:);
                Rp=R(indX(:,cntX),:,:); Rp=Rp(:);
                Zp=Z(indX(:,cntX),:,:); Zp=Zp(:);
                Dp=D(indX(:,cntX),:,:); Dp=Dp(:);
                
                scatter(Zp(indDet)+0.15*(Dp(indDet)-15),...
                    Rp(indDet)+0.01*(Dp(indDet)-15),...
                    Dp(indDet),repmat(co(cntX,:),sum(indDet),1),'filled')
                scatter(Zp(~indDet)+0.15*(Dp(~indDet)-15),...
                    Rp(~indDet)+0.01*(Dp(~indDet)-15),...
                    Dp(~indDet),repmat(co(cntX,:),sum(~indDet),1))    
            end
                
            xlabel('z [mm]'), ylabel('r_o [mm]')
            set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,...
                'XMinorGrid','on','YMinorGrid','on','FontSize',8,...
                'XLim',[15 70],'YLim',[0 4.5])
            print(f,[plotpath,filesep,'v',num2str(round(vU(cnt4)*1e3),'%02d'),...
                '_ang',num2str(round(angU(cnt5)*10),'%02d')],...
                ['-d',frm],['-r',num2str(res)]);
                       
        end
    end
    
end
   

%% main plot
rcU=sqrt(xycU(:,1).^2+xycU(:,2).^2);
[Z,R,D]=meshgrid((zR-zU)*1e3,rcU*1e3,dU*1e6);
xyzdDet=sum(sum(mDet>0,5),4);
xyzdAll=sum(sum(~isnan(mDet),5),4);

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none','Box','on',...
    'Position',[0.1 0.09 0.86 0.88],'FontSize',8);
hold on
scatter(Z(:)+0.15*(D(:)-15),R(:)+0.01*(D(:)-15),D(:),xyzdDet(:)./xyzdAll(:)*100,'filled')
colormap jet, c=colorbar('FontSize',8); c.Label.String='% detected'; c.Label.FontSize=8;
xlabel('z [mm]'), ylabel('r_o [mm]')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,...
    'XMinorGrid','on','YMinorGrid','on','FontSize',8,...
    'XLim',[15 70],'YLim',[0 4.5])
if ifplot
    print(f,[plotpath,filesep,'all'],['-d',frm],['-r',num2str(res)]);
end
    
    

    
end


