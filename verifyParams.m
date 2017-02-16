% JN 2017 02 14

function verifyParams (sim,mes,zR,plotpath)

sim.diam=round(sim.diam,9);
sim.z=round(sim.pos(:,1),5);

sim.xy=unwrap2(sim.pos(:,2:3),sim.pos(:,2),sim.pos(:,3));
sim.xyc=round(sim.xy.*repmat(zR./(zR-sim.z),1,2),5);
sim.xy=round(sim.xy,5);

sim.x=sim.xy(:,1); sim.y=sim.xy(:,2);
sim.xc=sim.xyc(:,1); sim.yc=sim.xyc(:,2);

sim.vxy=unwrap2(sim.vel(:,2:3),sim.pos(:,2),sim.pos(:,3));
sim.v=round(sqrt(sim.vxy(:,1).^2+sim.vxy(:,2).^2),5);
sim.ang=round(atan(sim.vxy(:,2)./(sim.vxy(:,1)+eps)),5);
sim.vxy=round(sim.vxy,5);
sim.ang(sim.vxy(:,1)==0)=pi/2;
sim.ang(all(sim.vxy==0,2))=0;


mes.z=mes.pos(:,1);
mes.xy=unwrap2(mes.pos(:,2:3),sim.pos(:,2),sim.pos(:,3));
mes.xyc=mes.xy.*repmat(zR./(zR-mes.z),1,2);

mes.x=mes.xy(:,1); mes.y=mes.xy(:,2);
mes.xc=mes.xyc(:,1); mes.yc=mes.xyc(:,2);

mes.vxy=unwrap2(mes.vel(:,2:3),sim.pos(:,2),sim.pos(:,3));
mes.v=sqrt(mes.vxy(:,1).^2+mes.vxy(:,2).^2);
mes.ang=atan(mes.vxy(:,2)./(mes.vxy(:,1)+eps));



vars={'diam','z','x','y','xy','xc','yc','xyc','v','ang'};
xl={'d [\mu m]','z [mm]','x [mm]','y [mm]','r [mm]','x_c [mm]','y_c [mm]','r_c [mm]','v [mm/s]','\phi [rad]'};
factor=[1e6 1e3 1e3 1e3 1e3 1e3 1e3 1e3 1e3 1];

for i=1:numel(vars)
    myerrorstats(sim.(vars{i})*factor(i),mes.(vars{i})*factor(i),...
        ['\delta ',xl{i}],[plotpath,filesep,'err_',vars{i}]);
end


end

function vecN = unwrap2(vec,x,y)

rotMinv=cat(3,[1 0; 0 1],[0 1; -1 0],[-1 0; 0 -1],[0 -1; 1 0]);

L=size(vec,1);
vec=vec'; vecN=nan(size(vec));
for i=1:L
    if x(i)>0 && y(i)>=0
        vecN(:,i)=rotMinv(:,:,1)*vec(:,i);
    elseif x(i)<=0 && y(i)>0
        vecN(:,i)=rotMinv(:,:,2)*vec(:,i);
    elseif x(i)<0 && y(i)<=0
        vecN(:,i)=rotMinv(:,:,3)*vec(:,i);
    elseif x(i)>=0 && y(i)<0
        vecN(:,i)=rotMinv(:,:,4)*vec(:,i);
    else
        error('Unwrap error.')
    end
end
vecN=vecN';
end

function f=myerrorstats(xt,xm,xl,output)

frm='png'; res=300;
nbins=20;

L=size(xt,1); N=size(xt,2);
xtU=unique(xt,'rows'); LU=length(xtU);

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 14 9]);
ax=axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.88 0.88],'FontSize',8);
hold on
indx=~isnan(xm);
y=xm(indx)-xt(indx);
if N>1, y=sqrt(sum(y.^2,2)); end
histogram(ax,y,nbins,'DisplayStyle','stairs','Normalization','pdf','LineWidth',2);
fprintf('%s\n all  mean=%+5.3f  std=%5.3f  meanAbs=%5.3f  stdAbs=%5.3f  maxAbs=%5.3f\n',...
    xl,mean(y),std(y),mean(abs(y)),std(abs(y)),max(abs(y)))
if LU<11
    for i=1:LU
        indx=all([xt==repmat(xtU(i,:),L,1) ~isnan(xm)],2);
        y=xm(indx)-xt(indx);
        if N>1, y=sqrt(sum(y.^2,2)); end
        histogram(ax,y,nbins,...
            'DisplayStyle','stairs','Normalization','pdf','LineWidth',1);
        fprintf('%4.1f  mean=%+5.3f  std=%5.3f  meanAbs=%5.3f  stdAbs=%5.3f  maxAbs=%5.3f\n',...
            sqrt(sum(xtU(i,:).^2,2)),mean(y),std(y),mean(abs(y)),std(abs(y)),max(abs(y)))
    end
    legend(cat(1,'all',num2cell(num2str(sqrt(sum(xtU.^2,2)),'%4.1f'),2)))
end
xlabel(xl)
ylabel('PDF')



if exist('output','var')
    print(f,output,['-d',frm],['-r',num2str(res)]);
end

end