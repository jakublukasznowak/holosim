
output='/home/pracownicy/jnowak/holo/resolution/plots';
frm='png';
res=300;

% constants
lambda=532e-9;
zR=69e-3;
px=3.1843e-6;
Dh=2048*px;
k0=3.83;
Nbit=8;

% mean fitted values
w=6.2607e-3;
xR=-0.4550e-3; xC=-sign(xR)*3.25e-3;
yR=0.8615e-3; yC=-sign(yR)*3.25e-3;
rmax=sqrt((xR-xC)^2+(yR-yC)^2);


%% 1st, 2nd and approx 3rd limit
z=(25:0.1:65)'*1e-3; zo=zR-z;
r=(0.2:0.2:0.8)'*Dh/2; % distance from the hologram center in hologram plane

dlim1=1.22*2*lambda/Dh*zo;

[Z,R]=meshgrid(zo,r);
dlim2=1.22*2*lambda/Dh*Z./(1-2/Dh*R);

%betaG=10;
%dlim3_approx=betaG/pi*2*lambda/Dh*zo;

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
plot(z*1e3,dlim1*1e6)
plot(z*1e3,dlim2*1e6)
%plot(z*1e3,dlim3_approx*1e6)
legend(cat(1,'Rayleigh',...
    num2cell(num2str(r*1e3,'r_o=%4.2f mm'),2)),...
    'Location','northeast')
xlabel('z [mm]')
ylabel('Resolution [\mum]')
title('Aperture resolution limits')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,'XLim',[min(z) max(z)]*1e3)
print(f,[output,filesep,'res1-3'],['-d',frm],['-r',num2str(res)]);



%% 4th limit
z=(25:0.1:65)'*1e-3; zo=zR-z;
k=(0.25:0.01:1)'; % length with respect to 1st bessel zero


% n=n(zo)
nlim4=zeros(length(z),1);
fitoptions=optimset('TolX',1e-3);
for cnt1=1:length(z)
    fun4limN=@(nn) abs(sqrt(2*nn+2.5)-sqrt(2*nn+0.5)-2*px*sqrt((1/zo(cnt1)-1/zR)/lambda));
    nlim4(cnt1)=fminsearch(fun4limN,10,fitoptions);
end
zo_check=1./((sqrt(2*nlim4+2.5)-sqrt(2*nlim4+0.5)).^2*lambda/4/px^2+1/zR);
z_check=zR-zo_check;
nlim4_approx=lambda./(1./zo-1/zR)/8/px^2-3/4;

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
%plot(z*1e3,nlim4,z_check*1e3,nlim4,z*1e3,nlim4_approx)
plot(z*1e3,nlim4)
xlabel('z [mm]')
ylabel('Number of resolved fringes')
title('Fringes resolvability')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,'XLim',[min(z) max(z)]*1e3)
print(f,[output,filesep,'res4a'],['-d',frm],['-r',num2str(res)]);


% d=d(zo,k)
dlim4=zeros(length(k),length(z));
fitoptions=optimset('TolX',1e-8);
for cnt1=1:length(k)
    for cnt2=1:length(z)
        fun4limN=@(nn) abs(sqrt(2*nn+2.5)-sqrt(2*nn+0.5)-2*px*sqrt((1/zo(cnt2)-1/zR)/lambda));
        fun4_d2n=@(dd) 0.5*(k0*k(cnt1)/pi/dd)^2*lambda*zo(cnt2)*(1-zo(cnt2)/zR)-0.25;
        fun4limD=@(dd) fun4limN(fun4_d2n(dd));
        dlim4(cnt1,cnt2)=fminsearch(fun4limD,3e-6,fitoptions);
    end
end

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
contourf(z*1e3,k,dlim4*1e6,30)
colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
xlabel('z [mm]')
ylabel('Bessel normalized coordinate \kappa')
title('Fringe frequency resolution limit')
print(f,[output,filesep,'res4bCon'],['-d',frm],['-r',num2str(res)]);

[Zo,K]=meshgrid(zo,k);
dlim4_approx=2/pi*k0*px*K.*(1-Zo/zR);

f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
contourf(z*1e3,k,dlim4_approx*1e6,30)
colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
xlabel('z [mm]')
ylabel('Bessel normalized coordinate \kappa')
title('Fringe frequency resolution limit approx')
print(f,[output,filesep,'res4bConAp'],['-d',frm],['-r',num2str(res)]);


f=figure('Color','white',...
    'PaperUnits','centimeters',...
    'PaperSize',[21 29.7],...
    'PaperPosition',[2.5 2.5 12 10]);
ax=axes('Parent',f,'Color','none',...
    'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
hold on
plot(z*1e3,dlim4(1:10:end,:)'*1e6)
legend([repmat('\kappa=',length(1:10:length(k)),1),num2str(k(1:10:end),'%4.2f\n')],'Location','northwest')
%plot(z*1e3,dlim4_approx(1:10:end,:)'*1e6,'r--')
%plot(z*1e3,dlim4(1:10:end,:)'*1e6)
%legend(num2str(k(1:10:end),'k=%4.2f\n'))
xlabel('z [mm]')
ylabel('Resolution [\mum]')
title('Fringe frequency resolution limit')
set(ax,'XGrid','on','YGrid','on','GridAlpha',0.5,'XLim',[min(z) max(z)]*1e3)
print(f,[output,filesep,'res4bLin'],['-d',frm],['-r',num2str(res)]);



%% 5th limit

% parmeters grid
z=(25:0.5:65)'*1e-3; zo=zR-z;
r=(0:1e-4:rmax)'; % distance from the beam center in hologram plane
n=(1:15)'; % number of fringes to be resolved
k=(0.25:0.05:1)'; % length with respect to 1st bessel zero

%% d=d(r,zo,n)
dlim5a=zeros([length(r) length(z) length(n)]);
for cnt1=1:length(r)
    rt=r(cnt1);
    contrast=(1-exp(-2*rmax^2/w^2))*exp(2*rt^2/w^2)*2^-Nbit;
    
    for cnt2=1:length(z)
        zt=zo(cnt2);
        mo=1/(1-zt/zR);
        magnif=mo^2*zt;
    
        for cnt3=1:length(n)
            nt=n(cnt3);
            
            rp=sqrt((2*nt+0.5)*mo*zt*lambda);
            dlim5a(cnt1,cnt2,cnt3)=pi^2/lambda/2*contrast^2/magnif*rp^3*exp(-2*rp^2/w^2);    
        end
        
    end
    
end

for cnt1=1:length(n)
    f=figure('Color','white',...
        'PaperUnits','centimeters',...
        'PaperSize',[21 29.7],...
        'PaperPosition',[2.5 2.5 12 10]);
    ax=axes('Parent',f,'Color','none',...
        'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
    hold on
    contourf(z*1e3,r*1e3,dlim5a(:,:,cnt1)*1e6,20)
    %caxis([min(dlim5a(:)) max(dlim5a(:))]*1e6)
    colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
    xlabel('z [mm]')
    ylabel('r_o [mm]')
    title(sprintf('Fringe contrast resolution limit for n=%d',round(n(cnt1))))
    print(f,[output,filesep,sprintf('res5a_n%03d',round(n(cnt1)))],['-d',frm],['-r',num2str(res)]);
end

% approx
[Z,R,N]=meshgrid(zo,r,n);
contrast=(1-exp(-2*rmax^2/w^2))*exp(2*R.^2/w^2)*2^-Nbit;
mo=1./(1-Z/zR);
dlim5a_approx=pi^2/2*contrast.^2.*sqrt(lambda*Z./mo).*(2*N+0.5).^1.5;
for cnt1=1:length(n)
    f=figure('Color','white',...
        'PaperUnits','centimeters',...
        'PaperSize',[21 29.7],...
        'PaperPosition',[2.5 2.5 12 10]);
    ax=axes('Parent',f,'Color','none',...
        'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
    hold on
    contourf(z*1e3,r*1e3,dlim5a_approx(:,:,cnt1)*1e6,20)
    %caxis([min(dlim5a(:)) max(dlim5a(:))]*1e6)
    colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
    xlabel('z [mm]')
    ylabel('r_o [mm]')
    title(sprintf('Fringe contrast resolution limit n=%d',round(n(cnt1))))
    print(f,[output,filesep,sprintf('res5aAp_n%03d',round(n(cnt1)))],['-d',frm],['-r',num2str(res)]);
end


%% d=d(r,zo,k)
dlim5b=zeros([length(r) length(z) length(k)]);
for cnt1=1:length(r)
    rt=r(cnt1);
    contrast=(1-exp(-2*rmax^2/w^2))*exp(2*rt^2/w^2)*2^-Nbit;
    
    for cnt2=1:length(z)
        zt=zo(cnt2);
        mo=1/(1-zt/zR);
        magnif=mo^2*zt;
    
        for cnt3=1:length(k)
            kt=k(cnt3);
            
            rp=@(d) lambda*zt/pi/d*k0*kt;
            fun5=@(d) abs(d-pi^2/lambda/2*contrast^2/magnif*rp(d)^3);%*exp(-2*rp(d)^2/w^2));
            
            fitoptions = optimset('TolX',1e-8);
            dlim5b(cnt1,cnt2,cnt3)=fminsearch(fun5,1e-6,fitoptions);            
        end
        
    end
    
end

for cnt1=1:length(k) 
    f=figure('Color','white',...
        'PaperUnits','centimeters',...
        'PaperSize',[21 29.7],...
        'PaperPosition',[2.5 2.5 12 10]);
    ax=axes('Parent',f,'Color','none',...
        'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
    hold on
    contourf(z*1e3,r*1e3,dlim5b(:,:,cnt1)*1e6,20)
    %caxis([min(dlim5a(:)) max(dlim5a(:))]*1e6)
    colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
    xlabel('z [mm]')
    ylabel('r_o [mm]')
    title(['Fringe contrast resolution limit for \kappa=',num2str(k(cnt1),'%4.2f')])
    print(f,[output,filesep,sprintf('res5b_k%03d',round(k(cnt1)*100))],['-d',frm],['-r',num2str(res)]); 
end

% approx
[Z,R,K]=meshgrid(zo,r,k);
contrast=(1-exp(-2*rmax^2/w^2))*exp(2*R.^2/w^2)*2^-Nbit;
mo=1./(1-Z/zR);
dlim5b_approx=(2*pi)^-0.25*sqrt(contrast*lambda.*Z./mo).*(k0*K).^0.75;
for cnt1=1:length(k) 
    f=figure('Color','white',...
        'PaperUnits','centimeters',...
        'PaperSize',[21 29.7],...
        'PaperPosition',[2.5 2.5 12 10]);
    ax=axes('Parent',f,'Color','none',...
        'Position',[0.1 0.1 0.85 0.82],'FontSize',8);
    hold on
    contourf(z*1e3,r*1e3,dlim5b_approx(:,:,cnt1)*1e6,20)
    %caxis([min(dlim5a(:)) max(dlim5a(:))]*1e6)
    colormap jet, c=colorbar; c.Label.String='Resolution [\mum]';
    xlabel('z [mm]')
    ylabel('r_o [mm]')
    title(['Fringe contrast resolution limit for \kappa=',num2str(k(cnt1),'%4.2f')])
    print(f,[output,filesep,sprintf('res5bAp_k%03d',round(k(cnt1)*100))],['-d',frm],['-r',num2str(res)]); 
end
