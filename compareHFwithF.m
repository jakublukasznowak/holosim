% Jakub Nowak 2017 01 13

HFpath='/home/pracownicy/jnowak/holo/simulations/databank/stillHuygensFresnel';


z=65e-3; x=3e-3; y=2e-3; diam=25e-6;
zR=-69e-3; zh=z+zR;
ksi=x*abs(z/zR); eta=y*abs(z/zR);
%HFfile=[sprintf('d%02d_zh%02d_x%02d_y%02d',...
%        round(diam*1e6),round(abs(zh)*1e3),round(ksi*1e4),round(eta*1e4)),...
%        '.mat'];
HFfile='d15_zh24_x20_y13.mat';

load([HFpath,filesep,HFfile])



[X,Y]=meshgrid(xv-x0,yv-y0);
Rp=sqrt(X.^2+Y.^2);

Rb=1.22*Rb;

[~,indx]=min(abs(xv-x0));
[~,indy]=min(abs(yv-y0));

[~,fieldR,fS]=Fraunhofer(z0,ksi0,eta0,diam,zR,xR,yR,Nx,dx,Ny,dy);
I_R=abs(fieldR).^2;

field_HF=fieldR-fieldS;
field_HF(Rp>=Rb)=nan;
I_HF=abs(field_HF).^2;

field_F=fieldR-fS;
field_F(Rp>=Rb)=nan;
I_F=abs(field_F).^2;



%% distribution comparison
figure
imagesc(xv,yv,I_R)
colormap jet
colorbar
title('Background intensity distribution')

figure
imagesc(xv,yv,I_HF-I_R)
colormap jet
colorbar
title('Background corrected Huygens-Fresnel intensity')

figure
imagesc(xv,yv,I_F-I_R)
colormap jet
colorbar
title('Background corrected Fraunhofer intensity')

figure
imagesc(xv,yv,I_HF-I_F)
colormap jet
colorbar
title('Difference between Huygens-Fresnel and Fraunhofer')


%% crossection comparison
figure
plot(xv,I_HF(indy,:)-I_R(indy,:),xv,I_F(indy,:)-I_R(indy,:))
legend({'HF','F'})
title('X intensity crossection')

figure
plot(yv,I_HF(:,indx)-I_R(:,indx),yv,I_F(:,indx)-I_R(:,indx))
legend({'HF','F'})
title('Y intensity crossection')