% Jakub Nowak 2017 01 04

function [i,iR,iS,Q,argSin] = FraunhoferI (z0,ksi0,eta0,diam,zR,xR,yR,Nx,dx,Ny,dy,outputPath)

if nargin<10
    Ny=Nx;
    dy=dx;
end


% beam
w0=1.75e-6;
B=1;
lambda=532e-9;

Rayleigh=pi*w0^2/lambda;
q1=zR-z0+1j*Rayleigh;
q2=zR+1j*Rayleigh;
m0=abs(q2/q1);

wR=w0*sqrt(1+zR^2/Rayleigh^2);

% object
a=diam/2;

% grid
xv=(-Nx/2+0.5:Nx/2-0.5)'*dx-xR;
yv=(-Ny/2+0.5:Ny/2-0.5)'*dy-yR;

x0=abs(q2/q1)*ksi0; y0=abs(q2/q1)*eta0;

[X,Y]=meshgrid(xv,yv); R=sqrt(X.^2+Y.^2);
Xp=X-x0; Yp=Y-y0; Rp=sqrt(Xp.^2+Yp.^2);


%% intensity direct calculation

% reference wave
iR=(B/Rayleigh*w0/wR*exp(-R.^2/wR^2)).^2;

% scattering
argBessel=2*pi*a*Rp/lambda/abs(z0);
Af=pi*a^2*(2*besselj(1,argBessel))./argBessel; % aperture Fourier transform
Q=m0/lambda/abs(z0)*exp(Rp.^2/wR^2).*Af;
%Q=m0/lambda/abs(z0).*Af;
argSin=pi*Rp.^2/m0/lambda/abs(z0);

iS=iR.*(2*Q.*sin(argSin)-Q.*Q);

% interference pattern
i=iR-iS;


%% save to file
if exist('outputPath','var')
    fileName=[outputPath,filesep,sprintf('d%02d_zh%02d_x%02d_y%02d',...
        round(diam*1e6),round(abs(z0)*1e3),round(ksi0*1e4),round(eta0*1e4))];
    save(fileName,'iS','Q','argSin','z0','ksi0','eta0','diam',...
        'zR','xR','yR','lambda','w0','q1','q2','m0',...
        'x0','y0','xv','yv','Nx','Ny','dx','dy')
end
end