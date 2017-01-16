% Jakub Nowak 2017 01 09

function [field,fieldR,fieldS] = HuygensFresnel (z0,ksi0,eta0,diam,zR,xR,yR,Nx,dx,Ny,dy,outputPath)

lambda=532e-9;
w0=1.75e-6;
B=1;

if nargin<10
    Ny=Nx;
    dy=dx;
end


xv=(-Nx/2+0.5:Nx/2-0.5)'*dx-xR;
yv=(-Ny/2+0.5:Ny/2-0.5)'*dy-yR;
Lx=length(xv); Ly=length(yv);

q1=zR-z0+1j*pi*w0^2/lambda;
q2=zR+1j*pi*w0^2/lambda;
k=2*pi/lambda;
a=diam/2;

m0=abs(q2/q1); x0=m0*ksi0; y0=m0*eta0;
[X,Y]=meshgrid(xv-x0,yv-y0); Rp=sqrt(X.^2+Y.^2);
%Rb=1.22*lambda/diam*abs(z0);
%indB=(Rp<1.2*Rb);
Rb=2.23*lambda/diam*abs(z0);
indB=(Rp<Rb);


NpxDone=0;
NpxTot=sum(sum(indB));


fieldR=nan(Ly,Lx);
fieldS=nan(Ly,Lx);
uS=etd(clock,NpxDone,NpxTot,60);

for i=1:Ly
    y=yv(i);

    for j=1:Lx
        x=xv(j);
        fieldR(i,j)=B/abs(q2)*exp(1j*k*abs(zR))*exp(-1j*k*(x^2+y^2)/2/q2);
        
        if indB(i,j)    
            
            % cos(th)=z0/r
            P=z0*1j*B/lambda/abs(q1)*exp(1j*k*abs(zR-z0));
            intfun = @(ksip,etap) apod2(ksip,etap,0,0,a,0.01*a).*...
                exp(-1j*k/2/q1*((ksip+ksi0).^2+(etap+eta0).^2)).*...
                exp(1j*k*sqrt(z0^2+(x-ksi0-ksip).^2+(y-eta0-etap).^2))./(z0^2+(x-ksi0-ksip).^2+(y-eta0-etap).^2);
            
            % cos(th)=1
%             P=1j*B*/lambda/abs(q1)*exp(1j*k*abs(zR-z0));
%             intfun = @(ksip,etap) apod2(ksip,etap,0,0,a,0.001*a).*...
%                 exp(-1j*k/2/q1*((ksip+ksi0).^2+(etap+eta0).^2)).*...
%                 exp(1j*k*sqrt(z0^2+(x-ksi0-ksip).^2+(y-eta0-etap).^2))./sqrt(z0^2+(x-ksi0-ksip).^2+(y-eta0-etap).^2);

            fieldS(i,j)=P*integral2(intfun,-1.2*a,1.2*a,-1.2*a,1.2*a,'AbsTol',1e-11,'RelTol',1e-3);
            
            NpxDone=NpxDone+1;
            uS=etd(uS,NpxDone);
        end
        
    end
    
end

field=fieldR-fieldS;


if exist('outputPath','var')
    fileName=[outputPath,filesep,sprintf('d%02d_zh%02d_x%02d_y%02d',...
        round(diam*1e6),round(abs(z0)*1e3),round(ksi0*1e4),round(eta0*1e4))];
    save(fileName,'fieldS','z0','ksi0','eta0','diam',...
        'zR','xR','yR','lambda','w0','q1','q2','m0',...
        'Rb','x0','y0','xv','yv','Nx','Ny','dx','dy')
end

end