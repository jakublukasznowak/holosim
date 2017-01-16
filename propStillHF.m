% JN 2017 01 09

Nworkers=3;
outpath=[pwd,filesep,'stillHuygensFresnel2'];

%% parameters

% droplets
d=[5 10 15 20 25]*1e-6;
z=[25 35 45 55 65]*1e-3;
x=[1 3]*1e-3;
y=[1 2]*1e-3;

% beam
zR=-69e-3; xR=0; yR=0;
zh=z+zR;

% sensor
Nx=2048; Ny=Nx;
dx=3.1843e-6; dy=dx;


%% prepare joblist

Nd=length(d);
Nz=length(z);
Np=length(x);
N=Nd*Nz*Np;
fprintf('Total number of iterations: %d\n',N);

dd=zeros(N,1);
zz=zeros(N,1);
xx=zeros(N,1);
yy=zeros(N,1);
for k=1:Nd
    for j=1:Nz
        for i=1:Np
            ind=(k-1)*Nz*Np+(j-1)*Np+i;
            dd(ind)=d(k);
            zz(ind)=zh(j);
            xx(ind)=x(i)*abs(z(j)/zR);
            yy(ind)=y(i)*abs(z(j)/zR);
        end
    end
end
    

%% parallel iteration

try
    pp=parpool('local',Nworkers);
catch
    disp('Cannot start parallel pool :(')
end

%uS=etd(clock,0,N,60);
parfor n=1:N
    HuygensFresnel(zz(n),xx(n),yy(n),dd(n),zR,xR,yR,Nx,dx,Ny,dy,outpath);
    fprintf('\n\n Iteration number %3d done !\n\n',n);
    %uS=etd(uS,n);
end