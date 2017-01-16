% Jakub Nowak 2017 01 05

beamNoiseLevel=0.02;
cameraNoiseLevel=1;

path0='/home/pracownicy/jnowak/holo/simulations/databank/stillHuygensFresnel';
path1='/home/pracownicy/jnowak/holo/simulations/holograms/stillHuygensFresnel';
%path2='/home/pracownicy/jnowak/holo/simulations/holograms/stillHuygensFresnelBG';


fileList=dir([path0,filesep,'*.mat']);
Nf=length(fileList);
uS=etd(clock,0,Nf,30);

for i=1:Nf
    name=fileList(i).name;
    load(fullfile(path0,name));
    [Ny,Nx]=size(fieldS);
    [~,fieldR,~]=Fraunhofer(z0,ksi0,eta0,diam,zR,xR,yR,Nx,dx,Ny,dy);
    
    fieldS(isnan(fieldS))=0;
    field=fieldR-fieldS;
    I=abs(field).^2;
    IR=abs(fieldR).^2;
    
    makeHologram(I,IR,beamNoiseLevel,cameraNoiseLevel,...
        [path1,filesep,name(1:end-4)],[]);%,[path2,filesep,name(1:end-4)]);
    
    uS=etd(uS,i);
end