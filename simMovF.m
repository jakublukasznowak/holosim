% Jakub Nowak 2017 01 13


%% parameters

outputPath='/home/pracownicy/jnowak/holo/simulations/movF2/holograms';

beamNoise=0.02;
cameraNoise=1;

dt=0.1e-3;
expT=5e-3;
lambda=532e-9;
zR=-69e-3;
beam=[zR 0 0];
Nx=2048; Ny=Nx;
dx=3.1843e-6; dy=dx;
method='Fraunhofer';

xy=[0.25 1 1 1 2 2 2 3 3 3;...
    0.25 0 0.5 1 0 1 2 0 1.5 3]'*1e-3;
z=[20 25 30 35 40 45 50 55 60 65]'*1e-3;
%z=[25 35 45 55 65]'*1e-3;

%d=(5:2.5:27.5)'*1e-6;
%d=[10 15 20]'*1e-6;
d=(5:5:25)'*1e-6;

v=(0:5:30)'*1e-3;
%v=[0 10 20 30]'*1e-3;
phi=[0 pi/4 pi/2 3/4*pi]';
theta=0;


%% prepare param mesh

% position vectors
z=zR+z;
xyz=cell2mat(mymesh(xy,z));

% velocity vectors
v_sph=cell2mat(mymesh(v,mymesh(phi,theta))); % spherical coordinates
v_car=[v_sph(:,1).*sin(v_sph(:,3)),...
    v_sph(:,1).*cos(v_sph(:,3)).*cos(v_sph(:,2)) ,...
    v_sph(:,1).*cos(v_sph(:,3)).*sin(v_sph(:,2))];
[v_car,iun]=unique(v_car,'rows','stable');
v_sph=v_sph(iun,:);
v_car=round(v_car,6);

% mesh all
params=cell2mat(mymesh(xyz,cell2mat(mymesh(d,v_car))));
pos=params(:,1:3);
diam=params(:,4);
vel=params(:,5:7);


% wrap position matrix to reduce number of holograms
Nxy=size(xy,1); Ndrop=size(pos,1);
Nholo=ceil(Ndrop/4/Nxy);
NE=Nholo*4*Nxy;
posExy=[pos(:,1:2); nan(NE-Ndrop,2)];
velExy=[vel(:,2:3); nan(NE-Ndrop,2)];

rotM=cat(3,[1 0; 0 1],[0 -1; 1 0],[-1 0; 0 -1],[0 1; -1 0]);
posExy=reshape(posExy',[2 Nxy NE/Nxy]);
velExy=reshape(velExy',[2 Nxy NE/Nxy]);
for i=1:NE/Nxy/4
    for j=1:4
        posExy(:,:,(i-1)*4+j)=rotM(:,:,j)*posExy(:,:,(i-1)*4+j);
        velExy(:,:,(i-1)*4+j)=rotM(:,:,j)*velExy(:,:,(i-1)*4+j);
    end
end
posExy=reshape(posExy,[2 NE])';
velExy=reshape(velExy,[2 NE])';

pos(:,1:2)=posExy(1:Ndrop,:);
vel(:,2:3)=velExy(1:Ndrop,:);


% scale horizontal distance and shift order to (z,x,y)
pos(:,1)=pos(:,1).*abs(pos(:,3)-zR)/abs(zR);
pos(:,2)=pos(:,2).*abs(pos(:,3)-zR)/abs(zR);
pos=circshift(pos,1,2);


% prepare parameter lists for consequtive holos
posEx=permute(reshape([pos; nan(NE-Ndrop,3)]',[3 Nxy*4 Nholo]),[2 1 3]);
diamEx=permute(reshape([diam; nan(NE-Ndrop,1)]',[1 Nxy*4 Nholo]),[2 1 3]);
velEx=permute(reshape([vel; nan(NE-Ndrop,3)]',[3 Nxy*4 Nholo]),[2 1 3]);

posC=squeeze(mat2cell(posEx,4*Nxy,3,ones(1,Nholo)));
posC{end}=posC{end}(~isnan(posC{end}(:,1)),:);
diamC=squeeze(mat2cell(diamEx,4*Nxy,1,ones(1,Nholo)));
diamC{end}=diamC{end}(~isnan(diamC{end}(:,1)));
velC=squeeze(mat2cell(velEx,4*Nxy,3,ones(1,Nholo)));
velC{end}=velC{end}(~isnan(velC{end}(:,1)),:);

holoparams=struct('pos',posC,'diam',diamC,'vel',velC,...
    'dt',dt,'expT',expT,'beam',beam,'sensor',[Nx dx Ny dy],'lambda',lambda,...
    'method',method,'beamNoiseLevel',beamNoise,'cameraNoiseLevel',cameraNoise);


%% generate and save holos
% try
%     pp=parpool('local',Nworkers);
% catch
%     disp('Cannot start parallel pool :(')
% end

fprintf('Total number of holograms : %d \n',Nholo)
uS=etd(clock,0,Nholo,30);
for i=1:Nholo
    [I,IR]=longExposure(holoparams(i).pos,holoparams(i).diam,holoparams(i).vel,...
        holoparams(i).dt,holoparams(i).expT,holoparams(i).beam,holoparams(i).sensor,...
        holoparams(i).method);
    makeHologram(I,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
        [outputPath,filesep,sprintf('holo%04d',i)],[]);
    
    makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
        [outputPath,filesep,sprintf('holo%04d_bg1',i)],[]);
    makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
        [outputPath,filesep,sprintf('holo%04d_bg2',i)],[]);
    makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
        [outputPath,filesep,sprintf('holo%04d_bg3',i)],[]);
    
    fprintf('Hologram %d out of %d done.\n\n',i,Nholo)
    uS=etd(uS,i);
end

makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
    [outputPath,filesep,sprintf('holo%04d_bg1',0)],[]);
makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
    [outputPath,filesep,sprintf('holo%04d_bg2',0)],[]);
makeHologram(IR,IR,holoparams(i).beamNoiseLevel,holoparams(i).cameraNoiseLevel,...
    [outputPath,filesep,sprintf('holo%04d_bg3',0)],[]);


save([outputPath,filesep,'holoparams.mat'],'holoparams')