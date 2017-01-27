
% params
path='/home/pracownicy/jnowak/holo/ProSilica/28/holograms/';
output='/home/pracownicy/jnowak/holo/beam/beamParams_28.mat';

beamshape='plane';
cutoff=[0 0.3 0.5 0.7 0.9]*1e-3;
dx=3.1843e-6;


% prepare
Nfilters=length(cutoff);
fileList=dir([path,'*.png']);
Nfiles=length(fileList);

xR=zeros(Nfiles,Nfilters);
yR=zeros(Nfiles,Nfilters);
wR=zeros(Nfiles,Nfilters);
cR=zeros(Nfiles,Nfilters);
eA=zeros(Nfiles,Nfilters);
th=zeros(Nfiles,Nfilters);


% fit beam shape to pictures
uS=etd(clock,0,Nfiles,30);
for i=1:Nfiles
    im=single(imread([path,fileList(i).name]));
    
    for j=1:Nfilters
        [xR(i,j),yR(i,j),wR(i,j),cR(i,j),eA(i,j),th(i,j)]=beamFitGaussian(im,cutoff(j),dx,beamshape);
    end
    
    uS=etd(uS,i);
end


% save to file
if strcmp(beamshape,'tilted')
    save(output,'xR','yR','wR','cR','eA','th','cutoff')
else
    save(output,'xR','yR','wR','cR','cutoff')
end




% params
path='/home/pracownicy/jnowak/holo/ProSilica/29/holograms/';
output='/home/pracownicy/jnowak/holo/beam/beamParams_29.mat';

beamshape='plane';
cutoff=[0 0.3 0.5 0.7 0.9]*1e-3;
dx=3.1843e-6;


% prepare
Nfilters=length(cutoff);
fileList=dir([path,'*.png']);
Nfiles=length(fileList);

xR=zeros(Nfiles,Nfilters);
yR=zeros(Nfiles,Nfilters);
wR=zeros(Nfiles,Nfilters);
cR=zeros(Nfiles,Nfilters);
eA=zeros(Nfiles,Nfilters);
th=zeros(Nfiles,Nfilters);


% fit beam shape to pictures
uS=etd(clock,0,Nfiles,30);
for i=1:Nfiles
    im=single(imread([path,fileList(i).name]));
    
    for j=1:Nfilters
        [xR(i,j),yR(i,j),wR(i,j),cR(i,j),eA(i,j),th(i,j)]=beamFitGaussian(im,cutoff(j),dx,beamshape);
    end
    
    uS=etd(uS,i);
end


% save to file
if strcmp(beamshape,'tilted')
    save(output,'xR','yR','wR','cR','eA','th','cutoff')
else
    save(output,'xR','yR','wR','cR','cutoff')
end