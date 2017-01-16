% Jakub Nowak 2017 01 05


% files
path='/home/pracownicy/jnowak/holo/ProSilica/27/holograms/';
fileList=dir([path,'*.png']);
Nfiles=length(fileList);

% fileList=fileList([1; Nfiles]);
% Nfiles=length(fileList);


% filters
cutoff=[0 0.3 0.5 0.7 0.9]*1e-3;
Nfilters=length(cutoff);

% grid
dx=3.1843e-6;

uS=etd(clock,0,Nfiles,30);
xR=zeros(Nfiles,Nfilters);
yR=zeros(Nfiles,Nfilters);
wR=zeros(Nfiles,Nfilters);
cR=zeros(Nfiles,Nfilters);
eA=zeros(Nfiles,Nfilters);
th=zeros(Nfiles,Nfilters);


for i=1:Nfiles
    im=single(imread([path,fileList(i).name]));
    
    for j=1:Nfilters
        [xR(i,j),yR(i,j),wR(i,j),cR(i,j),eA(i,j),th(i,j)]=fitGaussianBeam(im,cutoff(j),dx,'tilted');
    end
    
    uS=etd(uS,i);
end

save('beamParamsEx.mat','xR','yR','wR','cR','eA','th','cutoff')