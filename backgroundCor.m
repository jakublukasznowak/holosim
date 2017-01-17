
path='/home/pracownicy/jnowak/holo/simulations/holograms/movF_1';
output='/home/pracownicy/jnowak/holo/simulations/holograms/movF_1bc';

epsilon=0.001;
b1=20; b2=235; % brightness limits


load([path,filesep,'holoparams.mat'])
fileList=dir([path,filesep,'*.png']);
N=length(fileList);

uS=etd(clock,0,N,20);
for i=1:N
    imIn=double(imread([path,filesep,fileList(i).name]));

    if i>1
        par1=[holoparams(i).dt holoparams(i).expT holoparams(i).beam holoparams(i).sensor];
        par2=[holoparams(i-1).dt holoparams(i-1).expT holoparams(i-1).beam holoparams(i-1).sensor];
        if ~all(par1==par2,2)
            [~,IR]=longExposure([-1 1 1]*1e-3,10e-6,[0 0 0],...
                holoparams(i).dt,holoparams(i).expT,holoparams(i).beam,holoparams(i).sensor);
        end 
    else
        [~,IR]=longExposure([-1 1 1]*1e-3,10e-6,[0 0 0],...
             holoparams(i).dt,holoparams(i).expT,holoparams(i).beam,holoparams(i).sensor);
    end

    i1=min(IR(:)); i2=max(IR(:));
    imR=double(uint8((b2-b1)/(i2-i1)*IR+(b1*i2-b2*i1)/(i2-i1)));
    imR=imR/mean2(imR)*mean2(imIn);
    
    imBC=(imIn+epsilon)./(imR+epsilon)*mean2(imIn);
    
    i1=min(imBC(:)); i2=max(imBC(:));
    imOut=uint8((b2-b1)/(i2-i1)*imBC+(b1*i2-b2*i1)/(i2-i1));
    
    imwrite(imOut,[output,filesep,fileList(i).name]);
    uS=etd(uS,i);
end
