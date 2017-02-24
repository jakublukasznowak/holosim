
%% lowpass background
path='/home/pracownicy/jnowak/holo/ProSilica/28/holograms/';
outputn='/home/pracownicy/jnowak/holo/noise/noiseMapsLP_28.mat';
outputbg='/home/pracownicy/jnowak/holo/noise/bgMapsLP_28.mat';
%plots='/home/pracownicy/jnowak/holo/noise/plots28/lowpass/';

sectors=[20 50 100];
cutoffs=[0 0.5 0.9]*1e-3;

Ns=length(sectors); Nc=length(cutoffs);
fileList=dir([path,'*.png']); Nf=length(fileList);
noiseMaps=cell(Ns,Nc,Nf); bgMaps=cell(Ns,Nc,Nf);

disp('Low-pass background method ...')
uS=etd(clock,0,Nf,30);

for cnt3=1:Nf
    im=single(imread([path,fileList(cnt3).name]));
    for cnt2=1:Nc
        for cnt1=1:Ns
            [noiseMaps{cnt1,cnt2,cnt3},bgMaps{cnt1,cnt2,cnt3}]=...
                noiseMap(im,sectors(cnt1),cutoffs(cnt2));%,...
                %[plots,sprintf('sec%03d_lp%02d_holo%04d',sectors(cnt1),round(cutoffs(cnt2)*1e3*10),cnt3)]);
        end
        %close all
    end
    uS=etd(uS,cnt3);
end

save(outputn,'sectors','cutoffs','noiseMaps')
save(outputbg,'sectors','cutoffs','bgMaps')



%% median background
path='/home/pracownicy/jnowak/holo/ProSilica/28/holograms/';
outputn='/home/pracownicy/jnowak/holo/noise/noiseMapsMedian_28.mat';
outputbg='/home/pracownicy/jnowak/holo/noise/bgMapsMedian_28.mat';
%plots='/home/pracownicy/jnowak/holo/noise/plots28/median/';

sectors=[20 50 100];
medians=[7 15];

Ns=length(sectors); Nm=length(medians);
fileList=dir([path,'*.png']); Nf=length(fileList);
noiseMaps=cell(Ns,Nm,Nf); bgMaps=cell(Ns,Nm,Nf);

disp('Median background method ...')
uS=etd(clock,0,Nm*Nf,30);

for cnt2=1:Nm
    m=medians(cnt2);
    for cnt3=1:Nf
        im=imread([path,fileList(cnt3).name]);
        bg=zeros([m size(im)],'single');
        if cnt3<m/2, f1=1; f2=m;
        elseif cnt3>Nf-m/2+1, f1=Nf-m+1; f2=Nf;
        else f1=cnt3-floor(m/2); f2=cnt3+floor(m/2);
        end
        
        for i=1:m, bg(i,:,:)=imread([path,fileList(f1+i-1).name]); end
        bg=sort(bg,1);
        bg=squeeze(bg(ceil(m/2),:,:));
        
        for cnt1=1:Ns
            [noiseMaps{cnt1,cnt2,cnt3},bgMaps{cnt1,cnt2,cnt3}]=...
                noiseMap(im,sectors(cnt1),bg);%,...
                %[plots,sprintf('sec%03d_med%02d_holo%04d',sectors(cnt1),m,cnt3)]);
        end
        %close all
        uS=etd(uS,(cnt2-1)*Nf+cnt3);
    end
end

save(outputn,'sectors','medians','noiseMaps')
save(outputbg,'sectors','medians','bgMaps')
