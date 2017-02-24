
% pathHolo='/home/pracownicy/jnowak/holo/ProSilica/27/holograms/';
% pathBeam='/home/pracownicy/jnowak/holo/beam/beamParams_27.mat';
% outputmat='/home/pracownicy/jnowak/holo/noise/noiseEst27.mat';
% outputplots='/home/pracownicy/jnowak/holo/noise/plotsEst';
% prefix='27_';

M=7;
windows=[20 40 80 120];
Nw=length(windows);

fileList=dir([pathHolo,'*.png']);
Nf=length(fileList);

beam=load(pathBeam);
cR=beam.cR(:,1); wR=beam.wR(:,1);
xR=beam.xR(:,1); yR=beam.yR(:,1);
x=(-1023.5:1:1023.5)'*3.1843e-6;
y=(-1023.5:1:1023.5)'*3.1843e-6;
[X,Y]=meshgrid(x,y);


sigmaNa=nan(Nf,Nw); sigmaNb=nan(Nf,Nw);
sigmaNL2=nan(Nf,Nw);
sigmaCa2=nan(Nf,Nw); sigmaCb2=nan(Nf,Nw);

us=etd(clock,0,Nf,20);
for i=1:Nf
    iRaw=double(imread([pathHolo,fileList(i).name]));
    iFit=cR(i)*exp(-2*((X-xR(i)).^2+(Y-yR(i)).^2)/wR(i)^2);
    
    bg=zeros([M size(iRaw)]);
    if i<M/2, f1=1; f2=M;
    elseif i>Nf-M/2+1, f1=Nf-M+1; f2=Nf;
    else f1=i-floor(M/2); f2=i+floor(M/2);
    end
    for j=1:M, bg(j,:,:)=imread([pathHolo,fileList(f1+j-1).name]); end
    bg=sort(bg,1);
    iMedian=squeeze(bg(ceil(M/2),:,:));
    
    for j=1:Nw
        % sigma N a
        iMean=meanMap(iRaw,windows(j));
        Xm=meanMap(X,windows(j)); Ym=meanMap(Y,windows(j));
        iFitMean=cR(i)*exp(-2*((Xm-xR(i)).^2+(Ym-yR(i)).^2)/wR(i)^2);
        iMedianMean=meanMap(iMedian,windows(j));
        
        sigmaNamap=iMean./iFitMean-1;
        sigmaNa(i,j)=mean(sigmaNamap(:));
        
        % sigma N b
        sigmaNbmap=(iMean-iMedianMean)./iFitMean;
        sigmaNb(i,j)=mean(sigmaNamap(:));
        
        % sigma NL
        sigmaNL2map=varianceMap((iRaw-iMedian)./iFit,windows(j));
        sigmaNL2(i,j)=mean(sigmaNL2map(:));
        
        % sigma C a
        sigmaCa2map=varianceMap(iRaw-iFit,windows(j))-varianceMap(iRaw-iMedian,windows(j));
        sigmaCa2(i,j)=mean(sigmaCa2map(:));
        
        % sigma C b
        sigmaCb2map=varianceMap(iMedian-iFit,windows(j));
        sigmaCb2(i,j)=mean(sigmaCb2map(:));
    end
    us=etd(us,i);
end

save(outputmat,'sigmaNa','sigmaNb','sigmaNL2','sigmaCa2','sigmaCb2',...
    'M','windows','path*')



% plots
frm='png'; res=300;

sigmaNa2=sigmaNa.^2; sigmaNb2=sigmaNb.^2;
sigmaNL=nan(size(sigmaNL2)); sigmaNL(sigmaNL2>=0)=sqrt(sigmaNL2(sigmaNL2>=0));
sigmaCa=nan(size(sigmaCa2)); sigmaCa(sigmaCa2>=0)=sqrt(sigmaCa2(sigmaCa2>=0));
sigmaCb=nan(size(sigmaCb2)); sigmaCb(sigmaCb2>=0)=sqrt(sigmaCb2(sigmaCb2>=0));

vars={'sigmaNa','sigmaNb','sigmaNL2','sigmaCa2','sigmaCb2',...
    'sigmaNa2','sigmaNb2','sigmaNL','sigmaCa','sigmaCb'};
labels={'\sigma_{Na}','\sigma_{Nb}','\sigma_{NL}^2','\sigma_{Ca}^2','\sigma_{Cb}^2',...
    '\sigma_{Na}^2','\sigma_{Nb}^2','\sigma_{NL}','\sigma_{Ca}','\sigma_{Cb}'};
factor=[1 1 1 1 1,...
    1 1 1 1 1];
legStr=num2str(windows','window %3d');
ifleg=[1 0 0 0 0,...
    1 0 0 0 0];

time=(1:Nf)'/15;
for i=1:numel(vars)
    f=figure('Color','white',...
            'PaperUnits','centimeters',...
            'PaperSize',[21 29.7],...
            'PaperPosition',[2.5 2.5 16 7]);
    ax=axes('Parent',f,'Color','none','Box','off',...
        'Position',[0.10 0.13 0.87 0.83],'FontSize',8);
    hold on
    
    plot(time,eval(vars{i})*factor(i))
    xlabel('Time [s]')
    ylabel(labels{i})
    
    if ifleg(i), legend(legStr,'Location','northeast'); end
    set(ax,'XLim',[time(1) time(end)],'XGrid','on','YGrid','on','GridAlpha',0.5)
    
    print(f,[outputplots,filesep,prefix,vars{i}],['-d',frm],['-r',num2str(res)]);
end
    

        
       