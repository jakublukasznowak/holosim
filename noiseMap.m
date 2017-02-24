
function [noisemap,bgmap] = noiseMap (im,sectorSize,cutoff,output)

dx=3.1843e-6;
if nargin<4 || isempty(output), ifplot=false; else ifplot=true; end
if nargin<3 || isempty(cutoff), cutoff=0; end

if ~isfloat(im), im=single(im); end


% subtract background
if isscalar(cutoff)
    if cutoff>0
        bg=lowPass(im,struct('cutoffLenScale',cutoff,'dx',dx,'dy',dx));
    else
        bg=zeros(size(im));
    end
elseif all(size(im)==size(cutoff),2)
    bg=cutoff;
else
    error('Wrong cutoff.')
end
imb=im-bg;

% define sectors
[L1,L2]=size(imb);
if isscalar(sectorSize), ss1=sectorSize; ss2=sectorSize;
else ss1=sectorSize(1); ss2=sectorSize(2); end

N1=floor(L1/ss1); N2=floor(L2/ss2);
rest1=L1-N1*ss1; rest2=L2-N2*ss2;
rest1left=floor(rest1/2); rest2left=floor(rest2/2);
%rest1right=rest1-rest1left; rest2right=rest2-rest2left;


% noise map
noisemap=nan([N1,N2],'single'); bgmap=nan([N1,N2],'single');
for cnt1=1:N1
    for cnt2=1:N2
        sector=imb(1+rest1left+(cnt1-1)*ss1:rest1left+cnt1*ss1,1+rest2left+(cnt2-1)*ss2:rest2left+cnt2*ss2);
        sectorbg=bg(1+rest1left+(cnt1-1)*ss1:rest1left+cnt1*ss1,1+rest2left+(cnt2-1)*ss2:rest2left+cnt2*ss2);
        noisemap(cnt1,cnt2)=std(sector(:));
        bgmap(cnt1,cnt2)=mean(sectorbg(:));
    end
end


if ifplot
    frm='png'; res=300;
    
    x=(-N1/2+0.5:1:N1/2-0.5)'*dx*ss1;
    y=(-N2/2+0.5:1:N2/2-0.5)'*dx*ss2;
    
    f=figure('Color','white',...
        'PaperUnits','centimeters',...
        'PaperSize',[21 29.7],...
        'PaperPosition',[2.5 2.5 12 10]);
    ax=axes('Parent',f,'Color','none',...
        'Position',[0.08 0.1 0.92 0.89],'FontSize',8);
    hold on
    imagesc(x*1e3,y*1e3,noisemap)
    colormap jet, colorbar('FontSize',8);
    %c.Label.String='Time [s]'; c.Label.FontSize=8;
    xlabel('x [mm]'), ylabel('y [mm]')
    set(ax,'XLim',[-3.25 3.25],'YLim',[-3.25 3.25],'FontSize',8,...
        'YDir','reverse')
    print(f,output,['-d',frm],['-r',num2str(res)]);
    
    if nargout>1
        f=figure('Color','white',...
            'PaperUnits','centimeters',...
            'PaperSize',[21 29.7],...
            'PaperPosition',[2.5 2.5 12 10]);
        ax=axes('Parent',f,'Color','none',...
            'Position',[0.08 0.1 0.92 0.89],'FontSize',8);
        hold on
        imagesc(x*1e3,y*1e3,bgmap)
        colormap jet, colorbar('FontSize',8);
        %c.Label.String='Time [s]'; c.Label.FontSize=8;
        xlabel('x [mm]'), ylabel('y [mm]')
        set(ax,'XLim',[-3.25 3.25],'YLim',[-3.25 3.25],'FontSize',8,...
            'YDir','reverse')
        print(f,[output,'bg'],['-d',frm],['-r',num2str(res)]);
    end
end


end