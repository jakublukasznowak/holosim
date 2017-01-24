function [xR,yR,wR,cR,effAlpha,theta] = beamFitGaussian (im,cutoff,dx,beamshape,ifplot)

if nargin<5
    ifplot=false;
end

if nargin<4
    beamshape='plane';
end


% fminsearch options
if ifplot
    options=optimset('Display','iter','TolX',1e-5,'PlotFcns',@optimplotfval);
else
    options=optimset('Display','final','TolX',1e-5);
end


% grid
[Ny,Nx]=size(im);
x=(-Nx/2+0.5:1:Nx/2-0.5)'*dx;
y=(-Ny/2+0.5:1:Ny/2-0.5)'*dx;
[X,Y]=meshgrid(x,y);


% filtering
if cutoff>0
   imf=double(lowPass(im,struct('cutoffLenScale',cutoff,'dx',dx,'dy',dx)));
else
   imf=double(im);
end


% fit beam profile
if strcmp(beamshape,'plane')
    
    fun1 = @(b) sum(sum(abs(...
        b(4)*exp(-2*((X-b(1)).^2+(Y-b(2)).^2)/b(3)^2)...
        -imf)));
    b=fminsearch(fun1,[-0.5e-3 1e-3 6e-3 double(mean(imf(:)))],options);
    xR=b(1); yR=b(2); wR=b(3); cR=b(4); effAlpha=nan; theta=nan;
    
elseif strcmp(beamshape,'tilted')
    
    fun2 = @(b) sum(sum(abs(...
        b(4)*exp(-2*((X-b(1))*cos(b(6))+(Y-b(2))*sin(b(6))).^2./...
        (b(3)*(1-((X-b(1))*cos(b(6))+(Y-b(2))*sin(b(6)))*b(5))).^2).*...
        exp(-2*(-(X-b(1))*sin(b(6))+(Y-b(2))*cos(b(6))).^2/b(3)^2)...
        -imf)));
    b=fminsearch(fun2,[-0.5e-3 1e-3 6e-3 double(mean(imf(:))) 1e-3 1e-3],options);
    xR=b(1); yR=b(2); wR=b(3); cR=b(4); effAlpha=b(5); theta=b(6);
    
else
    error('Invalid beam shape.')
end


% plot
if ifplot
    figure, imagesc(x,y,imf), colorbar, hold on
    plot([xR xR+0.1*wR*cos(theta)],[yR yR+0.1*wR*sin(theta)],'ro-')
    title(sprintf('Cutoff = %.1f mm',cutoff*1e3))
end
   
end
