% Jakub Nowak 2017 01 04

function [i,iR] = longExposure (pos,diam,vel,dt,expT,beam,sensor,method)

refExpT=0.1e-3; % [s]

if nargin<8
    method='Fraunhofer';
end

if length(sensor)<4
    sensor=[sensor sensor];
end

%t=(dt/2:dt:expT-dt/2)'-expT/2;
t=(0:dt:expT)'-expT/2;
Nt=length(t);
Np=length(pos(:,1));


if strcmp(method,'Fraunhofer')
    
    post=pos+vel*t(1);
    [f,fR,~]=manyParticles(post,diam,beam,sensor,'Fraunhofer');
    
    i=dt/refExpT*abs(f).^2;
    iR=dt/refExpT*abs(fR).^2;
    
    uS=etd(clock,1,Nt,30);
    
    for j=2:Nt
        post=pos+vel*t(j);
        [f,fR,~]=manyParticles(post,diam,beam,sensor,'Fraunhofer');
        
        i=i+dt/refExpT*abs(f).^2;
        iR=iR+dt/refExpT*abs(fR).^2;
        
        uS=etd(uS,j);
    end
    
elseif strcmp(method,'HuygensFresnel')
    
    post=pos+vel*t(1);
    [f,fR,~]=manyParticles(post,diam,beam,sensor,'HuygensFresnel');
    
    i=dt/refExpT*abs(f).^2;
    iR=dt/refExpT*abs(fR).^2;
    
    uS=etd(clock,1,Nt,30);
    
    for j=2:Nt
        post=pos+vel*t(j);
        [f,fR,~]=manyParticles(post,diam,beam,sensor,'HuygensFresnel');
        
        i=i+dt/refExpT*abs(f).^2;
        iR=iR+dt/refExpT*abs(fR).^2;
        
        uS=etd(uS,j);
    end
    
elseif strcmp(method,'Fraunhofer1stOrder')
    
    posList=repmat(pos,Nt,1)+repmat(vel,Nt,1).*repmat(t,Np,3);
    diamList=repmat(diam,Nt,1);
    
    [~,iR,~,Q,phi]=FraunhoferI(posList(1,1),posList(1,2),posList(1,3),diamList(1,1),...
        beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
    iS=2*Q.*sin(phi);
    
    uS=etd(clock,1,Nt*Np,30);
    
    for j=2:Nt*Np
        [~,~,~,Q,phi]=FraunhoferI(posList(j,1),posList(j,2),posList(j,3),diamList(j,1),...
            beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
        iS=iS+2*Q.*sin(phi);
        
        uS=etd(uS,j);
    end
    
    iS=iR.*iS*dt/refExpT;
    iR=Nt*dt/refExpT*iR;
    i=iR-iS; 

else
    error('Invalid method.')
end


end
            
