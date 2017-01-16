% Jakub Nowak 2017 01 04

function [field,fieldR,fieldS] = manyParticles (pos,diam,beam,sensor,method)

bankpath=['/home/pracownicy/jnowak/holo/simulations/databank/still',method];

if nargin<5
    method='Fraunhofer';
end

if length(sensor)<4
    sensor=[sensor sensor];
end

N=length(pos(:,1));
%uS=etd(clock,0,N,30);

if strcmp(method,'Fraunhofer')
    
    [~,fieldR,fieldS]=Fraunhofer(pos(1,1),pos(1,2),pos(1,3),diam(1,1),...
        beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
    
    for i=2:N
%         [state,fS]=getHoloFromBank(pos(i,1),pos(i,2),pos(i,3),diam(i,1),...
%             beam(1),sensor(1),sensor(2),sensor(3),sensor(4),bankpath);        
%         if state>0
%             newFieldS=fS;
%         else
            [~,~,newFieldS]=Fraunhofer(pos(i,1),pos(i,2),pos(i,3),diam(i,1),...
                beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
%         end
        
        fieldS=fieldS+newFieldS;
        %uS=etd(uS,i);
    end
    
elseif strcmp(method,'HuygensFresnel')
    
    [~,fieldR,~]=Fraunhofer(pos(1,1),pos(1,2),pos(1,3),diam(1,1),...
        beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
    fieldS=zeros(size(fieldR));
    
    for i=1:N
        [state,fS]=getHoloFromBank(pos(i,1),pos(i,2),pos(i,3),diam(i,1),...
            beam(1),sensor(1),sensor(2),sensor(3),sensor(4),bankpath);
        if state>0
            newFieldS=fS;
        else
            [~,~,newFieldS]=HuygensFresnel(pos(i,1),pos(i,2),pos(i,3),diam(i,1),...
                beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
            disp('Hologram calculated with HuygensFresnel propagation.')
        end
        
        fieldS=fieldS+newFieldS;
        %uS=etd(uS,i);
    end
    
elseif strcmp(method,'Fraunhofer1stOrder')
    
    [~,iR,~,Q,phi]=FraunhoferI(pos(1,1),pos(1,2),pos(1,3),diam(1,1),...
        beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
    iS=2*Q.*sin(phi);
    for i=2:N
        [~,~,~,Q,phi]=FraunhoferI(pos(i,1),pos(i,2),pos(i,3),diam(i,1),...
            beam(1),beam(2),beam(3),sensor(1),sensor(2),sensor(3),sensor(4));
        iS=iS+2*Q.*sin(phi);
        %uS=etd(uS,i);
    end
    fieldR=iR;
    fieldS=iR.*iS;
    
else
    error('Invalid method.')
end

field=fieldR-fieldS;

end