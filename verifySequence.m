
function [SIM,MESpaired,MESfalse] = verifySequence (maxDist,holoparams,mfile,classification,output,plotpath)

if ischar(mfile), mfile=load(mfile); end
if ischar(holoparams), load(holoparams); end
if nargin<4 || isempty(classification), classification=mfile.classification; end
if nargin<5, output=''; end
if nargin<6 || isempty(plotpath), ifplot=false; else ifplot=true; end

nubblyInd=find(classification(:,2)=='Particle_nubbly');
holonum=mfile.holonum(nubblyInd);
metrics=mfile.metrics(nubblyInd,:);
metricsNames=mfile.metricsNames;
Nmet=size(metrics,2);
traces=mfile.traces(nubblyInd,:);
tracesNames=mfile.tracesNames;
Ntr=size(traces,2);
clear mfile


inmetrNames={'zpos','xpos','ypos','minsiz','majsiz','orient'};
[~,inmetrInd]=ismember(inmetrNames,metricsNames);


SIM=struct('pos',[],'diam',[],'vel',[],'holonum',[],'indNubbly',[]);
MESpaired=struct('pos',[],'diam',[],'vel',[],'holonum',[],'indNubbly',[]);
MESfalse=struct('pos',[],'diam',[],'vel',[],'holonum',[],'indNubbly',[]);

Nh=length(holoparams);
us=etd(clock,0,Nh,20);
for n=1:Nh
   holonumInd=find(holonum==n);
   mes=cell2struct(num2cell(metrics(holonumInd,inmetrInd),1),inmetrNames,2);
   
   if ifplot
       close all
       [pairedIndHolo,pairedData,falseIndHolo,falseData]=verifyDetection(holoparams(n),mes,maxDist,...
           [plotpath,filesep,num2str(n,'verify%04d')]);
   else
       [pairedIndHolo,pairedData,falseIndHolo,falseData]=verifyDetection(holoparams(n),mes,maxDist);
   end
   
   pairedIndNubbly=zeros(size(pairedIndHolo));
   pairedIndNubbly(pairedIndHolo>0)=holonumInd(pairedIndHolo(pairedIndHolo>0));
   falseIndNubbly=holonumInd(falseIndHolo);

   SIM.pos=cat(1,SIM.pos,holoparams(n).pos);
   SIM.diam=cat(1,SIM.diam,holoparams(n).diam);
   SIM.vel=cat(1,SIM.vel,holoparams(n).vel);
   SIM.holonum=cat(1,SIM.holonum,n*ones(length(holoparams(n).diam),1));
   SIM.indNubbly=cat(1,SIM.indNubbly,pairedIndNubbly);
   
   MESpaired.pos=cat(1,MESpaired.pos,pairedData.pos);
   MESpaired.diam=cat(1,MESpaired.diam,pairedData.diam);
   MESpaired.vel=cat(1,MESpaired.vel,pairedData.vel);
   MESpaired.holonum=cat(1,MESpaired.holonum,n*ones(length(pairedData.diam),1));
   MESpaired.indNubbly=cat(1,MESpaired.indNubbly,pairedIndNubbly);

   MESfalse.pos=cat(1,MESfalse.pos,falseData.pos);
   MESfalse.diam=cat(1,MESfalse.diam,falseData.diam);
   MESfalse.vel=cat(1,MESfalse.vel,falseData.vel);
   MESfalse.holonum=cat(1,MESfalse.holonum,n*ones(length(falseData.diam),1));
   MESfalse.indNubbly=cat(1,MESfalse.indNubbly,falseIndNubbly);
   
   us=etd(us,n);
end

nzsim=(SIM.indNubbly>0);
Lsim=length(SIM.indNubbly);
SIM.indAll=zeros(Lsim,1);
SIM.indAll(nzsim)=nubblyInd(SIM.indNubbly(nzsim));
if all(SIM.pos(:,1)<0)
    SIM.pos(:,1)=-SIM.pos(:,1);
end

nzpaired=(MESpaired.indNubbly>0);
Lpaired=length(MESpaired.indNubbly);
MESpaired.indAll=zeros(Lpaired,1);
MESpaired.indAll(nzpaired)=nubblyInd(MESpaired.indNubbly(nzpaired));

MESpaired.metrics=zeros(Lpaired,Nmet);
MESpaired.metrics(nzpaired,:)=metrics(MESpaired.indNubbly(nzpaired),:);
MESpaired.traces=cell(Lpaired,Ntr);
MESpaired.traces(nzpaired,:)=traces(MESpaired.indNubbly(nzpaired),:);


MESfalse.indAll=nubblyInd(MESfalse.indNubbly);
MESfalse.metrics=metrics(MESfalse.indNubbly,:);
MESfalse.traces=traces(MESfalse.indNubbly,:);


if ~isempty(output)
    save(output,'SIM','MESpaired','MESfalse','maxDist',...
        'metricsNames','tracesNames','nubblyInd')
end

end