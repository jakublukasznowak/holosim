% Jakub Nowak 2017 01 04

function [im,imB] = makeHologram (i,iR,beamNoiseLevel,cameraNoiseLevel,output,outputB)

b1=20; b2=235; % brightness limits

if nargin<6, outputB=''; end
if nargin<5, output=''; end
if nargin<4, cameraNoiseLevel=1; end
if nargin<3, beamNoiseLevel=0.03; end


beamNoise=beamNoiseLevel*randn(size(i));
cameraNoise=cameraNoiseLevel*randn(size(i));

% raw hologram
in=i+beamNoise.*iR;
i1=min(in(:)); i2=max(in(:));
im=uint8((b2-b1)/(i2-i1)*in+(b1*i2-b2*i1)/(i2-i1)+cameraNoise);

if ~isempty(output)    
    imwrite(im,[output,'.png'])
end


% background corrected hologram
inB=i./iR+beamNoise+cameraNoise./iR;

i1=min(inB(:)); i2=max(inB(:));
imB=uint8((b2-b1)/(i2-i1)*inB+(b1*i2-b2*i1)/(i2-i1));

if ~isempty(outputB)
    imwrite(imB,[outputB,'.png']);
end

end

