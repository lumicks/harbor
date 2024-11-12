% Licensed under the BSD 2-Clause License
% 
% Copyright 2024 Archishman Ghosh
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function D = dropletbox(file,level)
%Read an image file with a certain binarization threshold level between 0
%and 1. 0.2 to 0.4 are typical values. This function is called by
%fusiondist.
image=imread(file);
format short
%Binarize image with a certain threshold level
BW=im2bw(image,level);
%Find negative image for contrast.
BW = ~BW;
BBs = [];
%Find a box encapsulating fusion droplets. The length is the fusion
%distance. If observation of this box is necessary uncomment rectangle
%function
measurements = regionprops(BW, 'BoundingBox', 'Area');
for k = 1 : length(measurements)
  thisBB = measurements(k).BoundingBox;
  %rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  %'EdgeColor','r','LineWidth',2 )
BBs = [BBs;thisBB];
end

BBwidths = BBs(:,3);
[sortedwidths,sortingindices] = sort(BBwidths,'descend');
fusionindex = sortingindices(1);
D = BBs(fusionindex,3);
%Convert pixels to microns
D = D*86.6/1000;

end
