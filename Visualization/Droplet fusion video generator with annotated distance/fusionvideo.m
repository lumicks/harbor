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

function [] = fusionvideo(framerate,threshold)
% Create a video writer object
writerObj = VideoWriter('Fusionvideo.avi');
% Set frame rate
writerObj.FrameRate = framerate;
% Open video writer object and write frames sequentially
open(writerObj)
format short
% Read the tif image sequence exported from imageJ
Files=dir('*.tif');
for k=1:length(Files)
   Filename=Files(k).name;
   image=imread(Filename);
% Convert images to binary with a particular threshold
    BW=im2bw(image,threshold);
    BW = ~BW;
    im = imshow(image);
    BBs = [];
% Create a bounding box containing the fusing drops
    measurements = regionprops(BW, 'BoundingBox', 'Area');
    for k = 1 : length(measurements)
        thisBB = measurements(k).BoundingBox;
        BBs = [BBs;thisBB];
    end
% Sort all the widths of these bounding boxes in descending order    
    BBwidths = BBs(:,3);
    [sortedwidths,sortingindices] = sort(BBwidths,'descend');
% Choose the largest bounding box width. Everything else is noise.
    fusionindex = sortingindices(1);
    D = BBs(fusionindex,3);
% Convert to microns
    D = D*86.6/1000;
% Draw the bounding box and report the width onto the image
    x = BBs(fusionindex,1); y = BBs(fusionindex,2); w = BBs(fusionindex,3); h = BBs(fusionindex,4); 
    X = [x x+w];
    Y = [y+h y+h];
    line(X,Y,'Color',[0.0 1.0 0.0],'LineWidth',4) 
    text(120,225,[num2str(round(D,2)) '\mum'],'Color',[0.0 1.0 0.0],'FontSize',25,'Fontweight','bold')
% Export new tif series annotated with the bounding box and width value
    frame2 = [num2str(Filename) '-out.tif'];
    saveas(im,frame2)
% Write a new video with the annotated tif series
    input = imread(frame2);
    writeVideo(writerObj, input);
end
close(writerObj);
end