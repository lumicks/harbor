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

function output = dropletrupture

% Input the filename (as in 'rupturefile.h5', the color for plotting,
% the moving average window size n (default 5000),
% the limits of the rupture period from time = tstart to tend, and
% the point t0 in the x axis after which the force begins to show change.
% Orientation is 'lt' or 'rt'. A downward dipping curve for the recorded force with trap1 will
% always have a left orientation and a upward rising curve will have a
% right orientation. radius is the bead radius.

% This will prompt for the input
prompt = {'Enter your filename with extension:',...
    'Enter the color for plotting (e.g. r/g/b):',...
    'Enter the moving average window value',...
    'Enter the start time for the rupture experiment in s:',...
    'Enter the end time for the rupture experiment in s:',...
    'Enter the time point after which rupture begins in s:'...
    'Did you pull the droplet left (enter lt) or right (enter rt)?'...
    'Enter the bead radius'};
dlgtitle = 'Input filename, plotting and fitting parameters...';
dims = [1 80];
definput = {'Filename.h5','k.','5000','0','100','10','lt','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Filename = cell2mat(answer(1));
datacolor1 = cell2mat(answer(2));
n = cellfun (@ str2num, answer(3));
tstart = cellfun (@ str2num, answer(4));
tend = cellfun (@ str2num, answer(5));
t0 = cellfun (@ str2num, answer(6));
orientation = cell2mat(answer(7));
radius = cellfun (@ str2num, answer(8));

if ~exist(Filename,'file')
    msgbox('File not found. Check the name or directory.', 'Error','error');
    return
end


F1x = h5read(Filename,'/Force HF/Force 1x');

timeindex =[0:1:length(F1x)-1];
timeindex = timeindex';
t = timeindex./78125;

%moving average
F1x = movmean(F1x,n);

figure
hold
z = F1x(t>t0);
m2 = F1x(1);

if (orientation=='rt') 
    yline(max(z)-m2,'b','LineWidth',2);
    m1 = max(z);
elseif (orientation=='lt')
    yline(min(z)-m2,'b','LineWidth',2);
    m1 = min(z);
end
yline(F1x(1)-m2,'r','LineWidth',2);

xlim([0 tend-tstart])
plot(t,F1x-m2,datacolor1)
xlabel('Time (s)')
ylabel('Force 1x (pN)')
grid on

ST = abs((m1-m2)/(2*pi*radius*1.1)); % 1.1 is a correction factor

output = ['The surface tension is ',num2str(ST),' pN/',char(181),'m'];
msgbox(output, 'Output');

%~Archishman Ghosh

end