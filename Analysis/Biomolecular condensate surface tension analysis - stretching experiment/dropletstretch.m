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

function [output] = dropletstretch

% Input filename as in 'stretchfile.h5'
% Diameter is the edge to edge distance of the droplet in um before stretching
% bead radius is typically 1 micron
% Kt1 and Kt2 are the x components of the stiffnesses of traps 1 and 2 in
% pN/um
% n is the number (default 5000) of points for calculating moving average
% p and q are to shift the data such that 
% they center around 0 on the y axis

% This will prompt for the input
prompt = {'Enter your filename with extension:',...
    ['Enter the droplet edge-to-edge distance in ',char(181),'m:'],...
    ['Enter the bead radius in ',char(181),'m:'],...
    ['Enter the trap 1x stiffness in pN/',char(181),'m:'],...
    ['Enter the trap 2x stiffness in pN/',char(181),'m:'],...
    'Enter the moving average window value',...
    'For plotting, shift Force 1x by:'...
    'For plotting, shift Force 2x by:'};
dlgtitle = 'Input filename, fitting and plotting parameters...';
dims = [1 80];
definput = {'Filename.h5','10','1','500','500','5000','100','100'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Filename = cell2mat(answer(1));
Diameter = cellfun (@ str2num, answer(2));
beadradius = cellfun (@ str2num, answer(3));
Kt1 = cellfun (@ str2num, answer(4));
Kt2 = cellfun (@ str2num, answer(5));
n = cellfun (@ str2num, answer(6));
p = cellfun (@ str2num, answer(7));
q = cellfun (@ str2num, answer(8));

if ~exist(Filename,'file')
    msgbox('File not found. Check the name or directory.', 'Error','error');
    return
end

h = waitbar(0,'Extracting force and trap position data');
set(h,'Name','Calculating...');
pause(1)

%parsing
F1x = h5read(Filename,'/Force HF/Force 1x');
F2x = h5read(Filename,'/Force HF/Force 2x');
trapposition = h5read(Filename,'/Trap position/1X');

timeindex = [0:1:length(F1x)-1];
timeindex = timeindex';
t = timeindex./78125;

waitbar(.25,h,'Making a moving average of force data');
pause(1)


%moving average
F1x = movmean(F1x,n);
F2x = movmean(F2x,n);
trapposition = trapposition-trapposition(1);
trapposition = trapposition/0.249; %trap position in microns

waitbar(.5,h,'Fitting the stetching data');
pause(1)


%fitting
options1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

options2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
options3 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf],...
               'Upper',[Inf,Inf],...
               'Startpoint',[1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
g1 = fittype('f1*t + b1',...
    'dependent',{'F1x'},'independent',{'t'},...
    'coefficients',{'f1','b1'});

g2 = fittype('f2*t + b2',...
    'dependent',{'F2x'},'independent',{'t'},...
    'coefficients',{'f2','b2'});

g3 = fittype('v*t + x0',...
    'dependent',{'trapposition'},'independent',{'t'},...
    'coefficients',{'v','x0'});

[c1,gof1,output1] = fit(t,F1x,g1,options1);
f1 = c1.f1;
b1 = c1.b1; 
F1xfit = f1*t + b1;

[c2,gof1,output1] = fit(t,F2x,g2,options2);
f2 = c2.f2;
b2 = c2.b2; 
F2xfit = f2*t + b2;

[c3,gof1,output1] = fit(t,trapposition,g3,options3);
v = c3.v;
x0 = c3.x0; 
trapfit = v*t + x0;

waitbar(.75,h,'Calculating the surface tension');
pause(1)


chi_sys0 = (f2-f1)/(2*-v);

chi_0 = 1/(1/chi_sys0-1/Kt1-1/Kt2);

conversion = (log(Diameter/(2*beadradius)-1)+0.68)/pi;

ST = conversion*chi_0;

waitbar(1,h,'Plotting');
pause(1)

close(h)

%plotting
figure
yyaxis left
plot(t,F1x+p,'k.','MarkerSize',1)
grid on
hold
plot(t,F1xfit+p,'r','LineWidth',2)
plot(t,F2x-q,'k.','MarkerSize',1)
plot(t,F2xfit-q,'c','LineWidth',2)
xlabel('Time (s)')
ylabel('Force 1x & 2x (pN)')

yyaxis right
plot(t,trapposition,'g.','MarkerSize',1)
plot(t,trapfit,'g','MarkerSize',2)
ylabel('Trap 1x position (\mum)')
grid on

ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'g';


output = ['The surface tension of a droplet measuring ',...
    num2str(Diameter),' ',char(181),'m is ',...
    num2str(ST),' pN/',char(181),'m'];
msgbox(output, 'Output');

%~Archishman Ghosh
end