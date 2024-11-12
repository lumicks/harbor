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

function output = beadoscillationfit
% Input the filename (e.g. 'forcefile.h5'), an initial value for the force
% signal, an initial value for the amplitude, the frequency of
% oscillation, the bead radius r in um, and the trap stiffness Kt in pN/um.

% This will prompt for the input
prompt = {'Enter your filename with extension:',...
    'Enter the initial force value in pN:',...
    'Enter the initial amplitude value in pN',...
    'Enter the oscilation frequency in Hz:',...
    ['Enter the bead radius in ',char(181),'m:'],...
    ['Enter the trap 1x stiffness in pN/',char(181),'m:']};
dlgtitle = 'Input filename, initial values, parameters...';
dims = [1 80];
definput = {'Filename.h5','50','50','1','1','1000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Filename = cell2mat(answer(1));
initforce = cellfun (@ str2num, answer(2));
initamp = cellfun (@ str2num, answer(3));
frequency = cellfun (@ str2num, answer(4));
r = cellfun (@ str2num, answer(5));
Kt = cellfun (@ str2num, answer(6));

if ~exist(Filename,'file')
    msgbox('File not found. Check the name or directory.', 'Error','error');
    return
end


h = waitbar(0,'Extracting force and trap position data');
set(h,'Name','Calculating...');
pause(1)

Force1x = h5read(Filename,'/Force HF/Force 1x');

trapposition = h5read(Filename,'/Trap position/1X');

%The time index is converted to time in seconds with the sampling rate =
%78125 Hz.
timeindex = [0:1:length(Force1x)-1];
timeindex = timeindex';
t = timeindex./78125;

% Both the x component of force signal from trap 1 and the trap 1 position are fit to a cosine function.
options1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf,-Inf,0],...
               'Upper',[Inf,Inf,Inf,2*pi],...
               'Startpoint',[initforce initamp frequency pi],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

options2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-Inf,0],...
               'Upper',[Inf,Inf,Inf,2*pi],...
               'Startpoint',[4.2 0.125 frequency pi],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));
           
f1 = fittype('a0 + abs(a).*cos(2*pi.*f1.*t + abs(phi1))',...
    'dependent',{'Force1x'},'independent',{'t'},...
    'coefficients',{'a0','a','f1','phi1'});

f2 = fittype('b0 + abs(b).*cos(2*pi.*f2.*t + abs(phi2))',...
    'dependent',{'trapposition'},'independent',{'t'},...
    'coefficients',{'b0','b','f2','phi2'});

waitbar(.25,h,'Fitting the force oscillations');
pause(1)


[c1,gof1,output1] = fit(t,Force1x,f1,options1);
a0 = c1.a0;
a = c1.a; 
f1 = c1.f1;
phi1 = c1.phi1;
a = abs(a);
phi1 = abs(phi1);
Force1xfit = a0 + a.*cos(2*pi.*f1.*t + phi1);

waitbar(.5,h,'Fitting the trap oscillations');
pause(1)

[c2,gof2,output2] = fit(t,trapposition,f2,options2);
b0 = c2.b0;
b = c2.b; 
f2 = c2.f2;
phi2 = c2.phi2;
b = abs(b);
phi2 = abs(phi2);
trappositionfit = b0 + b.*cos(2*pi.*f2.*t + phi2);

phasediff = rad2deg(phi1-phi2);
if phasediff > 90
    msgbox(['Phase difference cannot be greater than ',char(960),'/4'], 'Error','error');
    close(h)
    return
end
    

waitbar(.75,h,'Calculating the complex moduli');
pause(1)

%Calculation of G' and G" at the given frequency. a = Ft0; b = Xt0;
%r = bead radius; Kt = trap stiffness

Ft0 = a;
Xt0 = b/0.249; % Conversion factor for voltage to microns.
Delta = phi1-phi2;
Y = Ft0/(Xt0*Kt);

G1 = (Ft0./(6*pi*r*Xt0)).*((cos(Delta)-Y)./((cos(Delta)-Y).^2+sin(Delta).^2));
G2 = (Ft0./(6*pi*r*Xt0)).*sin(Delta)./((cos(Delta)-Y).^2+sin(Delta).^2);

% Both the force signal response and the trap oscillation are plotted on
% the same graph.

waitbar(1,h,'Plotting');
pause(1)

close(h)

figure
yyaxis left
plot(t,Force1x,'k.','MarkerSize',1)
hold
plot(t,Force1xfit,'r','LineWidth',2)
xlabel('Time (s)')
ylabel('Force 1x (pN)')
grid on

trapposition = normalize(trapposition/0.249,'center','mean');
yyaxis right
plot(t,trapposition,'k.','MarkerSize',1)
hold
trappositionfit = normalize(trappositionfit/0.249,'center','mean');
plot(t,trappositionfit,'g','LineWidth',2)
xlabel('Time (s)')
ylabel('Trap position (\mum)')
grid on

ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'g';

% Output the complex moduli
output = ['At an oscillation frequency of ',...
    num2str(frequency),' Hz, the storage modulus is ',...
    num2str(G1),' Pa, and the loss modulus is ',num2str(G2),' Pa'];
msgbox(output, 'Output');

%~Archishman Ghosh

end