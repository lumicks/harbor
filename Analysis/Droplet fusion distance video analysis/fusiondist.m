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

function [Ds,Dnorm] = fusiondist(threshold)
% This function generates a fusion plot from images with edge to edge
% distance vs time, fitted to a stretched exponential function.
format short
% Directory should be empty of all other image files but the tif series
% exported from imageJ.
Files=dir('*.tif');
Ds = [];
Ts = [];
% Make an array of all edge to edge distances by calling dropletbox
% function.
for k=1:length(Files)
   Filename=Files(k).name;
   D = dropletbox(Filename,threshold);
   Ds = [Ds;D];
   t = str2num(Filename(1:5));
   Ts = [Ts;t];
end

Dscaled = Ds;
Dnorm = (Dscaled(1)-Dscaled)./(Dscaled(1)-Dscaled(end));
T = Ts-Ts(1);

% Fit edge to edge distance vs time to a stretched exponential function.
% Exponent is set at 1.5 for reasons explained in https://doi.org/10.1002/anie.202006711
options = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'Startpoint',[1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

           
f = fittype('1-exp(-((T./Tau).^1.5))',...
    'dependent',{'Dnorm'},'independent',{'T'},...
    'coefficients',{'Tau'});

[c,gof,output] = fit(T,Dnorm,f,options); 
Tau = c.Tau

Dnormfit = 1-exp(-((T./Tau).^1.5));

figure
plot(T,Ds,'k.')
xlabel('Time (s)')
ylabel('Edge to edge distance (\mum)')
grid on
figure
plot(T,Dnorm,'k.',T,Dnormfit,'r')
xlabel('Time (s)')
ylabel('Normalized Edge to edge distance')
grid on

end