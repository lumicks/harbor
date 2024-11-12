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

function [Tau] = forcefit(Filename,datacolor,fitcolor)
% read hd5 file
info = h5info(Filename);
Force_Trap_2x = h5read(Filename,'/Force HF/Force 2x');
% fetch real time value from frequency
time=[0:1:length(Force_Trap_2x)-1];
time=time';
t=time./78125;
% normalize the force data for trap 2 on x axis
F2 = Force_Trap_2x;
F = F2 - min(F2(:));
Fnorm = F ./ max(F(:));
% Normalized force vs time is fit to a stretched exponential function. 
% Exponent is set at 1.5 for reasons explained in https://doi.org/10.1002/anie.202006711
options = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0,0],...
               'Upper',[Inf,Inf,Inf,Inf],...
               'Startpoint',[0.01 0.01 0.01 0.01],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

f = fittype('a + (b-a).*(1-exp(-((t-t0)./Tau).^1.5)).*heaviside(t-t0)',...
    'dependent',{'Fnorm'},'independent',{'t'},...
    'coefficients',{'a', 'b', 'Tau','t0'});
% obtain parameters
[c,gof,output] = fit(t,Fnorm,f,options); 
offset = c.a;
plateau = c.b;
Tau = c.Tau;
t0 = c.t0;
% offset the time axis and plot fitted curve and data. 
ts = t-t0;
Ffit = 0 + (1-0).*(1-exp(-((t-t0)./Tau).^1.5)).*heaviside(t-t0);
Foffset = (Fnorm-offset)./(plateau-offset);
figure
plot(ts,Foffset,datacolor,ts,Ffit,fitcolor)
grid on
xlabel('Time (s)')
ylabel('Normalized Force')

end