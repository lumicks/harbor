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

function [] = burgerfit

prompt = {'Enter the excel file name containing the viscoelastic moduli and their respective standard errors vs frequency:'};
dlgtitle = 'Input filename with extension';
dims = [1 80];
definput = {'Filename.xlsx'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
datafile = cell2mat(answer(1));

if ~exist(datafile,'file')
    msgbox('File not found. Check the name or directory.', 'Error','error');
    return
end



d = xlsread(datafile);
w = d(:,1);
G1 = d(:,2);
G2 = d(:,3);
G1err = d(:,4);
G2err = d(:,5);
%close all
h = waitbar(0,'Parsing data data');
set(h,'Name','Calculating...');
pause(1)

options1 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0,0],...
               'Upper',[Inf,Inf,Inf,Inf],...
               'Startpoint',[1 1 1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));


options2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0,0],...
               'Upper',[Inf,Inf,Inf,Inf],...
               'Startpoint',[1 1 1 1],'algorithm','Trust-Region',...
               'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
               'TolFun',(1.0E-6),'TolX',(1.0E-6));

waitbar(0.5,h,'Fitting the moduli data');
pause(1)
           
f1 = fittype('w.^2*tau0*eta0./(1+(w*tau0).^2)+w.^2*tau1*eta1./(1+(w*tau1).^2)',...
    'dependent',{'G1'},'independent',{'w'},...
    'coefficients',{'tau0','eta0','tau1','eta1'});

f2 = fittype('w*eta0./(1+(w*tau0).^2)+w*eta1./(1+(w*tau1).^2)',...
    'dependent',{'G2'},'independent',{'w'},...
    'coefficients',{'tau0','eta0','tau1','eta1'});

wfit = [w(1):0.001:w(end)]';

[c1,gof1,output1] = fit(w,G1,f1,options1);
tau0 = c1.tau0;
eta0 = c1.eta0; 
tau1 = c1.tau1;
eta1 = c1.eta1;

G1fit = wfit.^2*tau0*eta0./(1+(wfit*tau0).^2)+wfit.^2*tau1*eta1./(1+(wfit*tau1).^2);

[c2,gof2,output2] = fit(w,G2,f2,options2);
tau0 = c2.tau0;
eta0 = c2.eta0; 
tau1 = c2.tau1;
eta1 = c2.eta1;

waitbar(1,h,'Plotting');
pause(1)

G2fit = wfit*eta0./(1+(wfit*tau0).^2)+wfit*eta1./(1+(wfit*tau1).^2);

close(h)

figure
errorbar(w,G1,G1err,'bo','MarkerFaceColor','b')
hold
plot(wfit,G1fit,'b')
errorbar(w,G2,G2err,'ro','MarkerFaceColor','r')
plot(wfit,G2fit,'r')

xlabel('\omega (Hz)')
ylabel('G*(\omega) (Pa)')

axis([2.5 1250 0.01 100])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on

legend('Storage Modulus','','Loss Modulus','','Location','northwest')

output = [char(964),'0 = ',num2str(tau0),' s, ',char(951),'0 = ',num2str(eta0),' Pa.s, ',char(964),'1 = ',num2str(tau1),' s, ',char(951),'1 = ',num2str(eta1),' Pa.s '];
msgbox(output, 'Output');

end