function[ FundExciteHz LifteredData ] = ...
    estimatePTperiod(xWind,idxFr,Fs,LiftHighQue,file,diag,hFig,use,MP)

Ns = length(xWind);

switch use
    case 'Freq', SearchQ = [30,.2]; PlotQ=[30,.2];
    case 'Alt', SearchQ = [3,.1]; PlotQ=[3,.1];
    case 'Speech',SearchQ = [400,60]; PlotQ = [2000,50];
end

StartQueSearch=round(Fs/SearchQ(1)); %start quefrency INDEX for  Search                          
StopQueSearch=round(Fs/SearchQ(2)); % stop quefrency INDEX for  Search

%pre allocation
C = nan(Ns,1); CL = C; FundExciteHz = nan;
LifteredData = C;

for i = 1:1
%X=fft(xWind);%.*hamming(length(xWind)));
%C=real(ifft(log(abs(X))));
C(:,i) = cceps(xWind); %used Matlab's instead to get unwrapped phase
if all(~isnan(C(:,i)))
%% lifter
%remove high-time quefrencies, so glottal impulse train period can be estimated
 CL(:,i) = lifter(C(:,i),Fs,LiftHighQue,'hicut',false); 

 LifteredData(:,i) = icceps(CL(:,i));
 %find fundamental
 %pick the quefrency with the highest amplitude -- since high-time
 %quefrencies have been liftered out, only low frequency impulse train is left
 [~,fundFreq]=max(CL(StartQueSearch:StopQueSearch,i)); 
 FundExciteHz(i,1) = Fs/(StartQueSearch+fundFreq-1);

else % no valid data this frame
    LifteredData(:,i) = NaN;
    FundExciteHz(i,1) = nan;
end
end


if diag
    
    
    
StartAx=round(Fs/PlotQ(1)); %start quefrency for plot                            
StopAx=round(Fs/PlotQ(2)); % stop quefrency for plot                 
que=(StartAx:StopAx)/Fs;

try clf(hFig) 
catch
    if hFig == 10
    figure(hFig);
    else
        pp = get(10,'pos'); 
        figure(hFig)
        set(hFig,'pos',[pp(1) pp(2)-480, 560 420])
    end
end
ax = axes('parent',hFig,'nextplot','add'); 
stem(ax,que,C(StartAx:StopAx),'r','displayname','un-liftered cepstrum')
stem(ax,que,CL(StartAx:StopAx),'displayname','Liftered Cepstrum')
title(ax,{['High-time Liftered Cepstrum into Period Estimator at Frame #: ',int2str(idxFr)'.'],...
    ['1/Quefrency cutoff: ',num2str(LiftHighQue,'%03.1f'),' [s]   File: ',file]})
xlabel(ax,'Quefrency [s]'),ylabel(ax,'Cepstral Amplitude (relative)'),legend(ax,'show')

 display(['For Frame #: ',int2str(idxFr),', Period Freq Estimate =',num2str(FundExciteHz,'%03.1f'),' Hz.']);

end
end


