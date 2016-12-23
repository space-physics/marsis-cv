%function [ hHst hMP hHst2 hCCf hPl] = makeplot(PixRng,ImgD,histOfFreq,histFreq,histCoord,...
%    loadMatFile,ShowCCfig,showHist,showPlasma,PlYLim)
%[ hHst hMP hHst2 hCCf hPl ] = makeplot(PixRng,ImgD,PrC(1),histFreq,histCoord,...
%                    PrC(2),PrC(5),PrC(9),PrC(8),PlYLim,PrC(11)); 
function [hMP, hCCf, hPl] = makeplot(PixRng,InterpImgD,PrC,...
                                  PlYLim,MS,hPl,fn,nRow,nCol,pAlt,pFreq)
          

%nRow = size(ImgD,1); nCol = size(ImgD,2);


freqLin = hPl.freqLin;
AltInt = hPl.AltInt;
%% get indicies for data plot region
hPl.IpAlt  = find(hPl.AltInt  >= pAlt(1)  & hPl.AltInt  <= pAlt(end));
hPl.IpFreq = find(hPl.freqLin >= pFreq(1) & hPl.freqLin <= pFreq(end));
IpAlt = hPl.IpAlt; IpFreq = hPl.IpFreq;
%% make transparency mask for image data (to give visual confirmation of data selection)
ImgMask = nan(nRow,nCol);
cMask = IpFreq(end)+1:nCol;     rMask = IpAlt(end)+1:nRow;
ImgMask(:,cMask) = 0;           ImgMask(rMask,:) = 0;
%%
freqTicks = [hPl.freqMHz(1) hPl.freqMHz(end)];

try
  MovCmap = load('MovCmap.mat') ; MovCmap = MovCmap.MovCmap;
catch
  MovCmap = 'jet';
end
try
  devCmap = load('devCmap.mat','devCmap'); devCmap = devCmap.devCmap;
catch
  devCmap = 'jet';
end
try
  BWcmap = load('BWcmap.mat','BWcmap'); BWcmap = BWcmap.BWcmap;
catch
  BWcmap = 'jet';
end
%% control box

%% Movie Image
hMP.f = figure('Position',[30 MS(1,4)-580 560 560],...
                'Toolbar','none',...
                'MenuBar','figure','Name','Data Playback','NumberTitle','off');
                
try
  set(hMP.f,'pointer','FullCrosshair')
end

hMP.ax = axes('parent',hMP.f,...
            'xlimmode','manual',...
            'ylimmode','manual',...
            'zlimmode','manual',...
            'climmode','manual',...
            'alimmode','manual','visible','off');

hMP.img = imagesc(freqLin,AltInt,nan(nRow,nCol),[PixRng(1),PixRng(end)]); %octave can't handle 'parent'

tmp = min(freqLin(end),pFreq(end));
hMP.pMask(1) = line([tmp, tmp],...
                    [AltInt(1), pAlt(end)],...
                    'color','red');
tmp = min(AltInt(end),pAlt(end));
hMP.pMask(2) = line([freqLin(1), pFreq(end)],...
                    [tmp,tmp], 'color','red');
        
   %+     hMP.mask = imshow(zeros(nRow,nCol));
   %     set(hMP.mask,'AlphaData',ImgMask)
        
%set(hMP.ax,'xtick',freqTicks)

try %matlab
    colormap(hMP.ax,MovCmap)
    hC = colorbar('peer',hMP.ax);
    set(hC,'pos',[0.91,0.15,0.0125,.75])
catch %octave
	colormap(MovCmap)
end

xlabel(hMP.ax,'Frequency (MHz)'),ylabel(hMP.ax,'Altitude (km)')
if length(fn)>16
    orbNum = fn(13:16);
else
    orbNum = [];
end

hMP.t = title(hMP.ax,['Interpolated 2D data from Orbit #',...
            orbNum,'. Upsampled by ',int2str(hPl.OversampleFactor),'x']);
%% Connected components figure 
if PrC(5)
   hCCf.f = figure('Position',[580 550 560 580],...
                'Toolbar','none',...
                'MenuBar','figure','Name','Connected Components','NumberTitle','off');
   hCCf.ax = axes('parent',hCCf.f);
   
   set(hCCf.ax, 'xlimmode','manual', 'ylimmode','manual',...
            'zlimmode','manual', 'climmode','manual',...
            'alimmode','manual','visible','off');
        
    hCCf.bgImg = imagesc(nan(size(InterpImgD,1),size(InterpImgD,2)),...
        'parent',hCCf.ax,[ PixRng(1),PixRng(end) ]);

    zeroInit = (zeros(size(InterpImgD,1),size(InterpImgD,2)));

    hCCf.ccImg = imagesc(zeroInit,[ 0,1 ]); %,'parent',hCCf.ax
    % colormap set at bottom of this function
    xlabel(hCCf.ax,'Frequency (MHz)'),ylabel(hCCf.ax,'Altitude (km)')
else
    hCCf = [];
end

%% plasma
if PrC(8)

end

%% 2D FFT
if PrC(11)    
    fwd = 0.28;%0.4; 
    xp = 0.04;
   %============= 
    hPl.fFFT = figure('Position',[50, 10, 1300, 450],...
                'Toolbar','none',...
                'MenuBar','figure','Name','2D FFT','NumberTitle','off');
            
   hPl.axFFT = axes('parent',hPl.fFFT,...
    'pos',[0*fwd+1*xp 0.15 fwd 0.77]);%[0.05 0.15 fwd 0.77]);

   hPl.axFFTslcV = axes('parent',hPl.fFFT,...
    'pos',[1*fwd+2*xp 0.15 fwd 0.77],'ydir','reverse');%[0.55 0.15 fwd 0.77]

hPl.axFFTslcH = axes('parent',hPl.fFFT,...
    'pos',[2*fwd+3*xp 0.15 fwd 0.77]);


hPl.iFFT = imagesc(nan(nRow,nCol),'parent',hPl.axFFT);
   hPl.tFFT = title(hPl.axFFT,'');
   xlabel(hPl.axFFT,'Freq. Bin #'), ylabel(hPl.axFFT,'Altitude Bin #')
   
   hPl.pFFTslV = line(nan(nRow,1),1:nRow,'parent',hPl.axFFTslcV);
   xlabel(hPl.axFFTslcV,'log10|FFT|'),ylabel(hPl.axFFTslcV,'Alt. Bin #')
   hPl.tFFTslcV = title(hPl.axFFTslcV,'');
   
   hPl.pFFTslH = line(1:nCol,nan(1,nCol),'parent',hPl.axFFTslcH);
   ylabel(hPl.axFFTslcH,'log10|FFT|'),xlabel(hPl.axFFTslcH,'Freq. Bin #')
   hPl.tFFTslcH = title(hPl.axFFTslcH,'');
  
   
   
    %============
end   

%% display summation over rows and columns
if PrC(19)
    % [73.8 47.2 434 342.3] default axis position in pixels
   hPl.fc = figure('numbertitle','off','name','Summation over Rows and Columns',...
       'position',[610,MS(1,4)-620,1000,420]);
   hPl.axr = subplot(1,2,1,'parent',hPl.fc);
   set(hPl.axr,'ylim',PlYLim,'xlim',[freqTicks(1),freqTicks(end)],...
                  'pos',[0.05 0.15 0.4 0.8] );
   hPl.pAlt = line(freqLin,nan(1,nCol),'linestyle','-','marker','none','parent',hPl.axr); %initialize line
   hPl.tPalt = title(hPl.axr,'');
   
   
   xlabel(hPl.axr,'Frequency (MHz)'), ylabel(hPl.axr,'\Sigma (altitudes)')
   %===========
   
   %============
   hPl.axc = subplot(1,2,2,'parent',hPl.fc);
   set(hPl.axc,'xlim',PlYLim,'ylim',[AltInt(1),AltInt(end)],...
       'nextplot','add','ydir','reverse','pos',[0.55,0.15,0.4,0.8]); 
   hPl.pFreq = line(nan(1,nRow),AltInt,'parent',hPl.axc);
   hPl.tPfreq = title(hPl.axc,'');
   ylabel(hPl.axc,'Altitude (km)'), xlabel(hPl.axc,'\Sigma (Frequencies)')
   %=============
end
%% sum FFT
if PrC(13)
    %========
    hPl.fFT = figure('pos',[550, MS(1,4)-1150, 1000, 420],'numbertitle','off','name','FFT of 1D Sums');
    hPl.axFTalt = subplot(1,2,1,'parent',hPl.fFT);
    set(hPl.axFTalt,'xlim',[1,nCol],'ylim',PlYLim,...
        'pos',[0.05 0.15 0.4 0.8]);
    hPl.tFTalt = title(hPl.axFTalt,'');
    hPl.pFTalt = line(1:nCol,nan(nCol,1),'parent',hPl.axFTalt);
    ylabel(hPl.axFTalt,'|FFT(\Sigma(frequencies))|'),xlabel(hPl.axFTalt,'Altitude Bin #')
    %========   

    %========
    hPl.axFTfreq = subplot(1,2,2,'parent',hPl.fFT);
    set(hPl.axFTfreq,'xlim',PlYLim,'ylim',[1,nRow],'ydir','reverse',...
        'pos',[0.55,0.15,0.4,0.8]);
    hPl.tFTfreq = title(hPl.axFTfreq,'');
    hPl.pFTfreq = line(nan(nRow,1),1:nRow,'parent',hPl.axFTfreq);
    xlabel(hPl.axFTfreq,'|FFT(\Sigma(altitudes))|'),ylabel(hPl.axFTfreq,'Frequency Bin #')
    %========
end
    
%% summed, 1D cross-correlation
if PrC(14)
    %========
    hPl.f1DC = figure('pos',[550, MS(1,4)-1150, 1000, 420],'numbertitle','off',...
        'name','Cross-correlation of 1D summations');
    hPl.ax1DCalt = subplot(1,2,1,'parent',hPl.f1DC);
    set(hPl.ax1DCalt,'xlim',[-nCol,nCol],'ylim',[0 1],'pos',[0.05 0.15 0.4 0.8]);
    hPl.t1DCalt = title(hPl.ax1DCalt,'xcorr(\Sigma(altitudes))');
    hPl.p1DCalt = line(-nCol+1:nCol-1,nan(2*nCol-1,1),'parent',hPl.ax1DCalt);
    ylabel(hPl.ax1DCalt,'xcorr(\Sigma(frequencies))'),xlabel(hPl.ax1DCalt,'Lags')
    %========   

    %========
    hPl.ax1DCfreq = subplot(1,2,2,'parent',hPl.f1DC);
    set(hPl.ax1DCfreq,'ylim',[-nRow nRow],'xlim',[0 1],'ydir','reverse','pos',[0.55,0.15,0.4,0.8]);
    hPl.t1DCfreq = title(hPl.ax1DCfreq,'xcorr(\Sigma(Frequencies))');
    hPl.p1DCfreq = line(nan(2*nRow-1,1),-nRow+1:nRow-1,'parent',hPl.ax1DCfreq);
    xlabel(hPl.ax1DCfreq,'xcorr(\Sigma(altitudes))'),ylabel(hPl.ax1DCfreq,'Lags')
    %========


    
end
%% Periodgram of row and column, by clicking on image data
if PrC(17)
   %========
    hPl.fPgam = figure('pos',[560, MS(1,4)-500, 1000, 420],'numbertitle','off',...
        'name','Periodogram from user click of Image Data');
    hPl.axPgam(1) = axes('parent',hPl.fPgam,...
        'pos',[0.05 0.15 0.4 0.8]);
    hPl.tPgam(1) = title(hPl.axPgam(1),'Periodogram(\Sigmaaltitudes)');
    
    PGL = hPl.nFFTpgam/2+1;
    
    hPl.pPgam(1) = line(1:PGL,nan(PGL,1),'parent',hPl.axPgam(1));
    xlabel(hPl.axPgam(1),'PSD: Freq. Bin #'),ylabel(hPl.axPgam(1),'PSD Estimated Magnitude')
    %========   
    
    %========
    hPl.axPgam(2) = axes('parent',hPl.fPgam,...
        'ydir','reverse','pos',[0.55,0.15,0.4,0.8]);
    hPl.tPgam(2) = title(hPl.axPgam(2),'Periodogram(\SigmaFrequencies)');
    hPl.pPgam(2) = line(nan(PGL,1),1:PGL,'parent',hPl.axPgam(2));
    xlabel(hPl.axPgam(2),'PSD Estimated Magnitude'),ylabel(hPl.axPgam(2),'PSD: Altitude Bin #')
    %========
end
%% 2D filter mask
if PrC(18)
    fwd = 0.28*1300; xp = 0.04*1300;
   %========
    hPl.f2Dmask = figure('pos',[560, MS(1,4)-500, 1300, 420],'numbertitle','off',...
        'name','2D Filter Mask technique output');
    
    hPl.ax2Dmask(1) = axes('parent',hPl.f2Dmask,'units','pixels',...
        'pos',[0*fwd+1*xp 50 fwd 350]);
    hPl.t2Dmask(1) = title(hPl.ax2Dmask(1),'2D Filter Mask output');
          
    
    %========   
    
    %========
    hPl.ax2Dmask(2) = axes('parent',hPl.f2Dmask,'units','pixels',...
        'pos',[1*fwd+2*xp 50 fwd 350]);
    hPl.t2Dmask(2) = title(hPl.ax2Dmask(2),'Mask used');
    %========
    
    %========
    hPl.ax2Dmask(3) = axes('parent',hPl.f2Dmask,'units','pixels',...
        'pos',[2*fwd+3*xp 50 fwd 350]);
    hPl.t2Dmask(3) = title(hPl.ax2Dmask(2),'Noisy Data Input');
    %========
    
    
end
%% Line Detection
if PrC(15)
    
    fwd = 0.28; xp = 0.04;
    
    hPl.fEdge = figure('pos',[70, MS(1,4)-1150, 1300, 450],'numbertitle','off',...
        'name','2D Filtering','Toolbar','none','MenuBar','figure');
    
    hPl.axEdge(3) = axes('parent',hPl.fEdge,...
    'pos',[xp, 0.15 fwd 0.77]);
    
    hPl.iEdge(3) = imagesc(freqLin,AltInt,nan(nRow,nCol),...
         'parent',hPl.axEdge(3),[-hPl.maxC hPl.maxC]);
    
    set(hPl.axEdge(3),...
            'xlimmode','manual',...
            'ylimmode','manual',...
            'zlimmode','manual',...
            'climmode','manual',...
            'alimmode','manual');
        ylabel(hPl.axEdge(3),'Altitude (km)'),xlabel(hPl.axEdge(3),'Frequency (MHz)')

    hPl.tEdge(3) = title(hPl.axEdge(3),['Eroded Data, Upsampled by Factor of ',int2str(hPl.OversampleFactor)]);
    
%=======
    hPl.axEdge(1) = axes('parent',hPl.fEdge,...
    'pos',[fwd+2*xp 0.15 fwd 0.77]);     
    hPl.iEdge(1) = imagesc(freqLin,AltInt,nan(nRow,nCol),...
         'parent',hPl.axEdge(1),[-hPl.maxC hPl.maxC]);
         set(hPl.axEdge(1),...
            'xlimmode','manual',...
            'ylimmode','manual',...
            'zlimmode','manual',...
            'climmode','manual',...
            'alimmode','manual');

    ylabel(hPl.axEdge(1),'Altitude (km)'),xlabel(hPl.axEdge(1),'Frequency (MHz)')


    %========   

    %========

    hPl.axEdge(2) = axes('parent',hPl.fEdge,...
     'pos',[2*fwd+3*xp,0.15,fwd,0.77]);
    
%========
switch hPl.EdgeMode
    case 'filter'

    
    hPl.tEdge(1) = title(hPl.axEdge(1),'Horizontal Line Filtered Output');

     hPl.iEdge(2) = imagesc(freqLin,AltInt,nan(nRow,nCol),...
         'parent',hPl.axEdge(2),[-hPl.maxC hPl.maxC]);
     
     set(hPl.axEdge(2),...
            'xlimmode','manual',...
            'ylimmode','manual',...
            'zlimmode','manual',...
            'climmode','manual',...
            'alimmode','manual');
        
         hPl.tEdge(2) = title(hPl.axEdge(2),'Vertical Line Filtered Output');
     ylabel(hPl.axEdge(2),'Altitude (km)'),xlabel(hPl.axEdge(2),'Frequency (MHz)')
    
    %colorbar
   hC =  colorbar('peer',hPl.axEdge(2),'location','southoutside');
   set(hC,'pos',[0.7,0.07,0.25,0.015]) %move over colorbar a little
   %set(hPl.axEdge(1),'pos',[0.05 0.15 0.4 0.8]) %put left axis back where it was
    case 'edge'
    hPl.tEdge(1) = title(hPl.axEdge(1),'Edge Detection');
    
    set(hPl.axEdge(2),'nextplot','add',...
         'xlim',[1,nCol],'ylim',[1,nRow],...
         'ydir','reverse');
     hPl.tEdge(2) = title(hPl.axEdge(2),'Hough Transformed Line Detection');
     %ylabel(hPl.axEdge(2),'Altitude (km)'),xlabel(hPl.axEdge(2),'Frequency (MHz)')
    ylabel(hPl.axEdge(2),'Alt. Pixel #'),xlabel(hPl.axEdge(2),'Freq. Pixel #')
end
if PrC(16)
    ccmap = BWcmap;
else
    ccmap = devCmap;
end

if PrC(5), colormap(hCCf.ax,ccmap), end

for i = 1:3
   colormap(hPl.axEdge(1),ccmap)
end
end
end
