function hPl = plasmaPlay(src,eventdata,...
    iFr,InterpImgD,hMP,hPl,pAlt,pFreq,nCol,nRow,fn,PrC,MS)


%AmplThres = Inf;

%<verify> here's a rudimentary sample frequency for alt, freq
% from original 2011 program, range of time delays is defined as:
%timeDelay = linspace(0.25,7.5,80); %ms <--from original program 2011
%Altitude = c./2.*timeDelay./1e6; %<-- my modification to 2011 program
freqLin = hPl.freqLin; %[MHz]
AltInt = hPl.AltInt;

FsDelay = 6.7908125; %1/((timeDelaySEC(end)-timeDelaySEC(1))/nRow); %#ok<COLND> %Hz, obtained from 1e-3*(7.5-0.25)/80
FsFreq = 30; %samples/period   (period is 1MHz) %(freqLin(2)-freqLin(1))*1e6;

cutoffAQ = FsDelay/10000; %<update> cepstral display parameter
cutoffFQ = FsDelay/10000; %<update> cepstral display parameter

%display(['plasmaPlay.m: Fs for time hard-coded: ',num2str(1/90.625e-6),' Hz'])

sFr = get(hPl.FrS,'value'); %get selected Frame #

if ~isempty(sFr), 
    iFr = sFr; singFr = true; 
%display(['Plotted frame #',int2str(iFr)])
else singFr = false;
end %this means the user picked a frame
                  
    
if ~singFr %true = don't bother with initialization steps
    ImgPalt(1:nCol,iFr) = NaN;  ImgPaltInterp  = ImgPalt; %initialization
    ImgPfreq(1:nRow,iFr) = NaN; %initialization
    ImgFFT = zeros(nRow,nCol);
end
    for i = iFr
           %InterpImgD = interp2(hPl.freqMHz,hPl.AltitudeKM,squeeze(ImgD(:,:,i)),freqLin,AltInt);
           InterpImgD = InterpImgD(:,:,i); 
           set(hMP.img,'cdata',InterpImgD)
   if ~singFr, pause(0.1), end         
       
     %{  
       ImgPalt(:,i)   = sum(ImgD(pAlt,:,i),1);
       ImgPfreq(:,i)  = sum(ImgD(:,pFreq,i),2);
       ImgPaltInterp  = interp1(freqMHz,ImgPalt(:,i), freqLin)';
     %}  
if PrC(19)    %compute and display summation over freq, alt
       ImgPaltInterp = sum(InterpImgD(pAlt,:),1)';
       ImgPfreq(:,i) = sum(InterpImgD(:,pFreq),2);
     
       %ImgPaltInterp(ImgPaltInterp>AmplThres) = AmplThres; %testing, used
       %to truncate large amplitudes
    
       set(hPl.pAlt,'ydata',ImgPaltInterp) %update plasma density sum plot
       set(hPl.tPalt,'string',['\Sigma (altitudes), Frame #',int2str(i)])
       set(hPl.pFreq,'xdata',ImgPfreq(:,i)) %update cyclotron sum plot
       set(hPl.tPfreq,'string',['\Sigma (frequencies), Frame #',int2str(i)])
end
       %plot(hPl.axr,freqLin,ImgPaltInterp(:,i))
%% 2D FFT       
if PrC(11) 
    
ImgFFT = fft2(InterpImgD);
       
       set(hPl.iFFT,'cdata',log10(abs(fftshift(ImgFFT))))
       set(hMP.t,'String',['Data, Frame #',int2str(i)])
       set(hPl.tFFT,'string',['2D FFT, Frame #',int2str(i)])

hFFTaxis = ancestor(src,'axes'); % get the parent axes of the image
    picked = get(hFFTaxis,'CurrentPoint'); 

if ~isempty(picked)    
    picked = round([picked(1,1) picked(1,2)]); 
    FFTslcV = ImgFFT(:,picked(1));
    FFTslcH = ImgFFT(picked(2),:);

        set(hPl.pFFTslV,'xdata',log10(abs(fftshift(FFTslcV))))

        set(hPl.pFFTslH,'ydata',log10(abs(fftshift(FFTslcH))))
        
     set(hPl.tFFTslcV,'string',['Column of 2D FFT at (r,c) = (',...
        int2str(picked(2)),',',int2str(picked(1)),')'])
    set(hPl.tFFTslcH,'string',['Row of 2D FFT at (r,c) = (',...
        int2str(picked(2)),',',int2str(picked(1)),')'])
end
end
%% Cepstral
if PrC(12) % do cepstral analysis of summed data
       hFig = 10; use = 'Freq';
      [FundFreqHz.Alt LifteredData.Alt] = ...
          estimatePTperiod(ImgPaltInterp,iFr,FsFreq,cutoffAQ,fn,true,hFig,use,MS);
      hFig = 11; use = 'Alt';
      [FundFreqHz.Freq LifteredData.Freq] = ...
          estimatePTperiod(ImgPfreq(:,i),iFr,FsDelay,cutoffFQ,fn,true,hFig,use,MS);
end
%% summed FFT
if PrC(13) %do FFT of summed INTERPOLATED data
           FTfreq  = fftshift(abs(fft(ImgPaltInterp)));
           set(hPl.pFTalt,'ydata',FTfreq)
           set(hPl.tFTalt,'string',['|FFT(\Sigma(altitudes))|, Frame #',int2str(i)])
           
           FTalt = fftshift(abs(fft(ImgPfreq(:,i))));
           set(hPl.pFTfreq,'xdata',FTalt)
           set(hPl.tFTfreq,'string',['|FFT(\Sigma(Frequencies))|, Frame #',int2str(i)])
end
%% summed cross-correlation (1D)
if PrC(14) %do 1D cross-correlation with summed INTERPOLATED data
                %setup trial vectors to cross-correlate with
               
                FundFreq = 1:hPl.XCfundFreq;
                  TryVecFreq = zeros(FundFreq(end),nCol);
                  TryVecAlt  = zeros(FundFreq(end),nRow);
                  for j = FundFreq 
                      TryVecFreq(j,1:FundFreq(j)+1:nCol) = 1;
                      TryVecAlt (j,1:FundFreq(j)+1:nRow) = 1;
                  end
                %XCalt = nan(hPl.XCfundFreq,2*nCol-1); XCfreq = nan(hPl.XCfundFreq,2*nRow-1);
                %========== TESTING ONLY =======
%                 testDataAlt = TryVecFreq(5,:);
%                 testDataFreq = TryVecAlt(5,:);
                %============================
                  for j = FundFreq
                      %========= TESTING ONLY =============
%                       [XCalt(j,:) Alags] = xcorr(testDataAlt,TryVecFreq (j,:),'coeff');
%                       [XCfreq(j,:) Flags] = xcorr(testDataFreq,TryVecAlt(j,:),'coeff');
%                       set(hPl.ax1DCalt,'ylimmode','auto')
%                       set(hPl.ax1DCfreq,'xlimmode','auto')
                      %====================================
                   %[XCalt(j,:)  Alags] = xcorr(ImgPaltInterp,TryVecFreq(j,:),'coeff');
                   %[XCfreq(j,:) Flags] = xcorr(ImgPfreq(:,i),TryVecAlt(j,:),'coeff');
                   XCalt(j) = TryVecFreq(j,:) * ImgPaltInterp;
                   XCfreq(j) = TryVecAlt(j,:) * ImgPfreq(:,i);
                  end
%    jS = get(hPl.sel1DCxc,'value');
           stem(hPl.ax1DCalt,FundFreq,XCalt)
           title(hPl.ax1DCalt,['xcorr(\Sigma(altitudes)), Frame #',int2str(i)])
           %set(hPl.p1DCalt,'ydata',XCalt,'xdata',FundFreq)%(jS,:))%,'xdata',Alags)
           %set(hPl.t1DCalt,'string',['xcorr(\Sigma(altitudes)), Frame #',int2str(i)])
           
           stem(hPl.ax1DCfreq,FundFreq,XCfreq),view(hPl.ax1DCfreq,90,90)% set(hPl.ax1DCfreq,'ydir','reverse')
           title(hPl.ax1DCfreq,['xcorr(\Sigma(frequencies)), Frame #',int2str(i)])
           %set(hPl.p1DCfreq,'xdata',XCfreq,'ydata',FundFreq)%(jS,:))%,'ydata',Flags)
           %set(hPl.t1DCfreq,'string',['xcorr(\Sigma(Frequencies)), Frame #',int2str(i)])
end
%% Periodogram
if PrC(17)
% nothing to do yet, see PSDhere.m
end
%% Line Detection
if PrC(15)
% ====== amplitude threshold data before finding edges    
     if PrC(16), amplThres = hPl.BWthres;
        ImgL = InterpImgD > amplThres; %binary image 
        ImgLthres = zeros(hPl.OversampleFactor*nRow,hPl.OversampleFactor*nCol); 
        tmp = ones(size(ImgLthres));
        ImgLthres(ImgL) = tmp(ImgL);
     else  amplThres = hPl.FilterAmplThres; %low-level noise rejection
        ImgL = InterpImgD> amplThres; %ignore pixels below this
        ImgLthres = zeros(hPl.OversampleFactor*nRow,hPl.OversampleFactor*nCol); 
        ImgLthres(ImgL) = InterpImgD(ImgL);
     end  
%==============
      
%============== do image erosion and display

%setup erosion structuring elements
horizSE = strel('rectangle',hPl.vEL);
vertSE = strel('rectangle',hPl.hEL);
    
%ImgErodedHoriz = imerode(InterpImgD,horizSE);    
%ImgErodedVert  = imerode(InterpImgD,vertSE);

ImgErodedHoriz = imerode(ImgLthres,horizSE);    
ImgErodedVert  = imerode(ImgLthres,vertSE);

UserErodeSel = get(hPl.UserErodeSel,'value');
if UserErodeSel==1, ImgEroded = ImgErodedVert;
else                ImgEroded = ImgErodedHoriz;
end
    set(hPl.iEdge(3),'cdata',ImgEroded)    
 

    switch hPl.EdgeMode
        case 'edge'
            Edges = edge(ImgEroded,hPl.EdgeDetType); % "Canny" edge detection algorithm
            set(hPl.iEdge(1),'cdata',Edges)
            set(hPl.tEdge(1),'string',...
                [hPl.EdgeDetType,' Line Detection, Frame #',int2str(i),'. Ampl. Threshold: ',num2str(amplThres,'%2.1e')])
        % Hough Transform
        [HT, Htheta, Hrho] = hough(Edges);%'theta',[-90:0.1:-89 -1:0.1:1 89:.1:89.9]);
        Hpeaks = houghpeaks(HT,5,'threshold',ceil(0.3*max(HT(:))));
        HL = houghlines(Edges,Htheta,Hrho,Hpeaks,'Fillgap',5','MinLength',4);
        
        % draw Houghlines output
        max_len = 0;
        
        delete(findobj(hPl.axEdge(2),'type','line')) %get rid of old Hough Lines
for k = 1:length(HL)
   xy = [HL(k).point1; HL(k).point2];
   plot(hPl.axEdge(2),xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');

end
            
            
        case 'filter'
        horizFiltKernel = fspecial('sobel'); % just a 3x3 filter
        horizEdges = imfilter(ImgLthres,horizFiltKernel,'same');
        set(hPl.iEdge(1),'cdata',horizEdges)
        set(hPl.tEdge(1),'string',['Horizontal Line Detection, Frame #',int2str(i)])
        
        vertFiltKernel = horizFiltKernel'; %that's all it takes per "fspecial" docs
        vertEdges = imfilter(ImgLthres,vertFiltKernel,'same');
        set(hPl.iEdge(2),'cdata',vertEdges)
        set(hPl.tEdge(2),'string',['Vertical Line Detection, Frame #',int2str(i)])

    end
    

end
    end
end