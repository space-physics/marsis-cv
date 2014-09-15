function playFrames(ImgD,nFrames,hMP,hCCf,hPl,hHst,PixRng,HistRng,PrC,Hit)

%note, hHst.iPxPk was computed in pixelProcess
% hHst.iPxPk is freq, alt converted (interpolated) to row, col indicies for ImgD

%==== play the frames! ====
    for k = 1:nFrames
        ImgFrame = squeeze(ImgD(:,:,k));
            pause(0.02)
            
            %update main "data playback" window
            set(hMP.img,'cdata',ImgFrame)
            
            %update title with current frame #
   set(hMP.t,'String',['Data, Frame #',int2str(k),' / ',int2str(nFrames)]) %<-- this line actually consumes a lot of time!
            %set(hMP.t,'String',['Data, Frame #',int2str(k)])
%% update histograms (if present)
            if PrC(9) %user wants histograms

                
if ~PrC(1) % histogram of entire frame
    [yb,~] = hist( ImgFrame(:) , HistRng ); %compute histogram
     set(hHst.h.br,'ydata',yb) %update histogram
     set(hHst.h.HstTitle,'String',...
      ['Histogram for all frequencies. Frame #',int2str(k)])
%===========================================
else %histogram of frequency(ies)
%====== update 1st freq histogram (should make this a for loop)
    [yb,~] = hist( ImgFrame(:,hHst.iPxPk(1,2)) , HistRng ); %compute histogram
     set(hHst.h.br,'ydata',yb) %update histogram
     set(hHst.h.HstTitle,'String',...
      ['Histogram at: ',num2str(hHst.histFreq(1)),' MHz. Frame #',int2str(k)])
%====== update 2nd freq histogram (should make this a for loop)
%             [yb,~] = hist(ImgFrame(:,hHst.iPxPk(2,2)), HistRng);
%                 set(hHst2.br,'ydata',yb)
%=========
end
            end
%% update connected components plot
if PrC(5)
      set(hCCf.ccImg,'cdata',Hit(:,:,k))
end
%% update moving cursors on "Pixel Processes" plots            
          if PrC(4) 
            for i = 1:hHst.nPP
                set(hHst.pPxTk(i),'xdata',k,'ydata',ImgFrame(hHst.iPxPk(i,1),hHst.iPxPk(i,2)))  %ImgD(hHst.iPxPk(i,1),hHst.iPxPk(i,2),k))
                %set(pPxTk2(i),'xdata',k,'ydata',ImgD(PxPk2(i,1),PxPk2(i,2),k))
            end
          end
%=== done! ====          
    end
%% do tasks after video is done playing
% set histogram to be of ALL frames of data "for all time" so to speak
            if PrC(9) %user wants histograms
if ~PrC(1) % make histogram of entire frame
    [yb,~] = hist( ImgD(:) , HistRng ); %compute histogram
    set(hHst.h.br,'ydata',yb) %update histogram
    set(hHst.h.HstTitle,'String',...
      'Histogram for all frequencies, Over all Frames')
%===========================================
else %histogram of frequency(ies)
%====== update 1st freq histogram (should make this a for loop)
    [yb,~] = hist( ImgD(:,hHst.iPxPk(1,2)) , HistRng ); %compute histogram
    set(hHst.h.br,'ydata',yb) %update histogram
    set(hHst.h.HstTitle,'String',...
      ['Histogram at: ',num2str(hHst.histFreq(1)),' MHz. Over all Frames'])
%====== update 2nd freq histogram (should make this a for loop)
%             [yb,~] = hist(ImgFrame(:,histFreq(2)),PixRng);
%                 set(hHst2.br,'ydata',yb)
%=========
end
            end
end