function hHst = pixelProcess(ImgD,hMP,hHst,nFrames,hPl,hCCf,PrC,grim,PcYLim,PixRng,HistRng)


    %each row: [freq, alt] <--this is to keep commonality with the x,y
    %ordering that mouse-clicking uses, since some code is shared in common
        hHst.PxPk(1:length(hHst.histAlt),1) = hHst.histFreq(1); 
        hHst.PxPk(:,2) = hHst.histAlt;
%         PxPk2 = [histAlt(1),histFreq(2);...
%                  histAlt(2),histFreq(2);...
%                  histAlt(3),histFreq(2)];

%% find indicies of ImgD corresponding to selected alt, freq
hHst.iPxPk = findPixelIndex(hPl.AltInt,hPl.freqLin,hHst.PxPk);

%% setup         
hHst.nPP = size(hHst.PxPk,1); %number of pixel processes to plot

PxImCmap = [ 1 0 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
         
  fPxPc(1) = figure('name','Pixel Processes','numbertitle','off',...
      'position',[50,50,560,420],'toolbar','none');
%   fPxPc(2) = figure('name','Pixel Processes','numbertitle','off',...
%       'position',[575+524,5,560,420]);
  for k = 1:hHst.nPP
      %get data for all frames at this specific pixel
      ImgAtPixel = squeeze( ImgD( hHst.iPxPk(k,1),  hHst.iPxPk(k,2), : ) );
      
      %make axes to plot data
      aPxPc(k,1) = subplot(hHst.nPP,1,k,'parent',fPxPc(1)); 
%       aPxPc(k,2) = subplot(nPP,1,k,'parent',fPxPc(2));
      
%% plot process at a pixel
      pPxPc(k,1) = plot(1:nFrames, ImgAtPixel,...
          'parent',aPxPc(k,1),'linestyle','-','marker','.','markersize',2);
%      pPxPc(k,2) = plot( tmp2,'parent',aPxPc(k,2),'linestyle','-','marker','.',...
%           'markersize',2);
      
%% label process parameters
      text(0.8,0.80,{['\mu= ',num2str(mean(ImgAtPixel),'%2.1e'),' '],...
                     ['\sigma= ',num2str(std(ImgAtPixel),'%2.1e'),' ']},'parent',aPxPc(k,1))
%       text(0.8,0.80,{['\mu= ',num2str(mean(tmp2),'%2.1e'),' '],...
%                      ['\sigma= ',num2str(std(tmp2),'%2.1e'),' ']},'parent',aPxPc(k,2))
                 
      set(aPxPc(k,1),'ylim',PcYLim)
%       set(aPxPc(k,2),'ylim',PcYLim)
      
%% create sliding marker for pixel process
      hHst.pPxTk(k) = line(nan,nan,'color',PxImCmap(k,:),'marker','+',...
          'markersize',5,'parent',aPxPc(k,1));
%       pPxTk2(k) = line(nan,nan,'color',PxImCmap(k,:),'marker','+',...
%           'markersize',5,'parent',aPxPc(k,2));
      
%% show pixel x,y markers on main "Data Playback" image
      line(hHst.PxPk(k,1),hHst.PxPk(k,2),'color',PxImCmap(k,:),'marker','.',...
          'parent',hMP.ax,'linestyle','none','markersize',8);
%       line(PxPk2(k,2),PxPk2(k,1),'color',PxImCmap(k,:),'marker','.',...
%           'parent',hMP.ax,'linestyle','none');
      
      title(aPxPc(k,1),['Pixel Process at: ',num2str(hHst.PxPk(k,1)),' MHz , ',...
          num2str(hHst.PxPk(k,2)),' km'])
%       title(aPxPc(k,2),['Pixel Process at: (r,c)=(',...
%           int2str(PxPk2(k,1)),',',int2str(PxPk2(k,2)),')',])
  end
  xlabel(aPxPc(end,1),'Frame #'),ylabel(aPxPc(end,1),'Pix Value')
  
  %% histogram
if PrC(9)

        histCoord = [568     600 560 420;...
             568+524 600 560 420];
         
    hHst.h = makehist(HistRng,ImgD,PrC(1),hHst.histFreq,histCoord(1,:));
    %hHst(2) = makehist(PixRng,InterpImgD,PrC(1),histFreq,histCoord(2,:));
else 
end
%==========================================================================
%% Step 1: Gaussian initialization


%initialize Gaussians
initFactors = grim.initFactors;
K = length(initFactors);
display(['Initialization GMM factors for ',hPl.fn,' : ',num2str(initFactors)])
muG(1,:)  = initFactors .* mean(mean(ImgD(:,:,1)));

sigG = nan(nFrames,K); 
sigG(1,:) = initFactors .* std(std(ImgD(:,:,1)));

rho = nan(nFrames,1); B = rho;
wgtG = nan(nFrames,K);    wgtG(1,:) = repmat( 1/K ,1,K);

Hit = false(hPl.nRow,hPl.nCol,nFrames); %initialize motion-detected array
RealHit = Hit;
iWpixel = nan(hPl.nRow,hPl.nCol,nFrames); %initial weight discard indicies array
iWpixel(:,:,1) = 1; %<--this is OK since all processes are weighted equally to start--have
%to initialize to something since we update to "t+1" (t+0 needs to have a
%starting value)
  %% control row, col
  
if PrC(3) %learn @ single pixel
    iRow = hHst.iPxPk(1,1); %uses first altitude choice
    iCol = hHst.iPxPk(1,2); %uses first frequency choice
else %learn @ all pixels
    iRow = 1:hPl.nRow;
    iCol = 1:hPl.nCol;
end
  
%% Grimson method
Actual =false;

ccFN = [hPl.fn(1:end-4),'_cc.mat'];
if ~exist(ccFN,'file')
    hWt = waitbar(0,'please wait');
for row = iRow
    for col = iCol
        for t = 1:nFrames-1
            %compute pixel parameters for all frames
        [wgtG(t+1,:), muG(t+1,:), sigG(t+1,:), rho(t,1),...
         Hit(row,col,t+1), iWpixel(row,col,t+1), RealHit(row,col,t+1)] =...
            PixelWeights(ImgD(row,col,t+1), wgtG(t,:), muG(t,:), sigG(t,:),...
            K, grim.lrnR, grim.sigmaThres,grim.BackgroundThres);
        % determine which pixels are background for each frame
         %B(t+1,1) = BackgroundModes(wgtG(t+1,:),sigG(t+1,:),grim.BackgroundThres);
        end
    if Actual, HitFrm{row,col} = find(RealHit(row,col,:)); %get time indicies of Hits
    else       HitFrm{row,col} = find(    Hit(row,col,:)); %get time indicies of Hits
    end
    
    end
    waitbar(row/hPl.nRow,hWt,'Processing Gaussians')
end
    try close(hWt),end
  save(ccFN,'RealHit','Hit','HitFrm','wgtG','muG','sigG','rho','iWpixel','B','initFactors')
  
else load(ccFN)
end


%% plot Grimson results for full frame
%figure,plot(squeeze(ImgD(60,60,:)))

%% plot Grimson results for single location
if PrC(3)

 %configure legend text
for i = 1:K, legTitle{i} = num2str(muG(1,i),'%2.1e'); end

% ======= plot learning Rho at a pixel
hHst.fRW = figure('toolbar','none','numbertitle','off','name','Weights & Learning Param.',...
    'pos',[630 650 1050 840]);
hHst.axRho = axes('parent',hHst.fRW,'pos',[0.05,0.05,0.775/2,0.815/2]);
plot(rho,'parent',hHst.axRho),
title(hHst.axRho,['Learning parameter \rho at: (r,c) = (',...
    num2str(hHst.PxPk(k,1)),',',num2str(hHst.PxPk(k,2)),')']),
xlabel(hHst.axRho,'frame #'), ylabel(hHst.axRho,'\rho')

%hHst.fWgt = figure('pos',[630 650 560 420],'name','Weights','numbertitle','off');
% ========== plot weights
hHst.axWgt = axes('parent',hHst.fRW,'pos',[0.55,0.05,0.775/2,0.815/2]);

plot(wgtG,'-','marker','none','parent',hHst.axWgt),legend(legTitle,'location','best')

title(hHst.axWgt,['Weights at: (r,c) =  (',...
    num2str(hHst.PxPk(k,1)),',',num2str(hHst.PxPk(k,2)),')'])
xlabel(hHst.axWgt,'frame #'),ylabel(hHst.axWgt,'weight \omega')

%show hits on weights
%{
line(HitFrm{iRow,iCol},... <--frame number at Hit
        wgtG(HitFrm{iRow,iCol},iWdiscard(iRow,iCol,HitFrm{iRow,iCol})),... <--discarded weight value at Hit frame #
        'linestyle','none','marker','*','color','red','markersize',8,...
        'parent',hHst.axWgt)
%}

%show weight that pixel "belongs" to
iWpixelList = squeeze(iWpixel(iRow,iCol,:));
for jj = 1:K
    wgtLine = nan(nFrames,1);  
    wgtMatch = find(jj==iWpixelList);
    for kk = 1:length(wgtMatch)
    wgtLine(wgtMatch(kk)) = wgtG(wgtMatch(kk),iWpixelList(wgtMatch(kk)));
    end
line(1:nFrames,...
     wgtLine,...
    'linestyle','-','marker','.','linewidth',2,'parent',hHst.axWgt)
end
%====== plot stem of pixel values w/hits indicated
%figure('pos',[50 5 560 420],'name','Pixel Process at single pixel','numbertitle','off')
hHst.axPPstem = axes('parent',hHst.fRW,'pos',[0.05,0.55,0.775/2,0.815/2]);
stem(squeeze(ImgD(iRow,iCol,:)),'parent',hHst.axPPstem)

%create Hit markers
line(HitFrm{iRow,iCol},... <-- Frame number
    squeeze(ImgD(iRow,iCol,HitFrm{iRow,iCol})),... <-- data value at frame #
    'linestyle','none','marker','*','color','red','markerSize',8,...
    'parent',hHst.axPPstem)
set(hHst.axPPstem,'xgrid','on')
%title(hHst.axPPstem,['Pixel Process at: (Freq[MHz],Alt[km]) = (',num2str(iRow),',',num2str(iCol),')'])
title(hHst.axPPstem,['Pixel Process at: (Freq[MHz],Alt[km]) = (',...
    num2str(hHst.PxPk(k,1)),',',num2str(hHst.PxPk(k,2)),')'])

xlabel(hHst.axPPstem,'Frame #'),ylabel(hHst.axPPstem,'Pixel Value')

if Actual
    disp(['Weights replaced: ',int2str(sum(RealHit(iRow,iCol,:))),' times.'])
else
    disp(['Weights replaced: ',int2str(sum(    Hit(iRow,iCol,:))),' times.'])
end

disp(['for (r,c) = (',int2str(iRow),',',int2str(iCol),...
    ') = (',num2str(hHst.PxPk(k,1)),',',num2str(hHst.PxPk(k,2)),') was replaced at frame #''s:'])
display(num2str(HitFrm{iRow,iCol}))
end
%% play movie
% play movie
if PrC(10)
    if Actual
        playFrames(ImgD,nFrames,hMP,hCCf,hPl,hHst,PixRng,HistRng,PrC,RealHit)
    else
        playFrames(ImgD,nFrames,hMP,hCCf,hPl,hHst,PixRng,HistRng,PrC,Hit    )
    end
end


  end