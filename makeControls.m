function hPl = makeControls(nFrames,iFr,InterpImgD,...
                            hMP,hPl,nCol,nRow,PlYLim,c,fn,MS,PrC)



%% make Data Player controls


uicontrol(hMP.f,'Style','text','Pos',[5 20 50 26],'String',{'Select','Frame'},'fontsize',9)
for i = 1:nFrames
    tmp(i,:) = num2str(i,'%03d');
end

%frame select control
hPl.FrS = uicontrol(hMP.f,'Style','popup',...
    'Pos', [7 5 50 15],'String',tmp,'Value',1);
    % callback is defined at end of this function--else the callback will
    % 'forget' the hPl. struct parts defined later in this function.

%========= contrast adjust
minC = 1e-17; maxC = 1e-14;
hPl.ConSel = uicontrol(hMP.f,'Style','slider',...
                'min',minC,'max',maxC, 'sliderStep',[.02,0.1],...
                'units','normalized','pos',[.15 0 .2 .025],'value',1e-15);

   %have to do this in separate step since handle isn't yet assigned in first call
   try
     set(hPl.ConSel,'CallBack',{@ChngClimSlider,hMP.ax,false})
   end

   %label control
uicontrol(hMP.f,'Style','text','units','normalized','pos',[.175,.025,.15,.025],...
       'string','Brightness',...
       'fontsize',9);
%===========
%% put sliders on Summation panels
if PrC(19)
    %=========
   hPl.slAlt = uicontrol(hPl.fc,'Style','slider',...
                'min',1.5*min(min(InterpImgD(:,:,1))),'max',PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[0 .01 .015 .35],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slAlt,'CallBack',{@ChngMaxLimSlider,hPl.axr})
   %===========

   %=============
   hPl.slFreq = uicontrol(hPl.fc,'Style','slider',...
                'min',1.5*min(min(InterpImgD(:,:,1))),'max',PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[.505 .01 .2 .035],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slFreq,'CallBack',{@ChngMaxLimSlider,hPl.axc})
end
%% 2D FFT
%uicontrol(hPl.fFFT,'Style','text','Pos',[5 20 60 26],'String',{'Select','2Dfft Slice'},'fontsize',9)
%hPl.FFTslDir = uicontrol(hPl.fFFT,'Style','popup',...
%    'Pos', [7 5 50 15],'String',{'vert','horiz'},'Value',1);

% for Button Down Callback, see bottom of this fcn


%% 1D FFT of sums
   %=============
   if PrC(13)
   %=========
   hPl.slFTalt = uicontrol(hPl.fFT,'Style','slider',...
                'min',1.5*min(min(InterpImgD(:,:,1))),'max',10*PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[0 .01 .015 .35],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slFTalt,'CallBack',{@ChngMaxLimSlider,hPl.axFTalt})
   %===========

   %=========
   hPl.slFTfreq = uicontrol(hPl.fFT,'Style','slider',...
                'min',1.5*min(min(InterpImgD(:,:,1))),'max',10*PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[.505 .01 .2 .035],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slFTfreq,'CallBack',{@ChngMaxLimSlider,hPl.axFTfreq})
   %===========
   end
%% Cross correlation
   if PrC(14)
       %{
   %=========
   hPl.sl1DCalt = uicontrol(hPl.f1DC,'Style','slider',...
                'min',1.5*min(min(ImgD(:,:,1))),'max',10*PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[0 .01 .015 .35],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.sl1DCalt,'CallBack',{@ChngMaxLimSlider,hPl.ax1DCalt})
   %===========

   %=========
   hPl.sl1DCfreq = uicontrol(hPl.f1DC,'Style','slider',...
                'min',1.5*min(min(ImgD(:,:,1))),'max',10*PlYLim(2), 'sliderStep',[.005,0.1],...
                'units','normalized','pos',[.505 .01 .2 .035],'value',PlYLim(2));

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.sl1DCfreq,'CallBack',{@ChngMaxLimSlider,hPl.ax1DCfreq})
   %===========
   %}
       %{
   %=======
   for i = 1:hPl.XCfundFreq
    tmp2(i,:) = num2str(i,'%02d');
   end
   hPl.sel1DCxc = uicontrol(hPl.f1DC,'Style','popup',...
    'Pos', [20 5 50 15],'String',tmp2,'Value',1);
   %=======
       %}
   end
%% line detection contrast
if PrC(15)
if strcmp(hPl.EdgeMode,'filter')
%========= contrast adjust

hPl.slEdge(1) = uicontrol(hPl.fEdge(1),'Style','slider',...
                'min',hPl.minC,'max',hPl.maxC, 'sliderStep',[.01,0.1],...
                'units','normalized','pos',[.45 0 .2 .025],'value',3e-13);

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slEdge(1),'CallBack',{@ChngClimSlider,hPl.axEdge,true})

   %label control
uicontrol(hPl.fEdge(1),'Style','text','units','normalized','pos',[.475,.025,.15,.025],...
       'string','Brightness',...
       'fontsize',9);
%===========
end
% ========= Horiz/vert erosion display selection
hPl.UserErodeSel = uicontrol(hPl.fEdge,'Style','popup','units','normalized',...
    'Pos', [.01 .01 .07 .05],'String',{'Vertical','Horizontal'},'Value',1);
    % callback is defined at end of this function--else the callback will
    % 'forget' the hPl. struct parts defined later in this function.
       %label control
uicontrol(hPl.fEdge,'Style','text','units','normalized','pos',[.01,.06,.070,.06],...
       'string',{'Erosion Mode','Display'},...
       'fontsize',9);
%===========

if ~PrC(16)
%========= Eroded data contrast adjust

hPl.slEdge(2) = uicontrol(hPl.fEdge,'Style','slider',...
                'min',hPl.minC,'max',hPl.maxC, 'sliderStep',[.01,0.1],...
                'units','normalized','pos',[.10 0 .2 .025],'value',1e-13);

   %have to do this in separate step since handle isn't yet assigned in first call
   set(hPl.slEdge(2),'CallBack',{@ChngClimSlider,hPl.axEdge(3),false})

   %label control
uicontrol(hPl.fEdge,'Style','text','units','normalized','pos',[.125,.025,.15,.025],...
       'string','Brightness',...
       'fontsize',9);


   %===========
end
end

%% Callback must be the last lines in this function,
%so that all the hPl struct parts are included.
%% 2D FFT callback
if PrC(11)
    %set(hPl.FFTslDir,'callback',{@plasmaPlay,iFr,ImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,nCol,nRow,fn,PrC,MS});
    set(hPl.iFFT,'ButtonDownFcn',{@plasmaPlay,iFr,InterpImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,nCol,nRow,fn,PrC,MS});
end
%% 2D filter mask selection
if PrC(18)
uicontrol(hPl.f2Dmask,'Style','text','Pos',[5 20 50 26],'String',{'Select Mask','Spacing'},'fontsize',9)
MaxMaskSpacing = 10;
for i = 1:MaxMaskSpacing
    maskChoice(i,:) = num2str(i,'%03d');
end

%frame select control
hPl.MaskSel = uicontrol(hPl.f2Dmask,'Style','popup',...
    'Pos', [7 5 50 15],'String',maskChoice,'Value',1);
%==========
set(hPl.MaskSel,'callback',...
    {@MaskRun,InterpImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,nCol,nRow,fn,PrC,MS})
%==========
end
%% Erosion callback
if PrC(15)
set(hPl.UserErodeSel,...
    'callback',{@plasmaPlay,iFr,InterpImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,nCol,nRow,fn,PrC,MS})
end
%% popup periodgram for row, column of where clicked on main data image
if PrC(17)
set(hMP.img,...
    'ButtonDownFcn',{@PSDhere,MS,hMP,hPl,InterpImgD})
end
%% this must be last of callbacks!
set(hPl.FrS,...
    'callback',{@plasmaPlay,iFr,InterpImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,nCol,nRow,fn,PrC,MS});
end