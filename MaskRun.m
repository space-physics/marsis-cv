function MaskRun(src,eventdata,InterpImgD,hMP,hPl,IpAlt,IpFreq,nCol,nRow,fn,PrC,MS)

WF = 2.07; %from orig image plot, how much pixels are stretched horizontally for plotting purposes
HF = 4.35; % likewise, for vertical stretch

iFr = get(hPl.FrS,'value');
MaskRepNum = 4; %number of times mask pattern repeats

MaskSpacing = get(src,'value');

MaskRows = 20; %pixels

%MaskD = zeros(IpAlt(end),IpFreq(end)); %limits processing to user-selected ROI

   MaskD( 1:MaskRows, 1:MaskSpacing+1: MaskRepNum*(MaskSpacing+1)) = 1;
   
 %{
   NIMG = imnoise(InterpImgD(IpAlt,IpFreq,iFr).*1e12,'salt & pepper',0.1);
   imagesc(NIMG,'parent',hPl.ax2Dmask(3)), colormap(hPl.ax2Dmask(3),gray)
 %}
   
  TL = get(hPl.ax2Dmask(3),'pos');
   set(hPl.ax2Dmask(3),'pos',[TL(1) TL(2), WF*IpFreq(end) HF*IpAlt(end)])
     imagesc(InterpImgD(IpAlt,IpFreq,iFr),'parent',hPl.ax2Dmask(3),[0 1e-14])
     xlabel(hPl.ax2Dmask(3),'pixels'),ylabel(hPl.ax2Dmask(3),'pixels')
   colormap(hPl.ax2Dmask(3),gray)
   
   Fout = filter2(MaskD,InterpImgD(IpAlt,IpFreq,iFr),'same');
   
  TL = get(hPl.ax2Dmask(1),'pos');
  set(hPl.ax2Dmask(1),'pos',[TL(1) TL(2), WF*IpFreq(end) HF*IpAlt(end)])
   imagesc(Fout,'parent',hPl.ax2Dmask(1),[0 1e-13])
   xlabel(hPl.ax2Dmask(1),'pixels'),ylabel(hPl.ax2Dmask(1),'pixels')
  colormap(hPl.ax2Dmask(1),jet)
  
  
  TL = get(hPl.ax2Dmask(2),'pos');
  set(hPl.ax2Dmask(2),'pos',[TL(1) TL(2), WF*size(MaskD,2) WF*IpAlt(end)])
  imagesc(MaskD,'parent',hPl.ax2Dmask(2))
  xlabel(hPl.ax2Dmask(2),'pixels'),ylabel(hPl.ax2Dmask(2),'pixels')
  colormap(hPl.ax2Dmask(2),jet)

  title(hPl.ax2Dmask(1),'Filter Output')
  title(hPl.ax2Dmask(2),'Filter Mask')
  title(hPl.ax2Dmask(3),'Noisy Data')

end