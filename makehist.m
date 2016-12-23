function hHst = makehist(PixRng,ImgD,histOfFreq,histFreq,histCoord)
%% histogram
hHst.f = figure('Name','Histogram','numbertitle','off',...
    'position',histCoord);
hHst.ax = axes('parent',hHst.f);
hHst.br = bar(PixRng,nan(1,length(PixRng)),'parent',hHst.ax);
set(hHst.ax,'xlim',[PixRng(1) PixRng(end)],...
        'ylim',[0 size(ImgD,2)/4],...
            'xlimmode','manual', 'ylimmode','manual',...
            'zlimmode','manual', 'climmode','manual',...
            'alimmode','manual','visible','on','ygrid','on','xgrid','on');
if histOfFreq
    HstTitle = ['Histogram at: ',num2str(histFreq),' MHz'];
else
    HstTitle = 'Histogram of Entire Frame';
end
hHst.HstTitle = title(hHst.ax,HstTitle);

end