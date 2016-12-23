function hPl = MDproc(ProgMode,fn,EdgeMode)
% Made for use with outputs of MARSIS "UserGUI.m"
% Michael Hirsch
% Jan-Apr 2012
% designed/tested on Linux, Matlab 2012a 
% REQUIRES MATLAB Image Processing Toolbox for Plasma Detection/Measurement
% operations
close all

% temp patch ---------
try
  load('alts.mat')
  hPl.AltitudeKM = AltitudeKM;
catch
  hPl.AltitudeKM = 200;
end
% END temp patch ----------

if nargin<1, ProgMode = 'histo'; end
if nargin<2 || isempty(fn)
  fn = '../marsis-utils/data/RDR601X/frm_ais_rdr_6019.mat';  %2108.dat
end
if nargin<3, EdgeMode = 'edge'; end

c = 299792458; % [m/s] vacuum speed of light
pAlt = [0 1200]; %[km] altitude range to examine for plasma freq. harmonics
pFreq= [0 6]; %[Mhz] frequency range to examine for cyclotron oscillation

%try load('defFileTifPlay.mat'),catch, fn = 'frm_ais_rdr_2107_data.mat'; end
%uigetfile({'*.mat';'*.tif';'*.mj2';'*.avi'},'Choose TIF video',fn);
if ~isempty(fn)
  save('defFileTifPlay.mat','fn')
end
%display(['Filename: ',fn])

hPl.OversampleFactor = 1; %not yet working fully (2012-03-17)
hPl.nFFTpgam = 512; %for periodogram, number of FFT
hPl.fn = fn;
%% plot control
switch ProgMode
    %========================
    % Grimson cases
    case 'full',             PrC = logical([1 0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0]); %does not show learning data
        hHst.histFreq = [2.5, 2.6]; %MHz
        hHst.histAlt  = [800, 850, 900];
        grim.lrnR = 0.05; %learning rate alpha = 1/numPreviousFramesToConsider
        grim.BackgroundThres = 0.6; %chosen arbitrarily by user--Kaewtrakulpong used 0.6
        grim.sigmaThres = 2.5; %per page 249 of original Grimson paper 1999
        grim.initFactors = [0.1  0.5 1 1.5 ];
    case 'singlePixel',      PrC = logical([1 0 1 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0]); 
        hHst.histFreq = [2.2, 2.6]; %MHz
        hHst.histAlt  = [850, 900];
        grim.lrnR = 0.1 ; %learning rate alpha = 1/numPreviousFramesToConsider
        grim.BackgroundThres = 0.6; %chosen arbitrarily by user--Kaewtrakulpong used 0.6
        grim.sigmaThres = 2.5; %per page 249 of original Grimson paper 1999
        grim.initFactors = [.01  0.1 1 5 ];
    %===========================
    case 'xcorr',            PrC = logical([0 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1 0 1]); 
        hPl.XCfundFreq = 15; %max pixel spacing (fundamental frequency) to try xcorr
    case 'histo2',          PrC = logical([1 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0]); 
    case 'histo',           PrC = logical([1 0 0 1 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0]); 
        hHst.histFreq = [2.2, 2.6]; %MHz
        hHst.histAlt  = [850, 900];
        grim.lrnR = 0.1 ; %learning rate alpha = 1/numPreviousFramesToConsider
        grim.BackgroundThres = 0.6; %chosen arbitrarily by user--Kaewtrakulpong used 0.6
        grim.sigmaThres = 2.5; %per page 249 of original Grimson paper 1999
        grim.initFactors = [.01  0.1 1 5 ];
    %case 'PlasmaFFT',        PrC = logical([0 0 0 0 0 0 1 1 0 0 1 0 1 0 0 0]);
    case 'PlasmaFFT',        PrC = logical([0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 1 0 1]);
    case 'PlasmaCepstral',   PrC = logical([0 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 1 0 1]);
    case '2Dfilter',         PrC = logical([0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 1]); 
    case 'Hough',            PrC = logical([0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 1]);
        hPl.EdgeMode = EdgeMode;
                                   %PrC = logical([0 0 0 0 0 0 1 1 0 0 1]); 
        hPl.minC = 1e-17; hPl.maxC = 1e-12;
        %hPl.EdgeMode = 'edge';%'edge';
        hPl.FilterAmplThres = 1e-26;
        hPl.vEL = [2 20]; %length of vertical erosion structuring element (pixels)
        hPl.hEL = [8 1]; %length of horizontal erosion structuring element (pixels)
        hPl.EdgeDetType = 'Canny';
        hPl.BWthres = 5e-16; %greater than this = 1, less than this = 0; used for PrC(16)
    otherwise error('unknown program mode ProgMode chosen by user')
end
%histOfFreq %false = hist of whole frame <-- PrC(1)
%loadMatFile %true = show learning data <--PrC(2) NOT USED CURRENTLY
%singlePlot %true = only a single frame--show learning information <-- PrC(3)
%pixelPlot  %true = plot single pixel processes and HISTOGRAMS<--PrC(4)
%ShowCCfig  %true = show diagnostic connected components plot <--PrC(5)
%showGrimson %true = show Grimsen paper computations <--PrC(6)
%linearData  %true = linear scale for all plots %false = log10 plot scale <--PrC(7)
%PlamsaDet % [not currently used] true = detecting plasma frequency <--PrC(8)
%showHist  % true = show histograms <--PrC(9)
%showMovie % true = show main movie <--PrC(10)
%show2DFFT = true; %show 2D FFT of data <--PrC(11)
%Cepstrum = false; %do cepstral alanysis <--PrC(12)
%FFT of sums = true; %do FFT of sums <--PrC(13)
%1D cross-correlation; <-- PrC(14)
%2D morphological processing <--PrC(15)
% Binary image thresholding <--PrC(16)
% Periodogram at row,col of image data click <---PrC(17)
% 2D filter mask convolution/correlation <--PrC(18)
% display summation over rows and columns <--PrC(19)
%% control frames
LIMnFrames = nan; %stop after this many frames
iFr = 1; %indices of frames to playback in movie
%% other parameters
efType = 'sobel';%'canny';
imType = 'separate';
cThres = [0.0, .9]; %canny
sThres = .02; %sobel 



%% setup program
% Variable "PixRng" is used in:
% makeplot: only to set upper and lower bounds of plots
% pixelProcess: N/A
% makehist: SETS HISTOGRAM BIN LOCATIONS
% playFrames: SETS HISTOGRAM BIN LOCATIONS
%============
% Variable "PcYLim" is used in:
% pixelProcess: to set axes limits for "pixel process" line plots only
% ===========
% Variable "PlYLim" is used in:
% makeControls: to set bounds for slider controls (contrast)
% makeplot: to set bounds for axes on several types of plots
% ===========
if PrC(7)
    HistRng = linspace(1e-24,1e-12,50);
    PixRng = [1e-20,1e-15];%logspace(-25,-15,100);
    PcYLim = [0 1e-14];
    PlYLim = [0 1e-11];
else
    HistRng = -20:0.25:-10;
    PixRng = [-20,-10]; %dB
    PcYLim = [-17 -13];
    PlYLim = [-17 -10];
end


set(0,'DefaulttextUnits','normalized') %make units of Text normalized
Nfreq = 160;
%% load data
Img = load(fn); 
ImgD = Img.signal_z; %.ImgS;
hPl.nRow = hPl.OversampleFactor * Nfreq; 
hPl.nCol = hPl.OversampleFactor * size(ImgD,1); 
nFramesAvail = size(ImgD,2)/Nfreq; %size(ImgD,3);
nFrames = min(nFramesAvail,LIMnFrames);
disp(['Will use frames 1 to ',int2str(nFrames),' out of ',int2str(nFramesAvail),...
    ' frames available in file:',fn])
hPl.freqMHz = Img.frequency_y(1:Nfreq)/1e3; % FIXME frequency is always constant per recent article(?)

% list of frequencies corresponding to x-pixels
hPl.freqLin = linspace(hPl.freqMHz(1),hPl.freqMHz(end),hPl.nCol)';

% list of alittudes corresponding to y-pixels
hPl.AltInt = linspace(hPl.AltitudeKM(1),hPl.AltitudeKM(end),hPl.nRow);

InterpImgD(hPl.nRow,hPl.nCol,nFrames) = nan;
for i = 1:nFrames
  % 80 x 12480    %ImgD(:,:,i)
  InterpImgD(:,:,i) = interp2(hPl.freqMHz, hPl.AltitudeKM,...
                              ImgD(:,(i-1)*Nfreq+1:i*Nfreq),...
                              hPl.freqLin,hPl.AltInt);
end
%% setup first image
MS = MonitorSizes();
%[ hHst hMP hHst2 hCCf hPl ] = makeplot(PixRng,ImgD,PrC(1),histFreq,histCoord,...
%                    PrC(2),PrC(5),PrC(9),PrC(8),PlYLim,PrC(11)); 
[hMP, hCCf, hPl] = makeplot(PixRng,InterpImgD,PrC,...
                    PlYLim,MS,hPl,fn,hPl.nRow,hPl.nCol,pAlt,pFreq);
                
hPl = makeControls(nFrames,iFr,InterpImgD,hMP,hPl,hPl.nCol,hPl.nRow,...
    PlYLim,c,fn,MS,PrC);
        
switch fn(end-2:end)
    
    case 'mat'
%=================

%           if dropFirstFrame
%               if all(ImgD(:,:,1))==0, ImgD(:,:,1) = []; end
%           end
      %  if ~PrC(7), ImgD =  log10(ImgD); end
  %===
  


% ====== Pixel Process Plot
if PrC(4)
    tic
    [hHst] = pixelProcess(InterpImgD,hMP,hHst,nFrames,hPl,hCCf,PrC,grim,PcYLim,PixRng,HistRng);
    toc
end

%==========================================================================    
    case 'mj2'
        fvMov = VideoReader(fn);
        frmRate = fvMov.FrameRate;
        display(['using recorded frame rate of: ',num2str(frmRate),' fps'])
        nFrames = fvMov.NumberOfFrames;
        vidHeight = fvMov.Height; vidWidth = fvMov.Width;
                % Preallocate movie structure.
            mov(1:nFrames) = ...
             struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
                'colormap', []);
            
            % Read one frame at a time.
        for k = 1 : nFrames
         mov(k).cdata = read(fvMov, k);
        end
        
        % Play back the movie once at the video's frame rate.
            movie(hMP.f,mov, 1, frmRate);
    case 'tif'
    Info = imfinfo(fn);
    nFrames = length(Info);
    hMP.img = imshow(nan(Info(1).Height,Info(1).Width),...
        'InitialMagnification',100);
    iptsetpref('ImshowBorder', 'tight')

    for k = 1: nFrames
            Img = imread(fn,'Index',k,'info',Info);
            pause(0.1)
            set(hMP.img,'cdata',Img)
    end
    
    otherwise, warning('Data File type was not completely defined in variable "fn"')
end

%% Measure plasma Frequency
if PrC(8)
hPl = plasmaPlay([],[],iFr,InterpImgD,hMP,hPl,hPl.IpAlt,hPl.IpFreq,hPl.nCol,hPl.nRow,fn,PrC,MS);
end
%%
if ~nargout, clear, end
end
