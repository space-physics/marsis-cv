function [wgtNext, muNext, sigNext, rho, HitCand,iwNextPix,RealHit] = PixelWeights(NextPix,wgtCurr,muCurr,sigCurr,K,lrnR,sigmaThres,BackgroundThres)

HitCand=false; %initialize
lMatch = zeros(1,K); % Eq. 5 p.249 Grimson
wgtBackground(1:K) = false;


%% Background: determine which distributions are Background
% key point: distributions with higher weights occur more frequently-->more
% likely to be background
% also: distributions with lower std. dev. (per Grimson) are more likely to
% be background. 
%
% Thus, ordering from large to small values (descending order) of
% omega/sigma and cumsum'ing allows thresholding to determine which
% distributions must be background

% sort omega/sigma
%  from: p.249 sect. 2.2 of Grimson, also see Kaew p.2-3
[~, iJ] = sort(wgtCurr./sigCurr,2,'descend'); 

B = 1; wSum = wgtCurr(iJ(B)); %initialize
while wSum <= BackgroundThres && B+1<=K %keep cumsum'ing until Threshold is passed
wgtBackground(iJ(B)) = true; %this distrib. is declared background, since cumsum<threshold at this frame
B = B+1;
wSum = wSum + wgtCurr(iJ(B)); %wgtCurr(iJ(B)) picks the next weight to add, in descending order of wgt/sigma
end
wgtBackground(iJ(B)) = true; %<--this line captures the distrib that threw the cumsum over the threshold

%% update
        
%three possible cases for each new pixel:
% 1) NewPixel is NOT within sigmaThres std. dev. of any distrib --> Foreground 
% 2) NewPixel is closest to a BACKGROUND distrib --> Background 
% 3) NewPixel is closest to a FOREGROUND distrib --> Foreground
% Step 1: Is pixel within sigmaThres std. dev. of any distrib?

    %get distance to each distrib:
    DistribDist = abs(NextPix-muCurr); % (distance of NewPixel from means)
     
    % is this pixel within sigmaThres std. dev. of each distribution
    CloseEnough = DistribDist <= sigmaThres .* sigCurr;
        

    
        
if ~any(CloseEnough)
    %pixel was so far out that a distribution must be replaced --> pixel is foreground
    RealHit = true; 

    [~,iDiscard]= min(wgtCurr); % p.249 of Grimson "least probable distribution is replaced..."
    
    muCurr(iDiscard) = NextPix; % p.249 of Grimson "current value as its mean value"
    
    sigCurr(iDiscard) = sigmaThres*sigCurr(iDiscard); %arbitrarily chosen to give a "higher" variance as 
    %per Grimson--I made up that I'm multiplying by sigmaThres--could be some other value

           %{        
           %WRONG METHOD
           if wgtCurr(iDiscard) ~= min(wgtCurr) %the discarded weight was not the weakest to start with
              %Let's make the discard distribution the weakest weight (swap
              %weights with the weakest weight)
               [~,iJ] = sort(wgtCurr,'ascend'); 
               tmp = wgtCurr(iDiscard); %store the non-minimum weight of the Gaussian based on NextPix
               wgtCurr(iDiscard) = wgtCurr(iJ(1)); %puts the lowest weight to replacement Gaussian that was based on NextPix
               wgtCurr(iJ(1)) = tmp;  
               %Hit = true; %tracks where a pixel was so far out that a Gaussian model had to be replaced
           end
 %}          
 
    %since no match, lMatch(:) remains at 0 (from top of this function)
     
     rho = nan; %nothing to learn--we already updated muCurr, sigCurr
     
     iwNextPix = iDiscard; %the new pixel's weight index -- which distrib. pixel belongs to
 %===========================       
else %at least one distrib. was "close enough" -- is it FG or BG?
 %======         
    %see if pixel is closer to BG or FG distrib.
    [~,iNearestDistrib] = min(DistribDist(CloseEnough));
    if wgtBackground(iNearestDistrib)
        RealHit = false; %this is a background pixel 
    else
        RealHit = true; %this is a foreground pixel 
    end
 %============================= 
    %since we are in the "close enough" case, whether FG or BG, we
    %need to update the weights and models
          
    %there was a match, make matched distro 1 (p. 249 Grimson)
    lMatch(iNearestDistrib) = 1;
            
    %rho can be chosen simply as a constant, 
    %OR rho may be chosen as
    %N(X_t|mu_t,sigma_t)--which is a single number drawn from the
    %Gaussian distribution at time t.
    rho = lrnR.*(muCurr(iNearestDistrib) +...
          sigCurr(iNearestDistrib).*randn(1)); %using t-1 based on Kaew paper
             
    muNext(iNearestDistrib)  = ...
          (1-rho) .* muCurr(iNearestDistrib)  + rho .* NextPix; %this uses iNearestDistrib
       
    sigNext(iNearestDistrib) = sqrt(...
          (1-rho) .* sigCurr(iNearestDistrib)^2 + ...
          rho .* ( NextPix - muNext(iNearestDistrib) )' * (NextPix - muNext(iNearestDistrib) ) );
        
    iwNextPix = iNearestDistrib; %the new pixel's weight index -- which distrib. pixel belongs to
end
        
%ALL weights are updated (we observe that matched distribution gets an extra boost "lrnR .* lMatch")
wgtNext = (1-lrnR) .* wgtCurr  +  lrnR .* lMatch; 
        
      %  sum(wgtNext) %sum of weights should always be 1 !
        
%push forward "unmatched" parameters Mu and Sigma
muNext (~lMatch) = muCurr(~lMatch);

sigNext(~lMatch) = sigCurr(~lMatch);

end 