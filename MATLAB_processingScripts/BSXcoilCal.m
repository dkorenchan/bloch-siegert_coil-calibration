 function [w1_Hz,Ampls,Rsq] = BSXcoilCal 

% Initially written 2/10/2025 by DK
% Updated 2/13/2025 by DK

% This function takes data from a Bloch-Siegert X-nuclear coil calibration 
% (spin-echo experient with an off-resonance pulse applied during the 
% second half of the echo), and it calculates the nutation frequency as a 
% function of the pulse amplitude using the Bloch-Siegert shift dependence 
% of the signal phase
%
% INPUTS:   NONE 
%
% OUTPUTS:  
%   w1_Hz           --  Vector of the measured B1 value at each saturation
%                       offset, in units of Hz
%   Ampls           --  Vector containing the linear amplitude factors used
%                       for the Bloch-Siegert pulses in each measurement
%   Rsq             --  R-squared fitting statistic for each value of w1_Hz

% Load the data and the offsets in Hz
[rawdata,hdr,~,fname]=Read_Tecmag_DK;
disp('Reading in amplitudes and timings from table in header...')
tableVals=Read_Tecmag_Header_Table(fname,{'de0:2','rf0:3'});
HoffsetHz=abs(hdr.ob_freq(1)*1e6);
XoffsetHz=abs(hdr.ob_freq(3)*1e6);
delaysS=tableVals{1}';
Ampls=tableVals{2}';
% Check that rawdata has dimensions (spectral,ampl,TE)
if size(rawdata,2)~=length(Ampls)
    rawdata=permute(rawdata,[1,3,2]);
end

% Use the raw FIDs to obtain the phase at each offset and delay time.
% Average the phase of the FID points across all FIDs that are greater than
% the noise
pts=5:105;
% SNRthresh=5; %SNR threshold
% FIDendIdx=find(rawdata(:,1,1)~=0,1,'last');
% FIDnoise=std(reshape(real(rawdata(FIDendIdx-9:FIDendIdx,:,:)),[],1)); 
%     %last 10 points of each FID, ignoring zeros
% pts=find(prod(prod(abs(rawdata)>SNRthresh*FIDnoise,2),3));
% pts(pts<5)=[]; %remove the first 4 points
% phaseRad=squeeze(mean(unwrap(angle(rawdata(pts,:,:)),[],1),1));
phaseRad=squeeze(unwrap(angle(mean(rawdata(pts,:,:),1)),[],1));
    %this way, we average the points together to reduce the noise!
    %THEN we take the angle!

% Then go through all amplitude valuess, phase-unwrap so that we get a good  
% linear change in phase, and fit to obtain the Bloch-Siegert shift 
% (in rad/s)
disp('Fitting to determine the Bloch-Siegert shift at each B1 amplitude setting...')
w_BS=zeros(size(phaseRad,1),1);
CI_BS=zeros(size(phaseRad,1),2);
Rsq=zeros(size(phaseRad,1),1);
for ii=1:size(phaseRad,1)
    fitpts=unwrap(squeeze(phaseRad(ii,:))');
%     % Keep only the points where all FID pts contributing to the phase are
%     % above the SNR threshold
%     SNRthresh=3; %SNR threshold
%     FIDsignal=squeeze(abs(rawdata(pts,ii,:)));
%     FIDnoise=std(squeeze(rawdata(end-42:end-32,ii,:)),0,1); %last 10 points of each FID
%     SNRtPass=logical(prod(FIDsignal>SNRthresh*repmat(FIDnoise,length(pts),1),1));
%     fitpts=fitpts(SNRtPass);
%     fitdels=delaysS(SNRtPass);
    fitdels=delaysS;

    % See if further unwrapping is required, and perform
    dS=diff(fitpts);
    dT=diff(fitdels);
    ptwiseSlp=dS./dT;

    %% DK: MAY NEED TO FIX THIS PART
    % Detect jumps in the pointwise slope, defined as a change in slope 
    % within a certain percentage (thresh) of the initial slope
    thresh=0.5;
    slpJumps=(abs(mean(ptwiseSlp(1:3))-ptwiseSlp)>abs(mean(ptwiseSlp(1:3))...
        *(1+thresh)));
    % Fit the slope of the first 3 points 
%     preFitIdx=find(slpJumps,1,'first');
    preFitIdx=3;
    precurve=fit(fitdels(1:preFitIdx),fitpts(1:preFitIdx),'poly1');
    precoeffs=coeffvalues(precurve);
    preSlp=precoeffs(1);
    
    for jj=find(slpJumps)'
        % Shift all points after the slope jump to match the initial 
        % pointwise slope
        n2pis=round((preSlp*(fitdels(jj+1:end)-fitdels(jj))...
            -(fitpts(jj+1:end)-fitpts(jj)))./2./pi); 
        add2pis=n2pis*2*pi;  
        fitpts(jj+1:end)=fitpts(jj+1:end)+add2pis;
    end

    % Then, fit the curve to a linear function
    [curve,gof]=fit(fitdels,fitpts,'poly1');
    coeffs=coeffvalues(curve);
    w_BS(ii)=abs(coeffs(1));
    CIraw=confint(curve);
    CI_BS(ii,:)=abs(CIraw(:,1)); %95% confidence bounds of slope
    Rsq(ii)=gof.rsquare;
    
    % Plot each amplitude
    figure; plot(curve,fitdels,fitpts);
    title(['Bloch-Siegert fitting, F3 ampl = ' num2str(Ampls(ii))]);
    xlabel('Delay [s]');
    ylabel('Phase [rad]');
end

% Reorder CI_BS pairs so that the smaller value is in the 1st column
for ii=1:size(CI_BS,1)
    if CI_BS(ii,1)>CI_BS(ii,2)
        CI_BS(ii,:)=CI_BS(ii,[2,1]);
    end
end

% Finally, calculate omega_1 (in Hz) using the relevant equation
disp('Calculating the B1 in Hz at each amplitude setting...')

% Calculate offset difference between 1H and both the positive and negative
% X-nuclear frequencies
DoffPos=abs(HoffsetHz-XoffsetHz);
DoffNeg=abs(HoffsetHz+XoffsetHz);
DoffEff=(1./DoffNeg+1./DoffPos).^-1; %effective offset, as if from one field component

w_BS_Hz=w_BS./(2*pi);
% w1_Hz=sqrt(w_BS_Hz*2.*DoffNeg);
w1_Hz=sqrt(w_BS_Hz*2.*DoffEff);

CI_BS_Hz=CI_BS./(2*pi);
% CI_w1_Hz=sqrt(CI_BS_Hz*2.*repmat(abs(offsetHz),numel(Ampls),2));
CI_w1_Hz=sqrt(CI_BS_Hz*2.*repmat(DoffEff,numel(Ampls),2));
CI_w1_Hz_neg=w1_Hz-squeeze(CI_w1_Hz(:,1));
CI_w1_Hz_pos=squeeze(CI_w1_Hz(:,2))-w1_Hz;

% Identify points where the BS shift assumptions might be invalid
% BSthresh=5;
% BSvalid=(abs(w1_Hz*BSthresh)<=abs(offsetsHz));
% BSinvalid=(abs(w1_Hz*BSthresh)>abs(offsetsHz));
BSmaxangle=3; %maximum w_eff angle where Bloch-Siegert valid, in degrees
BSvalid=(abs(atand(w1_Hz./DoffPos))<=BSmaxangle);
% BSinvalid=(abs(atand(w1_Hz./DoffPos))>BSmaxangle);
if sum(~BSvalid)>0
    warning(['Some points in calibration curve may not fulfill '...
        'Bloch-Siegert assumptions! (i.e. w_eff angle > ' num2str(BSmaxangle) char(176)])
end

% Plot results
figure; errorbar(Ampls(BSvalid),w1_Hz(BSvalid),...
    CI_w1_Hz_neg(BSvalid),CI_w1_Hz_pos(BSvalid),'bo'); 
if sum(~BSvalid)>0
    hold on; errorbar(Ampls(~BSvalid),w1_Hz(~BSvalid),...
        CI_w1_Hz_neg(~BSvalid),CI_w1_Hz_pos(~BSvalid),'ro');
    leglbls={'Bloch-Siegert valid','Bloch-Siegert invalid'};
    legend(leglbls,'Location','southeast');
end
title('Measured B_1 vs Amplitude Setting'); 
xlabel('F3 ampl');
ylabel('B_1 (Hz)');

end