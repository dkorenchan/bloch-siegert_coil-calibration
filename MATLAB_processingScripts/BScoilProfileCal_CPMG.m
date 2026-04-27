 function [w1_Hz,offsetsHz,Rsq,correctedAmpls,ProfFitPars] = ...
     BScoilProfileCal_CPMG(targetB1,proffitflg,verbose) 

% Initially written 11/21/2024 by DK as BScoilProfileCal.m
% Updated 2/13/2025 by DK
% Ported to current file by DK on 3/26/26

% This function takes data from a Bloch-Siegert coil calibration (CPMG
% experient with an off-resonance pulse applied during every other spin-
% echo), and it calculates the nutation frequency as a function of the
% pulse offset using the Bloch-Siegert shift dependence of the signal phase
%
% OPTIONAL INPUTS: 
%   targetB1    --  Nominal B1 amplitude for all offsets (Hz)
%   proffitflg   --  Logical indicating whether to adjust corrected
%                   amplitudes using Lorentzian fitting or not
%   verbose     --  Logical indicating whether to plot individual fits to
%                   determine each Bloch-Siegert shift
%
% OUTPUTS:
%   w1_Hz           --  Vector of the measured B1 value at each saturation
%                       offset, in units of Hz
%   offsetsHz       --  Vector containing the offset frequencies from the 
%                       water peak, in units of Hz
%   Rsq             --  R-squared fitting statistic for each value of w1_Hz
%   correctedAmpls  --  (if input targetB1 is specified) Vector of corrected
%                       amplitudes from Ampls, based upon the measured B1 
%                       values and the desired targetB1
%   LorFitPars      --  (if input lorfitflg=true) Struct containing 
%                       Lorentzian fitting equation and fitted constant 
%                       values pertaining to it

if ~exist('verbose','var')
    verbose=false;
end

if ~exist('proffitflg','var')
    proffitflg=false;
end

% Load the data and the offsets in Hz
[rawdata,hdr,~,fname]=Read_Tecmag_DK;
disp('Reading in offsets and timings from table in header...')
tableVals=Read_Tecmag_Header_Table(fname,{'o1_0:2','rf0:2'});
HoffsetHz=abs(hdr.ob_freq(1)*1e6);
offsetsHz=tableVals{1}';
Ampls=tableVals{2}';

disp('Reading in Bloch-Siegert pulse duration and number of loops from header...')
[~,seqPars]=Read_TNT_seqPars(fname,true);
BStime = seqPars.Dashboard.Sequence.BSpw;
numechoes = seqPars.Dashboard.Sequence.Nloops;

% % Check that rawdata has dimensions (spectral,offset,TE)
% if size(rawdata,2)~=length(offsetsHz) && size(rawdata,2)==length(delaysS)
%     rawdata=permute(rawdata,[1,3,2]);
% end

% See if not all offsets were measured
if size(rawdata,3)<length(offsetsHz)
    offsetsHz=offsetsHz(1:size(rawdata,3));
    Ampls=Ampls(1:size(rawdata,3));
end

% % Get the time axis data 
% echo_str = inputdlg({'What is the B-S pulse duration, in ms?','How many loops with a B-S pulse?'});
% BStime = str2double(echo_str{1})/1000;
% numechoes = str2double(echo_str{2});

% Double-check size of data: should be (npoints,Nechoes,Noffsets) 
if size(rawdata,1) ~= hdr.acq_points
    warning('Data was not formatted as (points) x (echoes) x (offsets)! Reshaping...')
    rawdata = reshape(rawdata,hdr.acq_points,numechoes,[]);
end

delaysS = BStime:BStime:(BStime*numechoes);

% Remove 0's at the end and the first point

data = rawdata;
index = find(abs(data(:,1,1))==0,1);
data(index:end,:,:) = [];

data(1,:,:) = [];

% Unwrap the phase of the points in each echo
echodata = squeeze(angle(mean(data,1)));
phaseRad = unwrap(echodata,[],1);


% Then go through all offsets, phase-unwrap so that we get a good linear 
% change in phase, and fit to obtain the Bloch-Siegert shift (in rad/s)
disp('Fitting to determine the Bloch-Siegert shift at each offset...')

w_BS=zeros(size(phaseRad,2),1);
CI_BS=zeros(size(phaseRad,2),2);
Rsq=zeros(size(phaseRad,2),1);

% figure; 
% nsp=ceil(sqrt(size(phaseRad,1)));
fitdels=delaysS';

for ii=1:size(phaseRad,2)
    fitpts=unwrap(squeeze(phaseRad(:,ii)));
%     % Keep only the points where all FID pts contributing to the phase are
%     % above the SNR threshold
%     SNRthresh=3; %SNR threshold
%     FIDsignal=squeeze(abs(rawdata(pts,ii,:)));
%     FIDnoise=std(squeeze(rawdata(end-42:end-32,ii,:)),0,1); %last 10 points of each FID
%     SNRtPass=logical(prod(FIDsignal>SNRthresh*repmat(FIDnoise,length(pts),1),1));
%     fitpts=fitpts(SNRtPass);
%     fitdels=delaysS(SNRtPass);

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
%     fitpts=squeeze(phaseRad(:,ii));
%     
%     % Check whether more unwrapping is needed, looking for (1) sign change
%     % in slope; and (2) magnitude of phase change
%     dS=diff(fitpts);
%     dSshift=circshift(dS,1);
%     dSshift(1)=dS(1);
%     dSchange=(dS.*dSshift<0);
%     dSmagn=abs(dS)>pi/2;
%     dSeval=dSchange&dSmagn;
%     while sum(dSeval) > 0
%         % Find where to change, then change it
%         changeIdx=find(dSchange,1)+1; %need +1 since dS is 1 fewer length than fitpts!
%         phaseAddSign=(dS(1)>0)*2-1;
%         fitpts(changeIdx)=fitpts(changeIdx)+phaseAddSign*2*pi;
%         % Then recalculate dS, dSshift, and dSchange
%         dS=diff(fitpts);
%         dSshift=circshift(dS,1);
%         dSshift(1)=dS(1);
%         dSchange=(dS.*dSshift<0);
%         dSmagn=abs(dS)>pi;
%         dSeval=dSchange&dSmagn;        
%     end

    % Fit phase across each echo train to a linear function
    [curve,gof]=fit(fitdels,fitpts,'poly1');
    coeffs=coeffvalues(curve);
    w_BS(ii)=abs(coeffs(1));
    CIraw=confint(curve);
    CI_BS(ii,:)=abs(CIraw(:,1)); %95% confidence bounds of slope
    Rsq(ii)=gof.rsquare;
    
%     % Plot the fitting for each offset
% %     subplot(nsp,nsp,ii); 
    if verbose
        figure; plot(curve,fitdels,fitpts);
        title(['Bloch-Siegert fitting, offset = ' num2str(offsetsHz(ii))]);
        xlabel('Delay [s]');
        ylabel('Phase [rad]');    
    end
end

% Reorder CI_BS pairs so that the smaller value is in the 1st column
for ii=1:size(CI_BS,1)
    if CI_BS(ii,1)>CI_BS(ii,2)
        CI_BS(ii,:)=CI_BS(ii,[2,1]);
    end
end

% Finally, calculate omega_1 (in Hz) using the relevant equation
disp('Calculating the B1 amplitude at each offset...')
w_BS_Hz=w_BS./(2*pi);

% Calculate offset difference between 1H and both the positive and negative
% frequencies going through the coil
DoffPos=abs(offsetsHz);
DoffNeg=abs(HoffsetHz*2+offsetsHz);
DoffEff=(1./DoffNeg+1./DoffPos).^-1; %effective offset, as if from one field component

% w1_Hz=sqrt(w_BS_Hz*2.*abs(offsetsHz));
w1_Hz=sqrt(w_BS_Hz*2.*DoffEff);

CI_BS_Hz=CI_BS./(2*pi);
% CI_w1_Hz=sqrt(CI_BS_Hz*2.*repmat(abs(offsetsHz),1,2));
CI_w1_Hz=sqrt(CI_BS_Hz*2.*repmat(DoffEff,1,2));
CI_w1_Hz_neg=w1_Hz-squeeze(CI_w1_Hz(:,1));
CI_w1_Hz_pos=squeeze(CI_w1_Hz(:,2))-w1_Hz;

% Identify points where the BS shift assumptions might be invalid
% BSthresh=5;
% BSvalid=(abs(w1_Hz*BSthresh)<=abs(offsetsHz));
% BSinvalid=(abs(w1_Hz*BSthresh)>abs(offsetsHz));
BSmaxangle=3; %maximum w_eff angle where Bloch-Siegert valid, in degrees
BSvalid=(abs(atand(w1_Hz./offsetsHz))<=BSmaxangle);
% BSinvalid=(abs(atand(w1_Hz./offsetsHz))>BSmaxangle);
% Fill in all the points inbetween the Bloch-Siegert invalid points
disp('Removing B-S valid points inbetween invalid ones...')
BSvalid(offsetsHz<max(offsetsHz(~BSvalid))&offsetsHz>min(offsetsHz(~BSvalid)))=false;
if sum(~BSvalid)>0
    warning(['Some points in calibration curve may not fulfill '...
        'Bloch-Siegert assumptions! (i.e. w_eff angle > ' num2str(BSmaxangle) char(176)])
end

% Plot results
offsets_kHz=offsetsHz/1000;
maxOffIdx=find(offsets_kHz==max(offsets_kHz));
minOffIdx=find(offsets_kHz==min(offsets_kHz));

figure; errorbar(offsets_kHz(BSvalid),w1_Hz(BSvalid),...
    CI_w1_Hz_neg(BSvalid),CI_w1_Hz_pos(BSvalid),'bo'); 
hold on; errorbar(offsets_kHz(~BSvalid),w1_Hz(~BSvalid),...
    CI_w1_Hz_neg(~BSvalid),CI_w1_Hz_pos(~BSvalid),'ro'); 
title('B_1 Amplitude vs Tx Offset'); 
xlabel('Offset Frequency (kHz)');
ylabel('B_1 (Hz)');

% If user specified a nominal B1 amplitude, plot on graph, then make other
% plots
if exist('targetB1','var')
    plot(offsets_kHz([minOffIdx,maxOffIdx]),[1,1]*targetB1,'r--');
    leglbls={'Bloch-Siegert valid','Bloch-Siegert invalid','Target'};
else
    leglbls={'Bloch-Siegert valid','Bloch-Siegert invalid'};
end
legend(leglbls,'Location','southeast');

if exist('targetB1','var')
    % B1 amplitude relative to target
    w1_rel=w1_Hz./targetB1;
    CI_w1_rel=CI_w1_Hz./targetB1;
    CI_w1_rel_neg=w1_rel-squeeze(CI_w1_rel(:,1));
    CI_w1_rel_pos=squeeze(CI_w1_rel(:,2))-w1_rel;
    
    figure; errorbar(offsets_kHz(BSvalid),w1_rel(BSvalid),...
        CI_w1_rel_neg(BSvalid),CI_w1_rel_pos(BSvalid),'bo'); 
    hold on; errorbar(offsets_kHz(~BSvalid),w1_rel(~BSvalid),...
        CI_w1_rel_neg(~BSvalid),CI_w1_rel_pos(~BSvalid),'ro'); 
    plot(offsets_kHz([minOffIdx,maxOffIdx]),[1,1],'r--');
    title('B_1 Relative to Target vs Tx Offset'); 
    xlabel('Offset Frequency (kHz)');
    ylabel('B_{1,meas}/B_{1,target}');
    legend(leglbls,'Location','southeast');
    
    % Correct Ampls and issue warnings if values exceed 100 or are lower 
    % than 1
%     if nargin>1
%         % See if user specified to fit 1./Ampls to Lorentzian based upon
%         % lorfitflg; if not, do this by default
%         if exist('lorfitflg','var')
%             disp('Input lorfitflg not specified. Setting to true...')
%             proffitflg=true;
%         else
%             proffitflg=false;
%         end
        
        if proffitflg %fit to Lorentzian
            disp(['Correcting inputted off-resonance pulse amplitudes where '...
                'Bloch-Siegert is invalid using Lorentzian fitting...'])
        else %just calculate point-by-point
            disp('Correcting raw inputted off-resonance pulse amplitudes...')
        end
        
        correctedAmpls=Ampls./w1_rel;

        invAmpls=1./correctedAmpls;
%         % Lorentzian fitting function, start/upper/lower bounds 
%         % (try a Gaussian??)
%         fcn='a*real((b+1i*(x-c))/(b^2+(x-c)^2))';%+d*x';
%         st=[1,40,15];%,.0001];
%         ub=[5,100,22];%,.001];
%         lb=[0.2,1,-1];%,-.001];
        % RC circuit response curve (like a Lorentzian, but not quite). Fit
        % using x = offsets_kHz and y = inverted amplitudes. B1 is 
        % proportional to current which is proportional to voltage, which 
        % can be described using this equation (a = B1max,  
        % b = tuned freq, c = Q, d = constant offset of curve)
        Hoffset_kHz=HoffsetHz/1000;
        offsets_kHzAbs=offsets_kHz+Hoffset_kHz;
        fcn='a/sqrt(1+b^2*(x/c-c/x)^2)+d';
        st=[.02,4,Hoffset_kHz,0];
        ub=[1,200,HoffsetHz+30,0.5];
        lb=[0.001,0,Hoffset_kHz-30,-0.5];        
        % Exclude far offsets from fitting
        farOff=abs(offsets_kHz)>1000;

        % Fit to curve (only fit points where Bloch-Siegert assumption is 
        % valid!), then plot results of fitting
%         [profcurve,profgof]=fit(offsets_kHz(BSvalid&~farOff),invAmpls(BSvalid&~farOff),...
%             fcn,'StartPoint',st,'Upper',ub,'Lower',lb); 
%         figure; plot(profcurve,offsets_kHz(BSvalid),invAmpls(BSvalid));
%         hold on; scatter(offsets_kHz(~BSvalid),invAmpls(~BSvalid),'rx');
%         xlabel('Offset (kHz)'); ylabel('1/ampl')
        [profcurve,profgof]=fit(offsets_kHzAbs(BSvalid&~farOff),invAmpls(BSvalid&~farOff),...
            fcn,'StartPoint',st,'Upper',ub,'Lower',lb); 
        figure; plot(profcurve,offsets_kHzAbs(BSvalid),invAmpls(BSvalid));
        hold on; scatter(offsets_kHzAbs(~BSvalid),invAmpls(~BSvalid),'rx');
        xlabel('Frequency (kHz)'); ylabel('1/ampl')

        profcoeffs=coeffvalues(profcurve);
        
        if proffitflg
            title('Before profile correction');
%             lorvals=1./correctedAmpls;
        
            % Obtain correctedAmpls using best-fit function, and plot results    
%             lorvals(BSinvalid)=lorcoeffs(1)*abs((lorcoeffs(2)+...
%                 1i*(offsets_kHz(BSinvalid)-lorcoeffs(3)))...
%                 ./(lorcoeffs(2)^2+(offsets_kHz(BSinvalid)-lorcoeffs(3)).^2));%...
% %                 +lorcoeffs(4)*offsets_kHz(BSinvalid);
%             profvals=profcoeffs(1)*real((profcoeffs(2)+...
%                 1i*(offsets_kHz-profcoeffs(3)))...
%                 ./(profcoeffs(2)^2+(offsets_kHz-profcoeffs(3)).^2));%...
% %                 +lorcoeffs(4)*offsets_kHz(BSinvalid);
            profvals=profcoeffs(1)./sqrt(1+profcoeffs(2)^2*...
                (offsets_kHzAbs./profcoeffs(3)-profcoeffs(3)./offsets_kHzAbs).^2)+profcoeffs(4);
            figure; plot(profcurve,offsets_kHzAbs(BSvalid),profvals(BSvalid));
            hold on; scatter(offsets_kHzAbs(~BSvalid),profvals(~BSvalid),'rx');
            title({'After profile correction'});%,'(of B-S invalid values only)'});
            xlabel('Frequency (kHz)'); ylabel('1/ampl')
            
%             correctedAmpls=1./profvals;
            correctedAmpls=1./sqrt(profvals);
        end

        % Generate struct output containing relevant info for
        % profile fit
        ProfFitPars.equation=fcn;
        ProfFitPars.a=profcoeffs(1);
        ProfFitPars.b=profcoeffs(2);
        ProfFitPars.c=profcoeffs(3);
        ProfFitPars.d=profcoeffs(4);
        ProfFitPars.gof=profgof;
        ProfFitPars.x_units='kHz, absolute frequency';
            
        % Generate warnings if any corrected amplitudes are not practically
        % usable
        if sum(correctedAmpls>100)>0
            warning(['Corrected amplitudes include values >100! ' ...
                'You will need to DECREASE the F1 attenuation '...
                'and adjust these values accordingly.'])
        end
        if sum(correctedAmpls<1)>0
            warning(['Corrected amplitudes include values <1! ' ...
                'Consider INCREASING the F1 attenuation'...
                'and adjusting these values accordingly.'])
        end        
%     end

%     % B1 amplitude relative to target, dB power scale
%     w1_reldB=log10(w1_rel)*20;
%     figure; scatter(offsets_kHz(BSvalid),w1_reldB(BSvalid),'b'); 
%     hold on; scatter(offsets_kHz(~BSvalid),w1_reldB(~BSvalid),'r');
%     plot(offsets_kHz,zeros(size(w1_Hz)),'r--');    
%     title('dB Power Above Target vs Tx Offset'); 
%     xlabel('Offset Frequency (kHz)');
%     ylabel('dB power');   
%     legend(leglbls,'Location','southeast');
end
end