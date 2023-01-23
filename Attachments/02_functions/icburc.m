function ic = icburc(x,fs,hp,lp,dim,options)
%--------------------------------------------------------------------------
% Function for the estimation of initial contats (ICs) of the stride cycle 
% during walking activity, starting from the antero-posterior (AP) -with 
% reference to the walking motor task- trunk accelerations.
% This method tries to overcome the main Achille's heel of the method 
% proposed in [1] which is the missed or erroneous detection of a 
% zero-crossing (ZC), as starting point to identify an IC, in case of 
% gait cycles in which abnormal accelerating patterns are superimposed to 
% the standard human accelerating pattern (which for example may occur in 
% the gait of a pathological subject). Indeed in this case the 
% superimposition of abnormal oscillations to the standard step-cycle 
% middle envelope of the gait may not result in an AP acceleration ZC, 
% which is at the basis of both the methods proposed in [1].
% To overcome the problem, the implemented method first segments the signal
% into different step cycles, the segmentation into different step cycles 
% is performed on the basis of a step-related middle envelope easily 
% obtainable through a low-pass (LP) filtering of the AP acceleration 
% signal with a cut off frequency slightly higher (125%) than the 
% step-related harmonic of the gait*¹, adaptively estimated for the gait 
% under analysis -choosing a LP filtering cut off frequency slightly higher
% than the step cycle related harmonic allows to best preserve periodic 
% fluctuations related to the step cycle while smoothing all the other 
% fluctuations, in this manner the step related pattern in highlighted-. 
% Then for each step-cycle identified the method searches the greatest 
% deceleration pattern -defined as the treat of the step cycle along which 
% there is the greatest monotonous deceleration- which is associated to an 
% IC happening. Finally the instant is identified as the median time point 
% between the instant of the acceleration peak of this pattern and the 
% instant in which, during this pattern, the signal assumes values lower 
% than the envelope; choosing the median point of these two allows for 
% greater robustness and highness of the method's performance in terms of 
% accuracy with respect to choosing the instant in which the peak occurred 
% (if during the decelerating the pattern no such straddle happens then 
% simply the peak is taken as IC instant).
% Moreover the fact that the method does not rely on ZCs instants makes it
% almost independent from signals' slow oscillations affliction not related
% to the gait cycle itself (which reflects the non-stationarity of the 
% motor task or eventual instrumentation's bias instability -even if
% usually this is not the case of accelerometry signals-) and in general 
% much more robust than the method proposed in [1] to the parameters choice
% of signals pre-processing stage, such as the cut off frequencies of the 
% denoising LP (see lp input parameter paragraph below) and high-pass (HP, 
% see hp input parameter paragraph below) filtering stages.
% As this method guarantees, assumed to correctly estimate the step related
% harmonic as the proper one and not as the stride related harmonic, to 
% find an IC for each step cycle, if it happens that during gaits with 
% impairment some step segmentation is lost (not identifiable on the step 
% related middle envelope), a correction mechanism is implemented so that
% an additional IC is searched even if the related step cycle is not 
% identifiable. In particular a check on the length of the step cycle is 
% done and if it lasts more than the 175% of the average step cycle length 
% (calculated as the inverse of the step cycle armonic) than the 2 greatest
% deceleration patterns are selected, instead of only the greatest, for 
% that actually segmented step cycle.
% ||| *¹: usually to find the step-related armonic in healthy subjects' 
%         gait is quite an easy task as it is sufficient to search for the 
%         spectral component carrying the maximum power in the signal, but 
%         in case a gait with pathology is analized this may not be true 
%         anymore (especially in case of significantly asymmetric gaits). 
%         Indeed in this case the periodicity of the walking activity may 
%         become coarsely the one of the gait cycle instead than the one of
%         the step cycle. 
%         Thus, in case of plausible high asymmetry of the gait due to 
%         different status of the two hemiparts of the subject from which  
%         the signal has been probed, the step-related harmonic of each of 
%         the sequences in input is found, instead of simply looking for  
%         the armonic carrying the higher contribute of power, firstly
%         performing a power spectral density (PSD) estimation of the 
%         sequence, then the harmonic carrying the maximum power contribute 
%         is chosen as plausible step-related harmonic candidate, finally a
%         check is made looking for the local maxima on the power spectral 
%         estimate (PSE) in a range of frequencies nearly around the double
%         of the candidate harmonic: if the ratio between the candidate 
%         harmonic contribute and its subsequent harmonic is greater than a
%         crisp thresold of 3.5 then the candidate is taken as the 
%         step-related harmonic (whereas it had been the gait cycle related 
%         harmonic a greater quantity of power would have been carried by 
%         the latter); otherwise the candidate is discarded and the 
%         step-related harmonic is set as the value equal to the subsequent
%         harmonic armonic in correspondence of the local maxima. 
%         The indication of possible significantly asymmetric gait can be 
%         given as input to the function (see Name-Value Argument inputs 
%         below) to perform this check and eventually the aforementioned
%         correction mechanism.
%         Even if this mechanism is used to enhance the chance to 
%         automatically find the right step-related harmonic, if the gait
%         under analysis is compromised at such a level in which one of the
%         limbs does not apport motility during the motor task but rather 
%         is carried behind -thus being the active part of the walking 
%         activity only the other limb-, it may happen that the proper gait
%         cycle related harmonic is pretended to be the step-related one 
%         since all the contribute to motion during the gait cycle is given
%         by one limb only. This causing the loss of nearly half of the ICs
%         occurring during that gait, indeed even if the decelerating 
%         pattern due to IC is still detectable similarly to those of the 
%         non impaired limb, the problem in this case is that only one IC 
%         for each pretended step-cycle is detected, in particular the one 
%         associated to the greatest monotonous deceleration (the 
%         correction mechanism does not enter in action as the estimated 
%         step-cycle armonic is for real the gait-cycle related harmonic, 
%         i.e. lost step cycles do not lead to condition of having step 
%         cycle lasting longer than the pretended characteristic step cycle
%         period of the gait).
% ||| PLEASE NOTICE that to correctly read the help and comments of the 
% function also N-dimensional array with only the last dimension's size not
% equal to 1 will be improperly referred to as "vectors" even if N is 
% greater than 2; while when referring to "N-dimensional arrays" it is 
% meant to be any other type of N-dimensional array (i.e. when any other 
% than the last dimension's size is not equal to 1). This is to more easily
% manage the cases of N-dimensional "vector" arrays by the moment that 
% Matlab treats the boolean indexing on these as on proper vectors in 2-D 
% (i.e mantaining the original orientation) instead of returning by default
% a column vector (as in the case of generic N-dimensional arrays) and to 
% easily escape the more complex sizes check that should be otherwise 
% implemented on these |||
% ||| If any error is found in the function please contact and warn me at 
% e-mail adress burcak.pasqua@gmail.com so that I can fix it for future 
% users |||


% INPUT ARGUMENTS:
%   - x: input sequence containing the AP accelerations. It can be a 
%        vector, a matrix or a N-dimensional array. If x is not a vector, 
%        then icburc treats as independent sequences data along dimension 
%        dim rather than along the first non-singleton dimension if dim is 
%        not given as input. 
%        The case of sequences in x made by less than 10 points is not 
%        treated by the function (which will throw an error in case) but 
%        shall be managed by the calling entity: the method ,for how it is 
%        implemented, likely would not find any reasonable result.
%        Moreover it is important to notice that even if works on sequences
%        of 11 points or longer, outputs may be no sense and/or highly 
%        inaccurate for sequences not containing at least few gait cycles).
%   - fs (optional): sampling frequency of the input sequence x in Hz.
%                    To not insert fs and still be able to set the value 
%                    for the next optional input parameters, set the value 
%                    of input parameter fs to [] (empty array).
%   - hp (optional): proper 2-element vector (row or column) or scalar 
%                    containing respectively at least the cut off frequency  
%                    and eventually the order for the Butterworth  
%                    anticausal high-pass (HP) filter of the AP 
%                    acceleration in x prior to the other steps. 
%                    The cut off frequency of this filtering stage should 
%                    be set according to low frequency oscillations which
%                    may affect the gait signals and on which there is no 
%                    interest, usually this is the case of low frequency 
%                    oscillations occurring due to not perfect stationarity 
%                    of the walk under analysis.
%                    Indeed either in laboratory set up walking trials or 
%                    in real-world walking, the motor task may be affected 
%                    by accelerations and decelerations. For example an 
%                    experimental set up may consist in acquiring signals 
%                    continuosly during a trial made up of many walking
%                    bouts (WBs), thus accelerations and decelerations are
%                    necessarely present, and in particular if the WBs are
%                    always of the same length scale there will be 
%                    oscillations always with a certain characteristic
%                    (average) period. The case of real-world gait is more 
%                    variegated by this point of view, but in general it is
%                    presumable to not to be interested in short WBs, in 
%                    which the task is highly not stationary and likely to 
%                    be affected by border effects that do not relate with
%                    the health status of the subject, thus a suggested
%                    value range for the cut off frequency for this 
%                    filtering stage is [0.1,0.25] Hz, which does not lead 
%                    to a deterioration of the information related to gait 
%                    cycles (in which the fundamental armonic is around
%                    1.5-2.5 Hz in healthy subjects) and allows for 
%                    attenuation of oscillations in WBs of 2-4 s or longer.
%                    A reasonably good and robust -in terms of filter 
%                    stability- value for the order of this filter is 4.
%                    icburc treats entered cut off frequency as non
%                    normalized so it should be bounded in range
%                    (0,Nyquist frequency), if a valid hp is entered fs is 
%                    mandatory and no longer optional.
%                    If hp is not given the default behaviour is to not to 
%                    perform such HP filtering; if only a scalar is given 
%                    it is treated as if it is the cut off frequency and
%                    the related order used for the filter is 4. 
%                    To use the default behaviour -i.e. no such HP 
%                    filtering- and still be able to set the value of the 
%                    next input parameters, set the value of this parameter
%                    to [] (empty array).
%   - lp (optional): proper 2-element vector (row or column) or scalar 
%                    containing respectively at least the cut off frequency
%                    and eventually the order for the Butterworth 
%                    anticausal LP filtering of the AP acceleration in x 
%                    prior to other steps in order to attenuate eventual 
%                    high frequency noise corrupting the signals. 
%                    According to [1] and established that armonics proper
%                    of the walking motor task do not exceed 25-35 Hz in 
%                    humans, the suggested value of cut off frequency at 
%                    this filtering stage is 20 Hz, while a reasonably good
%                    and robust -in terms of filter stability- value for 
%                    the order of this filter is 4. 
%                    icburc treats entered cut off frequency as non
%                    normalized so it should be bounded in range
%                    (0,Nyquist frequency), if a valid lp is entered fs is 
%                    mandatory and no longer optional.
%                    If lp is not given the default behaviour is to not to 
%                    perform such LP filtering; if only a scalar is given
%                    it is treated as if it is the cut off frequency and 
%                    the related order used for the filter is 4. 
%                    To use the default behaviour -i.e. no such LP 
%                    filtering- and still be able to set the value of the 
%                    next input parameters, set the value of this parameter
%                    to [] (empty array).
%   - dim (optional): dimension of x along which the function will operate
%                     treating the elements as vectors.
%                     By default icburc operates along the first
%                     non-singleton dimension of x. To use the default 
%                     behaviour for this optional parameter and still be 
%                     able to set the value for the next optional input 
%                     parameters, set the value of input parameter dim to 
%                     [] (empty array).
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%  NAME-VALUE ARGUMENTS:
%  Specify optional pairs of arguments as Name1=Value1,...,NameN=ValueN,
%  where Name is the argument name and Value is the corresponding value. 
%  Name-Value arguments must appear after other arguments, but the order
%  of the pairs does not matter.
%  Before R2021a, use commas to separate each name and value, and enclose
%  Name in quotes.
%    - TimeConversion: string scalar or character vector, specifies whether
%                      to return ICs instant estimates in samples (default
%                      behaviour) if "off" or in seconds if "on".
%                      Conversion can be requested only if sampling 
%                      frequency fs is given.
%    - PSE: power spectra estimate (PSE) as PSD or power spectrum, 
%           associated to the input sequence from which to find the step 
%           cycle related harmonic.
%           It can be a vector, a matrix or a N-dimensional array. If it is
%           a vector, icburc uses this as common spectral estimate for each
%           of the independent sequences in x (eventually only one); while 
%           if it is not then all the sizes except the one along the 
%           operating dimension shall match those of x and icburc treats as
%           independent PSEs data along this dimension.
%           If this Name-Value argument is not given as input then the PSE
%           of x is obtained from x with the default behaviour of pwelch
%           built-in Matlab function.
%    - Frequency Ref: frequency reference points associated to each of the 
%                     PSEs entered by the user through PSE Name-Value input
%                     argument.
%                     It must be a vector of length equal to the size of 
%                     entered PSE Name-Value argument along the operating
%                     dimension and it is used as common frequency 
%                     reference points vector for each of the independent
%                     estimates.
%                     The function treats the frequency reference points as
%                     a non-normalized frequency reference in range 
%                     [0,Nyquist frequency], if a frequency reference is 
%                     entered fs is mandatory and no longer optional.
%                     The frequency reference can be accepted only if the 
%                     PSE Name-Value argument has been inserted by the 
%                     user; if only PSE is entered, then this is treated as 
%                     one sided PSE of a real signal, so the frequency 
%                     reference is created as a vector of N equidistant
%                     points in [0,Nyquist frequency] range, with N being
%                     the number of points of the given PSD estimates, also
%                     in this case fs is mandatory and no longer optional; 
%                     while if also PSE is not given the frequency 
%                     reference is obtained during the same call to pwelch 
%                     used to obtain the PSE (see previous Name-Value input
%                     argument).
%                     ----------------------------------------------------
%                     ||| PLEASE BE AWARE that in PSE is given but the 
%                     frequency reference is not, then the points in the 
%                     created frequency reference vector may not match the 
%                     proper ones if the PSE was obtained by the user on 
%                     not equidistant frequency points or with built-in
%                     Matlab functions having an odd NFFT (in this latter 
%                     case the 
%                     PSD estimate PSEs refer to frequencies in the 
%                     half-open interval [0,fs/2) or [0,π), while the 
%                     frequency reference interval used by this function 
%                     has both the endpoints included as aforementioned). 
%                     So in this case an error on the adaptive estimate of
%                     the enveloping LP cut off frequency can be introduced
%                     at this level, THEREFORE ALSO THE INPUT OF THE PROPER 
%                     FREQUENCY REFERENCE IS ALWAYS RECOMMENDED WHEN PSE IS
%                     GIVEN BY THE USER |||
%    - Type: string scalar or character vector, specifies if the input
%            signals in x may be acquired from significantly asymmetric
%            gait if "asymm" or not (default behaviour) if "symm".

% OUTPUT ARGUMENTS:
%   - ic: output ICs instants. It will have the same size of x except that 
%         along the operating dimension, which size will be equal to the 
%         number of maximum ICs detected in any of the independent 
%         sequences in x; for any of the sequences in which less ICs are 
%         detected then the remaining slots along this dimension will be 
%         filled with NaN.
%         ICs instants can be returned either in samples or in seconds
%         according to what requested from the user.

% Ref.:
%   - [1]: Wiebren Zijlstra, At L. Hof. "Assessment of spatio-temporal gait
%          parameters from trunk accelerations during human walking", Gait 
%          Posture 2003,18:1-10.

% Author:       Burçak Carlo Pasqua   
% Last update:  24/06/2022     
%--------------------------------------------------------------------------

%% input arguments check and assignments

arguments
        x {mustBeNumeric,mustBeFinite}
        fs {mustBeScalarOrEmpty,mustBeNumeric,mustBeFinite} = []
        hp {mustBeNumeric,mustBeFinite} = []
        lp {mustBeNumeric,mustBeFinite} = []
        dim {mustBeScalarOrEmpty,mustBeInteger,mustBePositive} = []
        options.TimeConversion {mustBeMember(options.TimeConversion, ...
            ["on","off"])} = "off"
        options.PSE {mustBeNumeric,mustBeFinite} = []
        options.FrequencyRef {mustBeNumeric,mustBeFinite} = []
        options.Type {mustBeMember(options.Type,["asymm","symm"])} = "symm"
end

% dim
shppre = size(x);  % shape of original x
dimg = false;
if ~isempty(dim)
    dimg = true;
    if dim > ndims(x)
        error(['input parameter dim cannot exceed the number of ' ...
            'dimensions of input sequence x'])
    end
    if shppre(dim) < 10
        error(error(['x size along operating dimension shall be ' ...
            'greater than 10']))
    end
end

% x
if  dimg && dim ~= 1
    x = shiftdim(x,dim-1);
elseif shppre(1) == 1  % shifting the leading non-singleton dimension to 
                       % the fore if there is a leading dimension of size 1
        [x,dim] = shiftdim(x);  
        dim = dim+1;
else  % case in which operating dimension is the first
    dim = 1;
end
shp = size(x);  % shape of adjusted x on which to work

% hp
if ~isempty(hp)
    if isempty(fs)
        error(['a valid hp input parameter can be entered only if fs ' ...
            'is given too'])
    end
    if length(hp)==1
        hp(2) = 4;
    elseif length(hp) >= 3 || ~isvector(hp)
        error('unacceptable hp input parameter')
    end
    if any(hp <=0 ) || hp(1) >= fs/2 || mod(hp(2),1) ~= 0
        error(['unacceptable hp input parameter, order of filter ' ...
            'shall be a positive integer and the cut off frequency ' ...
            'shall be in range (0,Nyquist frequency)'])
    end
    hp(1) = hp(1)/(fs/2);
end

% lp
if ~isempty(lp)
    if isempty(fs)
        error(['a valid lp input parameter can be entered only if fs ' ...
            'is given too'])
    end
    if length(lp)==1
        lp(2) = 4;
    elseif length(lp) >= 3 || ~isvector(lp)
        error('unacceptable lp input parameter')
    end
    if any(lp <=0 ) || lp(1) >= fs/2 || mod(lp(2),1) ~= 0
        error(['unacceptable lp input parameter, order of filter ' ...
            'shall be a positive integer and the cut off frequency ' ...
            'shall be in range (0,Nyquist frequency)'])
    end
    lp(1) = lp(1)/(fs/2);
end

% TimeConversion
if options.TimeConversion == "on" && isempty(fs)
    error('fs is mandatory if time conversion is requested')
end

% PSE
if ~isempty(options.PSE)
    pxx = options.PSE;
    shppxxpre = size(pxx);
    if numel(pxx) == length(pxx)  % if it is a vector (length return the 
                                  % size of the biggest dimension of the 
                                  % input array)
        if shppxxpre(1) == 1
            pxx = shiftdim(pxx);  % standardization to column vectors
        end
    else   % case of noy being a vector        
        w2m = true(1,length(shppre));
        w2m(dim) = false;  % logical indeces indicating where to look to
                           % check matching between shapes
        if ~isequal(shppxxpre(w2m),shppre(w2m))
            % if shape do not accord that of x
            error(['shape of x and entered PSD shall accord if the ' ...
                'latter is not a vector or a N-dimensional "vector"'])
        end
        if dim ~= 1
            pxx = shiftdim(pxx,dim);
        end
    end
    f = 0:1/(shppxxpre(1)-1):1;
else
    if ~isempty(options.FrequencyRef)
        error(['The power spectral estimate frequency reference can ' ...
            'be inserted only if the power spectral estimate is ' ...
            'inserted too'])
    end
    [pxx,f] = pwelch(x-mean(x)); f= f/pi;  % offset removal prior to PSD 
                                           % estimation as if it not
                                           % priorly done by the user and 
                                           % an offset is present in the 
                                           % signal, then this ievitably 
                                           % affects the estimate
    if ~isvector(x) && ~ismatrix(x)
    % if is was a N-dimensional array reshaping is done as pwelch matrixes
    % N-dimensional array inputs
        pxx = reshape(pxx,[length(pxx(:,1)),shp(2:end)]);
    end
end
shppxx = size(pxx);

% FrequencyRef
if ~isempty(options.FrequencyRef)
    if isempty(fs)
        error('frequency reference can be accepted only if fs is given')
    end
    f = options.FrequencyRef;
    if any(f<0 | f>(fs/2))
        error(['unnaceptable frequency refernce, it should be bounded ' ...
            'between constant component frequency and Nyquist frequency'])
    end
    if numel(f) == length(f)  % if it is a vector
        if numel(f) ~= shppxxpre(dim)
            error (['length of the frequency reference shall match ' ...
                'the size of the power spectral estimate along the ' ...
                'operating dimension'])
        end
        shpf = size(f);
        if shpf(1) == 1
            f = shiftdim(f);  % standardization to column vector
        end
    else
        error(['the power spectral estimate frequency reference shall' ...
            ' be a vector or a N-dimensional "vector"'])
    end
    f = f/(fs/2);
end

%% body

% LP filtering
if ~isempty(lp)
    [b,a] = butter(lp(2),lp(1));
    x = filtfilt(b,a,x);
end

% HP filtering
if ~isempty(hp)
    [b,a] = butter(hp(2),hp(1),"high");
    x = filtfilt(b,a,x);
end

% step-cycle (normalized) armonic estimation
[~,sca] = max(pxx);
if options.Type == "asymm"
    % check if what found is the step-related armonic and if not correct it
    chk = zeros(shppxx);
    inds = 1:shppxx(1); inds = repmat(inds,[1,shppxx(2:end)]);
    inds = inds > 0.8*2*sca & inds < 1.2*2*sca;
    chk(inds) = pxx(inds);
    [~,sca2] = max(chk);
    adj4p = reshape(0:shppxx(1):prod(shppxx)-1,[1,shppxx(2:end)]);
    correction = pxx(sca+adj4p)./pxx(sca2+adj4p) <= 3.5;
    sca(correction) = sca2(correction);
end
if all(sca==sca(1))  % efficency check
    sca = sca(1);
end
co = f(sca)*1.25;  % step related frequency of the gait of each independent 
                   % input sequence in x

% step-related mean envelope extraction
if numel(co) ~= 1
    env = zeros(shppxx);
    for i = 1:numel(co)
        [b,a] = butter(4,co(i));
        env(:,i) = filtfilt(b,a,x(:,i));
    end
else
    [b,a] = butter(4,co);
    env = filtfilt(b,a,x);
end

% segmentation into step cycles (mean envelope based)
segm = islocalmin(env);
segm(1,:) = true;
segm = cumsum(segm);  % each step cycle of a sequence has a different 
                      % incremental number associated

% ICs identification
ic = nan(shp);
s = 1;  % sequence counter
sc = 1;  % step cycle counter
i = 1;  % ICs counter
scp = 2./f(sca);  % average step cycle period (in samples)
while s <= prod(shp(2:end))
    if length(scp)~=1
    % average step cycle duration (in samples) of current sequence 
        cscp = scp(s);  
    else
    % average step cycle duration (in samples) of all sequences in x
        cscp = scp;  
    end
    nsc = max(segm(:,s));  % number of step cycles for current sequence
    while sc <= nsc
        curr = segm(:,s)==sc;  % current step cycle indices in current seq. 
        dt = sum(curr);  % duration (in samples) of current step cycle
        if sc > 1 && sc < nsc || dt > cscp*0.6
        % id the step cycle is the first or last of the sequence, ICs are 
        % searched only if the duration is longer than 60% of the average
        % step cycle duration for that sequence (to avoid uncurrect 
        % detection in case for that step cycle the IC occurred out of the 
        % acquisition window)
            currsc = nan(shp(1),1); 
            currsc(curr) = x(curr,s);  % current step cycle
            upper = find(islocalmax(currsc));
            lower = find(islocalmin(currsc));
            if ~isempty(upper) && ~isempty(lower)
            % if not at least a pair of maximum and minimum is found for
            % the current step cycle then no IC is identified for the
            % current step cycle (another idea coul be to set it equal to 
            % the median between the maximum of current step cycle envelope
            % and its end -conceptually similar to zero crossing method 
            % proposed in [1]-)
                if upper(1) > lower(1)
                % as we are interested in deceleration patterns, in each 
                % step cycle only variation patterns starting with a local  
                % maxima are kept
                    lower(1) = [];
                end
                if length(upper) > length(lower)
                % and only variation patterns ending with a local minima 
                % are kept
                    upper(end) = [];
                end
                artificious = zeros(shp(1),1);  % artificious vector of 
                                                % zeros in which only 
                                                % entries in correspondence
                                                % of local maxima of 
                                                % current step cycle are
                                                % set equal to the
                                                % deceleration entity
                artificious(upper) = currsc(upper)-currsc(lower); 
                if dt >= cscp*1.75 && length(upper) > 1
                % correction mechanism if the current step cycle is too 
                % long
                    [~,accpeak] = maxk(artificious,2);
                    accpeak = sort(accpeak);  % sorting if they are not 
                                              % found in the right sample
                                              % order
                    % indeces in upper-lower pair
                    [~,indslu] = maxk(currsc(upper)-currsc(lower),2);  
                    indslu = sort(indslu);
                elseif ~isempty(upper)
                    [~,accpeak] = max(artificious);
                    % index in upper-lower pair
                    [~,indslu] = max(currsc(upper)-currsc(lower));  
                end
                if ~isempty(upper)
                    for k = 1:length(accpeak)
                        pattern = x(accpeak(k):accpeak(k)+(lower( ...
                            indslu(k))-upper(indslu(k))),s) > env( ...
                            accpeak(k):accpeak(k)+(lower(indslu(k))- ...
                            upper(indslu(k))),s);
                        straddle = find(~pattern.*[pattern(1);pattern( ...
                            1:end-1)]);  % straddle detection (if any)
                        if isempty(straddle)
                        % case in which no straddle point is found
                            ic(i,s) = accpeak(k);
                        else
                        % case in which straddle point is found (below 
                        % index equal to 1 is explicited as it may happen 
                        % that more than one straddling point is present, 
                        % in case the first is taken)
                            ic(i,s) = accpeak(k)+floor(straddle(1)/2);
                        end
                        i = i+1;
                    end
                end
            end
        end
        sc = sc+1;
    end
    s = s+1;
    sc = 1;
    i = 1;
end
%   keeping operating dimension size equal to the greatest number of ICs
%   found in any of the sequences
ic(end-min(sum(isnan(ic)),[],"all")+1:end,:) = [];  
ic = reshape(ic,[numel(ic(:,1)),shp(2:end)]);

if options.TimeConversion == "on" 
% conversion from samples to time in seconds
    ic = ic/fs;
end

% dimension reordering if needed
if dimg  % dimension reordering if needed
    dbadj = length(shppre)~=length(shp);  % # deleted dims when adjusted x 
                                          % to operate along columns
    ri = 1:length(shp);
    ri = circshift(ri,dim-1-dbadj);  % indeces for reordering
    ic = permute(ic,ri);
    if dbadj > 0
    % in case shifting the dimension of dim positions has led to a loss of
    % one or more dimensions during computation (but this loss does not 
    % have to remain since that operating dimension was not the last one),
    % restoring of these is done by the moment that we want to return an 
    % array with the same shape of pxx with only the operating dimension 
    % contracted to size equal to 1
        ic = reshape(ic,shppre);
    end
elseif dim ~= 1  % same result as the previous conditional statement but 
                 % more fast in case user did not input any operating 
                 % dimension
    ic = shiftdim(ic,-(dim-1));  
end