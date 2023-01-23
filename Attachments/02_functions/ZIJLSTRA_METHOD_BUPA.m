function ic = ZIJLSTRA_METHOD_BUPA(x,varargin)
%--------------------------------------------------------------------------
% ic = ZIJLSTRA_METHOD_BUPA(x,dim,fs,hp,lpl,lpu,fp,convertflag).
% Function for the estimation of initial contats (ICs) of the stride cycle 
% during walking activity, starting from the antero-posterior (AP) -with 
% reference to the walking motor task- trunk accelerations according with 
% peak detection method or zero-crossing method (for more details read 
% paragraph on fp input parameter) proposed by Zijlstra et al. in [1].
% In case peak detection method is used and no peaks are found prior to a
% certain zero-crossing (and after the previous one), then the relative 
% zero-crossing instant is identified as IC.
% The function, according to what requested by the user, is able to return 
% the ICs both in samples or in seconds. 
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
%        then iczijl treats as independent sequences data along dimension 
%        dim rather than along the first non-singleton dimension if dim is 
%        not given as input. 
%        The case of sequences in x made by only 1 point is not treated by 
%        the function (which will throw an error in case) but shall be 
%        managed by the calling entity: the method ,for how it is 
%        implemented, cannot even find a zero-crossing in such sequences. 
%        Moreover it is important to notice that even if works on sequences
%        of 2 points or longer, outputs may be no sense and/or highly 
%        inaccurate for sequences not containing at least few gait cycles).
%   - dim (optional): dimension of x along which the function will operate
%                     treating the elements as vectors.
%                     By default iczijl operates along the first
%                     non-singleton dimension of x. To use the default 
%                     behaviour for this optional parameter and still be 
%                     able to set the value for the next optional input 
%                     parameters, set the value of input parameter dim to 
%                     [] (empty array).
%   - fs (optional): sampling frequency of the input sequence x in Hz.
%                    To not insert fs and still be able to set the value 
%                    for the next optional input parameters, set the value 
%                    of input parameter fs to NaN.
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
%                    Indeed either in laboratory set up walking trials and 
%                    in real-world walking, the motor task may be affected 
%                    by acceleartions and decelerations. For example an 
%                    experimental set up may consist in acquiring signals 
%                    continuosly during a trial made up of many walking
%                    bouts (WBs), thus accelerations and decelerations are
%                    necessarely present, and in particular if the WBs are
%                    always of the same length scale there will be 
%                    oscillations always with a certain characteristic
%                    (average) period of the WB. The case of real-world
%                    gait is more various by this point of view but in
%                    general it is presumable to not to be interested in 
%                    short WBs, in which the task is highly not stationary 
%                    and likely to be affected by border effects that does 
%                    not relate with the health status of the subject, thus
%                    a suggested value range for the cut off frequency in 
%                    this case is [0.1,0.25] Hz, which does not lead to a 
%                    deterioration of the information related to gait 
%                    cycles (in which the fundamental armonic is around
%                    1.5-2.5 Hz in healthy subjects) and allows for 
%                    attenuation of oscillations in WBs of 2-4 s or longer.
%                    A reasonably good and robust -in terms of filter 
%                    stability- value for the order of this filter is 4.
%                    If fs is not given iczijl treats the cut off frequency
%                    as a normalized cut off frequency in range (0,1) where
%                    0 is the costant component and 1 is the Nyquist 
%                    frequency.
%                    If hp is not given (and "high" is not given as one of
%                    the entries of the adapt) the default behaviour is to 
%                    not to perform such HP filtering; if only a scalar is 
%                    given it is treated as if it is the cut off frequency
%                    and the related order used for the filter is 4. 
%                    To use the default behaviour -i.e. no such HP 
%                    filtering- and still be able to set the value of the 
%                    next input parameters, set the value of this parameter
%                    to [] (empty array).
%   - lpl (optional): proper 2-element vector (row or column) or scalar  
%                     containing respectively at least the cut off 
%                     frequency and eventually the order for the
%                     Butterworth anticausal low-pass (LP) filtering of the
%                     AP acceleration in x prior to zero-crossings (ZCs) 
%                     detection, which indeed allow to the fundamental 
%                     frequency of the walking motor task to emerge. 
%                     According to [1] the suggested value of cut off 
%                     frequency at this filtering stage is 2 Hz for healthy
%                     subjects, while a reasonably good and robust -in 
%                     terms of filter stability- value for the order of 
%                     this filter is 4.
%                     If fs is not given iczijl treats the cut off 
%                     frequency as a normalized cut off frequency in range
%                     (0,1) where 0 is the costant component and 1 is the  
%                     Nyquist frequency.
%                     If lpl is not given the default behaviour is to   
%                     perform such LP filtering with a cut off frequency
%                     of 2 Hz and order 4 as suggested in [1], in this case  
%                     insertion of fs input parameter is mandatory and no 
%                     more optional (see fs input parameter above); if only
%                     a scalar is given it is treated as if it is the cut 
%                     off frequency and the related order used for the 
%                     filter is 4. 
%                     To use the default behaviour and still be able to set
%                     the value of the next input parameters, set the value  
%                     of this parameter to [] (empty array).
%   - lpu (optional): proper 2-element vector (row or column) or scalar 
%                     containing respectively at least the cut off 
%                     frequency and eventually the order for the
%                     Butterworth anticausal LP filtering of the AP 
%                     acceleration in x prior to peaks detection, this 
%                     allows to attenuate high frequency noise that may be
%                     corrupting the gait signals. 
%                     According to [1] and established that armonics proper
%                     of the walking motor task do not exceed 25-35 Hz in 
%                     humans, the suggested value of cut off frequency at 
%                     this filtering stage is 20 Hz, while a reasonably 
%                     good and robust -in terms of filter stability- value 
%                     for the order of this filter is 4. 
%                     If fs is not given iczijl treats the cut off 
%                     frequency as a normalized cut off frequency in range 
%                     (0,1) where 0 is the costant component and 1 is the 
%                     Nyquist frequency.
%                     If lpu is not given the default behaviour is to not  
%                     to perform such LP filtering; if only a scalar is 
%                     given it is treated as if it is the cut off frequency
%                     and the related order used for the filter is 4. 
%                     To use the default behaviour -i.e. no such LP 
%                     filtering- and still be able to set the value of the 
%                     next input parameters, set the value of this 
%                     parameter to [] (empty array).
%   - fp (optional): proper 2-element vector (row or column) or scalar 
%                    containing respectively at least 'MinPeakDistance' and
%                    eventually 'MinPeakHeight' parameters in a.u. used as
%                    findpeak built-in Matlab function parameters during 
%                    the peaks detection step.
%                    If 'MinPeakDistance' and 'MinPeakHeight' are not given
%                    the default behaviour is to not perform peaks 
%                    detection and ICs are identified directly as the
%                    zero-crossing instants (zero-crossing method is 
%                    performed instead of peak detection method, for more 
%                    see refernce [1], to perform peak detection method as 
%                    described in [1] simply insert a value of 
%                    'MinPeakDistance' of 0). 
%                    To use the default behaviour and still be able to set 
%                    the value of the next input parameters, set the value
%                    of this parameter to [] (empty array).
%   - convertflag (optional): string scalar or character vector, specifies
%                             whether to return ICs instant estimates in 
%                             samples (defaul behaviour) if "samp" or in
%                             seconds if "time". "time" can be inserted
%                             only if sampling frequency fs is given.
                    
% OUTPUT ARGUMENTS:
%   - ic: output ICs instants. It will have the same size of x except that 
%         along the operating dimension, which size will be equal to the 
%         number of maximum ICs detected in any of the sequences in x; for
%         any of the sequences in which less ICs are detected then the
%         remaining slots along this dimension will be filled with NaN.
%         ICs instants can be returned either in samples or in seconds
%         according to what requested from the user.

% Ref.:
%   - [1]: Wiebren Zijlstra, At L. Hof. "Assessment of spatio-temporal gait
%          parameters from trunk accelerations during human walking", Gait 
%          Posture 2003,18:1-10.

% Author:       Burçak Carlo Pasqua   
% Last update:  20/06/2022     
%--------------------------------------------------------------------------

%% assignments and checks

narginchk(1,8);

% dim
shppre = size(x);  % shape of original x
dimg = false;  % flag for dim parameter given
if nargin >= 2
    dim = varargin{1};
    try
        if shppre(dim) == 1
            error(['x2 size along operating dimension is equal to 1, ' ...
                'size along the operating dimension shall be greater ' ...
                'than 1 as the method cannot be applied to sequences ' ...
                'of 1 point'])
        end
        if ~isempty(dim)
            dimg = true;
        end
    catch
        error(['dim must be a positive integer and must not exceed ' ...
            'the number of dimensions of x'])
    end    
end

% x
if ~isnumeric(x)
    error('invalid input type for x')
end
if any(isnan(x)) 
    warning('be warned that x given as input contains some NaN values')
end
if dimg && dim ~= 1
    x = shiftdim(x,dim-1);
elseif shppre(1) == 1  % shifting the leading non-singleton dimension to 
                       % the fore if there is a leading dimension of size 1
        [x,dim] = shiftdim(x);  
        dim = dim+1;
else  % case in which operating dimension is the first
    dim = 1;
end
shp = size(x);  % shape of adjusted x on which to work

% fs
fsg = false;
if nargin >= 3
    fs = varargin{2};
    if ~isnumeric(fs) || ~(fs>0) % (last logical operation also includes 
                                 % check on fs being NaN )
        error(['unacceptable fs parameter, it shall be a number ' ...
            'greater than 0']);
    else
        fsg = true;
    end
end

% hp
hpg = false;
if nargin >= 4
    hp = varargin{3};
    if ~isnumeric(hp) || (~all(hp>0) || ~isvector(hp))  && ~isempty(hp)
        error(['unacceptable hp parameter, it shall be a vector with ' ...
            'elements greater than 0 or empty'])
    elseif  ~isempty(hp)
        if ~fsg && hp(1) >= 1
            error('normalized cut off frequency shall be in (0,1) range')
        elseif fsg && hp(1) >= fs/2
            error('cut off frequency shall be in (0,fs/2) range')
        else
            hpg = true;
            switch length(hp)            
                case 2
                    ch = hp(1);  % cut off freq. of the HP filtering
                    nh = hp(2);  % filter order of the HP filtering
                case 1
                    ch = hp;  % cut off freq. of the HP filtering
                    nh = 4;  % filter order of the HP filtering
                otherwise
                    error('hp shall contain at most 2 elements')
            end
            if fsg
                ch = ch/(fs/2);
            end
        end
    end
end

% lpl
if nargin >= 5 && ~isempty(varargin{4})
    lpl = varargin{4};
    if ~isnumeric(lpl) || (~all(lpl>0) || ~isvector(lpl)) && ~isempty(lpl)
        error(['unacceptable lpl parameter, it shall be a vector with ' ...
            'elements greater than 0 or empty'])
    elseif  ~isempty(lpl)
        if ~fsg && lpl(1) >= 1
            error('normalized cut off frequency shall be in (0,1) range')
        elseif fsg && lpl(1) >= fs/2
            error('cut off frequency shall be in (0,fs/2) range')
        else
            switch length(lpl)            
                case 2
                    cll = lpl(1);  % cut off freq. of the LP filtering
                    nll = lpl(2);  % filter order of the LP filtering
                case 1
                    cll = lpl;  % cut off freq. of the LP filtering
                    nll = 4;  % filter order of the LP filtering
                otherwise
                    error('lpl shall contain at most 2 elements')
            end
            if fsg
                cll = cll/(fs/2);
            end
        end
    end
else  % default behaviour
    cll = 2/(fs/2);
    nll = 4;
end

% lpu
lpug = false;
if nargin >= 6
    lpu = varargin{5};
    if ~isnumeric(lpu) || (~all(lpu>0) || ~isvector(lpu)) && ~isempty(lpu)
        error(['unacceptable lpu parameter, it shall be a vector with ' ...
            'elements greater than 0 or empty'])
    elseif  ~isempty(lpu)
        if ~fsg && lpu(1) >= 1
            error('normalized cut off frequency shall be in (0,1) range')
        elseif fsg && lpu(1) >= fs/2
            error('cut off frequency shall be in (0,fs/2) range')
        else
            lpug = true;
            switch length(lpu)            
                case 2
                    clu = lpu(1);  % cut off freq. of the LP filtering
                    nlu = lpu(2);  % filter order of the LP filtering
                case 1
                    clu = lpu;  % cut off freq. of the LP filtering
                    nlu = 4;  % filter order of the LP filtering
                otherwise
                    error('lpu shall contain at most 2 elements')
            end
            if fsg
                clu = clu/(fs/2);
            end
        end
    end
end

% fp
fpg = false;
if nargin >= 7 
    fp = varargin{6};
    if ~isnumeric(fp)
        error('invalid input type for fp')
    end
    if ~isempty(fp)
        if fp(1) < 0
            error('min peak distance cannot be set to a negative value')
        end
        switch length(fp)            
            case 2
                mpd = fp(1);  % min peak distance
                mph = fp(2);  % min peak height
            case 1
                mpd = fp;  % min peak distance
                mph = -inf;  % min peak height
            otherwise
                error('fp shall contain at most 2 elements')
        end
        fpg = true;
    end
end

% convertflag
convert = false;
if nargin == 8
    if varargin{7} == "time"
        if ~fsg
            error('ICs can be returned in seconds only if fs is given')
        end
        convert = true; 
    elseif varargin{7} ~= "samp"
        error('invalid convertflag parameter')
    end
end

%% body

% offset removal: eventual offset on the measures does not carry
% information about the gait cycle events and if present may significantly
% alterate the detection of the zerocrossings, thus being computationally
% low expensive a stage of offset removal is always performed (offset
% may be present due to bias of probing sensor rather than non perfect 
% stationarity of the gait). Even if eventually (hp given by user) a HP
% filtering is performed, if the offset is very high a certain ratio of
% this could remain inevitably worsening the correct detection of the ZCs
x = x - mean(x);

% % un/comment from here <------------------------------------------------------
% % number of subplots definition and plot of original input sequences:
% % graphic representation is allowed only if a single sequence is given in
% % input
% if isvector(x)  
%     spe = 3;  % number of suplots for pattern enhanced show off at each
%               % processing stage. Minimum number is 3 in case only lower LP
%               % filtering is performed as one suplot for each status of 
%               % signal and one for the superimposition of all the stages
%     spd = 2;  % number of suplots for pattern detection (peaks and
%               % zero-crossings) show off. Minimum number is 2 in case
%               % zero-crossing method is performed, and becames 3 in case 
%               % that peak detection method is used
%     spe = spe+1*hpg;
%     spe = spe+1*lpug;
%     spd = spd+1*fpg;
%     figpe = figure;  % figure object for pattern enhancement show off
%     sgtitle('signal enhancemente across the different filtering stages'),
%     subplot(spe,1,1),plot(x),
%     title('original sequences passed to the function'),
%     xlim tight,grid minor
%     axspe = subplot(spe,1,spe);plot(x),
%     title('superimposition'),
%     xlim tight, grid minor,
%     hold on
%     figpd = figure;  % figure object for pattern enhancement show off
%     sgtitle('patterns detection and ICs identification'),
%     axspd = subplot(spd,1,spd);
%     title('superimposition'),
%     xlim tight, grid minor,
%     hold on
%     cpe = 2;  % counter on pattern enhancement show off plots
%     cpd = 1;  % counter on pattern detection show off plots
% end  % to here <--------------------------------------------------------------

% HP (costants and very slow trends components) filtering
if hpg
    [b,a] = butter(nh,ch,"high");
    x = filtfilt(b,a,x);
%     if isvector(x)  % un/comment from here <----------------------------------
%         % figure('Name','freq. response of HP filtering'),freqz(b,a)
%         figure(figpe),
%         subplot(spe,1,cpe),plot(x),
%         title('HP filtered'),
%         xlim tight,grid minor
%         cpe = cpe+1;
%         subplot(axspe),plot(x)
%     end  % to here <----------------------------------------------------------
end

% LPU filtering
if lpug
    [b,a] = butter(nlu,clu);
    x = filtfilt(b,a,x);
%     if isvector(x)  % un/comment from here <----------------------------------
%         % figure('Name','freq. response of upper LP filtering'),freqz(b,a)
%         figure(figpe),
%         subplot(spe,1,cpe),plot(x),
%         title('LPU filtered'),
%         xlim tight, grid minor
%         cpe = cpe+1;
%         subplot(axspe),plot(x)
%     end  % to here <----------------------------------------------------------
end

% peaks detection according to 'MinPeakDistance' and 'MinPeakHeight'
if fpg
    % as findpeaks function only can accept a vector sequence eventual non
    % vector arrays have to be vectorized, the issue that could happen due
    % to this implementation is that some local maxima toward the start of 
    % the sequences other than the first or toward the end of the sequences
    % other than the last could not be recognized as local maxima if 
    % eventually under the threshold 'MinPeakDistance' -if given- from 
    % peaks coming from other sequences than the proper after the 
    % vectorization  
    pksi = find(islocalmax(x,'MinSeparation',mpd) & x>mph);
%     if isvector(x)  % un/comment from here <----------------------------------
%         figure(figpd),
%         subplot(spd,1,cpd+1),plot(1:length(x),x,pksi,x(pksi),'go'),
%         title('peaks'),
%         xlim tight, grid minor
%         subplot(axspd),plot(1:length(x),x,pksi,x(pksi),'go')
%     end  % to here <----------------------------------------------------------
end

% LPL filtering
[b,a] = butter(nll,cll);
x = filtfilt(b,a,x);
% if isvector(x)  % un/comment from here <--------------------------------------
%     % figure('Name','freq. response of upper LP filtering'),freqz(b,a)
%     figure(figpe),
%     subplot(spe,1,cpe),plot(x),title('LPL filtered'),xlim tight, grid minor
%     subplot(axspe),plot(x)
% end  % to here <--------------------------------------------------------------

% zero-crossing detections
zc = x > 0;  
zc = find(zc.*~[reshape(zc(2:end,:),[shp(1)-1,shp(2:end)]);false([1, ...
    shp(2:end)])]);  % indeces of zero-crossings from positive to negative
                     % (i.e. ICs with zero-crossings method)
% if isvector(x)  % un/comment from here <--------------------------------------
%     figure(figpd),
%     subplot(spd,1,cpd),plot(1:length(x),x,zc,x(zc),'co'),
%     title('zero-crossings'),
%     xlim tight, grid minor
%     subplot(axspd),plot(1:length(x),x,zc,x(zc),'co')
% end  % to here <--------------------------------------------------------------

% ICs identification
if fpg  % ICs with peak detection method
    ic = repmat(pksi',[length(zc),1]) < zc;
    ic = sum(ic,2);  
    npf = find(diff(ic) == 0) + 1;  % indeces of zero-crossings for which 
                                    % no prior peak has been found 
    ic(ic == 0) = [];  % in case the first zero-crossing is not preceeded 
                       % by any peak -this could happen as at the   
                       % acquisition could start at any point of the gait
                       % cycle (eventually no prior peaks) or as the start
                       % of the motor task oscillations may occur without 
                       % having a gait cycle before which lead to a
                       % zero-crossing without peaks before
    ic = sort([unique(pksi(ic));zc(npf)]);  
else  % ICs with zero-crossing method
    ic = zc;  
end
% if isvector(x) && ~isempty(ic)  % un/comment from here <----------------------
%     subplot(axspd),xline(ic,'k--')
% end  % to here <--------------------------------------------------------------

% re-structuring output to expected shape (find function linearizes)
if ~isvector(x)
    edges = 1:shp(1):prod(shp)+1;  % enpoint linear indexes for each 
                                   % independent sequence in x
    icpseq = reshape(histcounts(ic,edges),[1,shp(2:end)]);  % # ICs found
                                                            % per each seq.
    rs = repmat((1:max(icpseq(:)))',[1,shp(2:end)]);  % re-structured shape
    rs = rs < icpseq+1;
    t = zeros(size(rs)); t(rs) = ic; t(~rs) = NaN;  % temporary side var
    ic = t; 
end
adj = 0:shp(1):shp(1)*(prod(shp(2:end))-1);  % subtraction subtrahend to
                                             % pass from the linear indeces
                                             % to proper indeces referred
                                             % to each sequence (find
                                             % buil-in Matlab function
                                             % linearizes the indeces)
adj = reshape(adj,[1,shp(2:end)]);
ic = ic-adj;

if convert  % conversion from samples to time in seconds
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

