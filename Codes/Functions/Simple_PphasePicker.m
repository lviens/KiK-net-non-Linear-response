function [loc, snr_db] = Simple_PphasePicker(x, dt)
%   AN AUTOMATIC P-PHASE ARRIVAL TIME PICKER
%
%   Computes P-phase arrival time in windowed digital single-component
%   acceleration or broadband velocity record without requiring threshold
%   settings. Returns P-phase arrival time in second, and signal-to-noise
%   ratio in decibel. Input waveform must be an evenly spaced vector.
%
%   Syntax:
%      [loc, snr_db] = PphasePicker(x, dt, type, pflag, Tn, xi, nbins, o)
%
%   Input:
%            x = raw broadband velocity or acceleration data in
%                single-column format
%           dt = sampling interval in second (e.g., 0.005)

%   Output:
%          loc = P-phase arrival time in second
%       snr_db = signal-to-noise ratio in decibel
%
%   Example:
%          Let x be a single component strong-motion acceleration waveform
%          with 100 sample per second (dt = 0.01). The input for
%          P-phase picking will be
%
%          [loc, snr_db] = PphasePicker(x, 0.01);
%
%
%   IMPORTANT NOTE- 1: If sampling rate of input signal is lower than 100
%   samples-per-second, use Tn = 0.1 s instead of 0.01 s to avoid numerical
%   errors. If Tn is not specified, default values will be used based on
%   the sampling rate.
%
%   Kalkan, E. (2016). "An automatic P-phase arrival time picker", Bull. of
%   Seismol. Soc. of Am., 106, No. 3, doi: 10.1785/0120150111
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 16.0 $  $Date: 2017/02/17 18:14:00 $
%   Modified by Loic Viens Date: 2018/09/06 

if nargin == 2
    validateattributes(x,{'double'},{'real','finite','vector'}, ...
        'PphasePicker','X');
    validateattributes(dt,{'double'},{'real','finite','scalar'}, ...
        'PphasePicker','DT');
else
    error('Not enough inputs.  See help documentation.');
end


    if dt <= 0.01
        Tn = 0.01;
        nbins = round(2/dt);
    else
        Tn = 0.1;
        nbins = 200;
    end
    xi = 0.6;

o = 'to_peak';
x_d = detrend(x); % detrend waveform

switch o
    case {'to_peak'}
        ind_peak = find(abs(x_d) == max(abs(x_d)));
        xnew = x_d(1:ind_peak);
    otherwise
        xnew = x_d;
end

% Construct a fixed-base viscously damped SDF oscillator
omegan = 2*pi/Tn;           % natural frequency in radian/second
C = 2*xi*omegan;            % viscous damping term
K = omegan^2;               % stiffness term
y(:,1) = [0;0];             % response vector

% Solve second-order ordinary differential equation of motion
A = [0 1; -K -C]; Ae = expm(A*dt); AeB = A\(Ae-eye(2))*[0;1];
for k = 2:length(xnew); y(:,k) = Ae*y(:,k-1) + AeB*xnew(k); end

veloc = (y(2,:))';          % relative velocity of mass
Edi = 2*xi*omegan*veloc.^2; % integrand of viscous damping energy

% Apply histogram method
R = statelevel(Edi,nbins);
locs = find(Edi > R(1));
indx = find(xnew(1:locs(1)-1).*xnew(2:locs(1)) < 0); % get zero crossings
TF = isempty(indx);

% Update first onset
if TF == 0
    loc = indx(end)*dt;
else
    R = statelevel(Edi,ceil(nbins/2)); % 
    locs = find(Edi > R(1));
    indx = find(xnew(1:locs(1)-1).*xnew(2:locs(1)) < 0); % get zero crossings
    TF = isempty(indx);
    if TF == 0
        loc = indx(end)*dt; 
    else
        loc = -1; 
    end
end

% Compute SNR
if ~loc == -1
    snr_db = -1;
else
    snr_db = SNR(x,x(1:round(loc/dt)));
end
return

function [levels, histogram, bins] = statelevel(y,n)
ymax = max(y);
ymin = min(y)-eps;

% Compute Histogram
idx = ceil(n * (y-ymin)/(ymax-ymin));
idx = idx(idx>=1 & idx<=n);
histogram = zeros(n, 1);
for i=1:numel(idx)
    histogram(idx(i)) = histogram(idx(i)) + 1;
end

% Compute Center of Each Bin
ymin = min(y);
Ry = ymax-ymin;
dy = Ry/n;
bins = ymin + ((1:n)-0.5)*dy;

% Compute State Levels
iLowerRegion = find(histogram > 0, 1, 'first');
iUpperRegion = find(histogram > 0, 1, 'last');

iLow  = iLowerRegion(1);
iHigh = iUpperRegion(1);

% Define the lower and upper histogram regions halfway
% between the lowest and highest nonzero bins.
lLow  = iLow;
lHigh = iLow + floor((iHigh - iLow)/2);
uLow  = iLow + floor((iHigh - iLow)/2);
uHigh = iHigh;

% Upper and lower histograms
lHist = histogram(lLow:lHigh, 1);
uHist = histogram(uLow:uHigh, 1);

levels = zeros(1,2);
[~, iMax] = max(lHist(2:end));
[~, iMin] = max(uHist);
levels(1) = ymin + dy * (lLow + iMax(1) - 1.5);
levels(2) = ymin + dy * (uLow + iMin(1) - 1.5);

return

function [r] = SNR(signal,noise)
%  Compute signal-to-noise ratio 
aps = mean(signal.^2); % average power of signal
apn = mean(noise.^2);  % average power of noise
r = 10*log10(aps/apn); % signal-to-noise ratio in decibel
return