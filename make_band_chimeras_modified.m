%% Checking current program file directory (run this first)
%disp(pwd)
path = pwd;
%Changing directories (alter this section based on location of .wav files)
cd 'C:\Users\russe\Desktop\Purdue-First Year\Spring Semester\BME 51100\Cochlear Implant Project\Melodies'
ls  %Should see the .wav files after running this part
%% List part
pat = "2";
list = ls;
list_array = [];
% class(list)
for l = 1:length(list)
    tf = contains(list(l, :), pat);
    if tf == 1
        list_array = [list_array + " ", list(l, :)];
    end
end
% for l = 1:length(list_array)
%     disp(list_array(l))
% end

%% Script for Execution
for i = 1:2:length(list_array)
    make_band_chimeras(list_array(i), list_array(i+1), [1 8 16 32])
end
%% Function that synthesizes single/multi-band auditory chimeras
function make_band_chimeras(name1, name2, Nbands, do_play, do_plot)
% Synthesize and store single or multi-band auditory chimeras from pair of original signals.
% Usage: make_band_chimeras(name1, name2, Nbands, do_play, do_plot)
%	name1, name2	Original signals (WAV files), Mono or Stereo.
%				If the second signal is 'noise' or is not specified,
%				random noise with same power spectrum as NAME1 is used
%	Nbands		Number of frequency bands used for synthesis.
%				Can be scalar or vector, in which case one pair of
%				chimeras is synthesized for each element of Nbands.
%				By default Nbands = [1 8 32];
%	do_play		Set flag to play chimeras (default no play)
%	do_plot		Set flag to plot stimulus waveforms (default no plot)
%
%	The function creates two WAV files to store the chimeras
%	for each value of Nbands.  Files are named
%		'name1_env+name2_fts-nb#' and 'name2_env+name1_fts-nb#', 
%	where # is the number if bands.  
%	The files are stored in a folder specified by the Matlab global variable
%	'ChimeraFolder'.  By default, this is set to the current folder '.'
%
%	Copyright Bertrand Delgutte, 1999-2000
%
if nargin == 0, error('Give names of two WAV files'); end
if nargin < 3, Nbands = [1 8 32]; end
if nargin < 4, do_play = 1; end %Play the chimera sounds
if nargin < 5, do_plot = 0; end

refilter = 0;
if refilter, refilt_code = '_f2'; else refilt_code = ''; end

global ChimeraFolder;
if isempty(ChimeraFolder), ChimeraFolder = '.'; end

% read original waveforms
[orig1, Fs] = audioread(name1);
nchan = size(orig1,2);
if nargin >= 2 & ~strcmp(name2, 'noise'),
	[orig2, Fs2] = audioread(name2);
   if Fs ~= Fs2, error('Incompatible sampling rates'); end
   if size(orig2,2) ~= nchan, error('Inconsistent numbers of channels'); end
	if length(orig2) < length(orig1),
	    orig2 = [orig2; zeros(length(orig1)-length(orig2),nchan)];
	elseif length(orig2) > length(orig1),
	    orig1 = [orig1; zeros(length(orig2)-length(orig1),nchan)];
	end
else	% synthesize noise with same power spectrum as original
	name2 = 'noise';
	orig2 = psd_matched_noise(orig1);
   orig2 = orig2/max(abs(orig2(:)));
end

if do_play,	% play original sounds
	disp('Playing original sounds')
	tmp = [orig1; orig2];
	if strcmp(computer, 'SUN4'), sunsound(tmp, Fs); 
    else sound(tmp, Fs); end
    pause(ceil(length(tmp)/Fs))
end

if do_plot,	% plot waveforms of original signals

	t=[0:length(orig1)-1]*1000/Fs;

	figure(1)
	clf
	subplot(2,1,1)
	plot(t, orig1);
	title(sprintf('Original "%s" Sound', name1))
	subplot(2,1,2)
	plot(t, orig2)
	title(sprintf('Original "%s" Sound', name2))
	xlabel('Time (ms)')
	drawnow
end

Fmin = 80;	% lower frequency of filterbank in Hz
Fmax = .4*Fs;	% upper frequency of filterbank (.8 * Nyquist)

for nb = Nbands,

    % determine band cutoffs equally spaced on basilar membrane
    Fco = equal_xbm_bands(Fmin, Fmax, nb);	

    % compute multi-band chimeras
    [env1_fts2, env2_fts1] = multi_band_chimera(orig1, orig2, Fco, Fs, refilter);

    % normalize and save
    env1_fts2 = env1_fts2./max(abs(env1_fts2(:)));
    env2_fts1 = env2_fts1./max(abs(env2_fts1(:)));
    
    chimfileA = sprintf('%s_env+%s_fts-nb%d%s.wav', name1, name2, nb, refilt_code);
    chimfileB = sprintf('%s_env+%s_fts-nb%d%s.wav', name2, name1, nb, refilt_code);
    chimfileA_string = num2str(chimfileA);
    chimfileB_string = num2str(chimfileB);
    audiowrite([ChimeraFolder '\' chimfileA_string], env1_fts2, Fs, 'BitsPerSample', 16);
    audiowrite([ChimeraFolder '\' chimfileB_string], env2_fts1, Fs, 'BitsPerSample', 16);


    if do_play,		% play band chimeras
	   	disp(sprintf('Playing %d-band chimeras', nb))
    	tmp = [env1_fts2; env2_fts1];
    	if strcmp(computer, 'SUN4'), sunsound(tmp, Fs); 
      	else sound(tmp, Fs); end
      	pause(ceil(length(tmp)/Fs))
    end

    if do_plot,		% plot waveforms of band chimeras

       figure(2)
       clf
       
       subplot(2,1,1)
       plot(t, env1_fts2)
       xlabel('Time (ms)')
       title(sprintf('%d-band Chimera ("%s" envelope, "%s" fine structure)', ...
          nb, name1, name2))
       
       subplot(2,1,2)
       plot(t, env2_fts1)
       xlabel('Time (ms)')
       title(sprintf('%d-band Chimera ("%s" envelope, "%s" fine structure)', ...
          nb, name2, name1))
       
       drawnow
    end	% if do_plot

end	% loop over Nbands
end

function fco = equal_xbm_bands(fmin, fmax, N)
% EQUAL_XBM_BANDS - Divide frequency interval into N bands of equal width
% 		    along the human basilar membrane.
% Based on M.C. Liberman's cochlear frequency map for the cat 
% scaled to match human frequency range of hearing.
%	
% Usage: fco = equal_xbm_bands(fmin, fmax, N)
%	fmin	minimum frequency in Hz
%	fmax	maximum frequency in Hz
%	N	number of frequency bands
%	fco	Vector of band cutoff frequencies in Hz (size [1 X (N+1)])
%
%	Copyright Bertrand Delgutte, 1999-2000
%
xmin = cochlear_map(fmin, 20000);
xmax = cochlear_map(fmax, 20000);

dx = (xmax-xmin)/N;
x = xmin:dx:xmax;

fco = inv_cochlear_map(x, 20000);
end

function x = cochlear_map(f, Fmax)
% INV_COCHLEAR_MAP - Convert frequency to distance along the basilar membrane
% 		     using M.C. Liberman's cochlear frequency map for the cat.
%
% Usage: x = cochlear_map(f, Fmax)
%	f		Frequency in Hz
%	Fmax	Maximum frequency represented on the basilar membrane in Hz.
%			By default, this is 57 kHz, the value for the cat.
%			Setting Fmax to 20,000 Hz gives a map appropriate for the human cochlea.
%	x		Percent distance from apex of basilar membrane
%
%	Copyright Bertrand Delgutte, 1999-2000
%
% SEE ALSO: inv_cochlear_map			

if nargin > 1,
   	f = f * inv_cochlear_map(100)/Fmax;
end

x = log10(f/456 + .8)/.021;
end

function f = inv_cochlear_map(x, Fmax)
% INV_COCHLEAR_MAP - Convert distance along the basilar membrane to frequency
% 		     using M.C. Liberman's cochlear frequency map for the cat.
%
% Usage: f = inv_cochlear_map(x, Fmax)
%	x		Percent distance from apex of basilar membrane
%	Fmax	Maximum frequency represented on the basilar membrane in Hz.
%			By default, this is 57 kHz, the value for the cat.
%			Setting Fmax to 20 kHz gives a map appropriate for the human cochlea.
%	f		Frequency in Hz
%
%	Copyright Bertrand Delgutte, 1999-2000
%
% SEE ALSO: cochlear_map			

f = 456 * 10.^(.021 *x) - 364.8;

if nargin > 1,
   	f = f * Fmax/inv_cochlear_map(100);
end
end

function [e1_fs2, e2_fs1] = multi_band_chimera(x1, x2, Fco, Fs, refilter)
% MULTI_BAND_CHIMERA - Synthesize pair of multi-band "auditory chimeras"
% by dividing each signal into frequency bands, then interchanging 
% envelope and fine structure in each band using Hilbert transforms.
%
% Usage:  [e1_fs2, e2_fs1] = multi_band_chimera(x1, x2, Fco, Fs, refilter)
%	x1, x2	original signals
%	Fco	    band cutoff frequencies (Nbands+1), or filter bank
%	Fs	    sampling rate in Hz (default 1)
%   refilter  set to 1 to filter again after exchange operation (default 0)
% 	e1_fs2  chimera with envelope of x1, fine structure of x2 in each band
% 	e2_fs1  chimera with envelope of x2, fine structure of x1 in each band
%
%	Copyright Bertrand Delgutte, 1999-2000
%
if nargin < 3, error('Specify original signals and cutoff frequencies'); end
if nargin < 4, Fs = 1; end
if nargin < 5, refilter = 0; end

if min(size(x1)) == 1, x1 = x1(:); end	% make column vectors
if min(size(x2)) == 1, x2 = x2(:); end
nchan = size(x1, 2);

% pad with zeroes to match lengths
if length(x1) < length(x2),
	x1 = [x1; zeros(length(x2)-length(x1),nchan)];
elseif length(x2) < length(x1),
	x2 = [x2; zeros(length(x1)-length(x2),nchan)];
end
   
% Because the Hilbert transform and the filter bank are both linear operations,
% they commute and associate.  We create a bank of complex FIR filters whose
% real and imaginary parts are in quadrature (cos and sin).  This complex filter 
% is directly applied the original signals. The Hilbert envelope in each band
% is the absolute value of the complex filter output.
% This approach avoids computation of large FFTs as in Matlab's 'hilbert'.

if min(size(Fco)) > 1 | isstr(Fs), 
   b = Fco;	% kluge to specify filters
else b = quad_filt_bank(Fco, Fs);
end

e1_fs2 = zeros(size(x1));
e2_fs1 = zeros(size(x1));

% loop over filters in bank 
for k = 1:size(b,2),

	zfilt1 = fftfilt(b(:,k), x1);
	zfilt2 = fftfilt(b(:,k), x2);
	
	% interchange envelope and fine structure	
  	e1_fs2_filt = abs(zfilt1) .* cos(angle(zfilt2));
    e2_fs1_filt = abs(zfilt2) .* cos(angle(zfilt1));
      
    % refilter backward to avoid delay accumulation
    if refilter,
        len = length(e1_fs2_filt);
        e1_fs2_filt(len:-1:1,:) = fftfilt(real(b(:,k)), e1_fs2_filt(len:-1:1,:));
    	e2_fs1_filt(len:-1:1,:) = fftfilt(real(b(:,k)), e2_fs1_filt(len:-1:1,:));
	end
      
	% accumulate over frequency bands
  	e1_fs2 = e1_fs2 + e1_fs2_filt;
    e2_fs1 = e2_fs1 + e2_fs1_filt;
     
end
end

function b = quad_filt_bank(cutoffs, Fs, N)
% QUAD_FILT_BANK - Create bank of FIR complex filters
% 	whose real and imaginary parts are in quadrature
% Usage: b = quad_filt_bank(cutoffs, Fs, N)
%   cutoffs   band cutoff frequencies (Nbands + 1)
%   Fs        sampling rate (default 1)
%   N         filter order (default 8*Fs/(min(bandwidth)))
%   b         filter bank [N+1 X Nbands]
%
%	Copyright Bertrand Delgutte, 1999-2000
%
if nargin < 2, Fs = 1; end

w = 2*cutoffs/Fs;   % Matlab normalized frequency
if nargin < 3 | isempty(N), 
   N = round(8/min(diff(w))); 
	if rem(N,2) ~= 1, N = N+1; end	% make even
end

nbands = length(cutoffs) - 1;
b = zeros(N+1, nbands);

bw = diff(w)/2;    % lowpass filter bandwidths
fo = w(1:nbands) + bw;   % frequency offsets
t = [-N/2:N/2];     % time vector

for k = 1:nbands, 
   b(:,k) = (fir1(N, bw(k)) .* exp(j*pi*fo(k)*t)).'; 
end
end

