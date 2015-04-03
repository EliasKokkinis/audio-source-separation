function frames = owa(s, frame_length, hop_size, wnd)
% Analysis of a signal to (windowed) time-domain frames
% ------------------------------------------------------
%
% Syntax: frames = owa(s, frame_length, hop_size)
%
%   INPUTS:   s            - a vector containing the time-domain signal samples
%             frame_length - the window length in samples
%             hop_size     - the hop size in samples
%             wnd          - a vector (of length frame_length) that contains the window samples
%
%   OUTPUTS:  frames       - a matrix (no.frames x frame_length) that contains the (windowed) signal frames.
%
% ------------------------------------------------------
%
% Author: Elias Kokkinis
% Version: 3.2
% Last revision: 04/04/2015
%
% ======================================================
if nargin < 4
    wnd = hamming(frame_length, 'periodic');
end

% Ensure column vectors
s = s(:);
wnd = wnd(:);
% Zero pad
spad = [s;zeros(frame_length,1)]; 

% Number of frames
N = floor((length(spad) - frame_length)/hop_size) + 1; 

% Preallocate for speed
frames = zeros(frame_length, N);

fno = 1;
for idx = 1:hop_size:(length(spad) - frame_length)
    frames(:, fno) = wnd.*spad(idx:idx + frame_length - 1);
    fno = fno + 1;
end