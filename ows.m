function resignal = ows(frames, hop_size, fadewin)
% Reconstruction of a signal from its time-domain frames
% ------------------------------------------------------
%
% Syntax: resignal = ows(frames, hop_size, recon)
%
%   INPUTS:   frames   - a matrix (frame_length x no. frames) that contains the signal frames
%             hop_size - the hop size in samples used during the analysis
%             fadewin  - reconstruction window
%
%   OUTPUTS:  resignal - the reconstructed signal
%
% ------------------------------------------------------
%
% Author: Elias Kokkinis
% Version: 5.1
% Last revision: 04/04/2015
%
% ======================================================
frame_length = size(frames, 1);

if nargin < 3
    fadewin = hann(frame_length)./hamming(frame_length);
end

% Length of output signal
L = size(frames, 2)*hop_size + frame_length;
% Pre-allocate for speed
resignal = zeros(L,1);

for i = 1:size(frames, 2)
 resignal(1 + ((i - 1)*hop_size):frame_length + ((i - 1)*hop_size))=resignal(1 + ((i - 1)*hop_size):frame_length + ((i - 1)*hop_size)) + fadewin.*frames(:, i);   
end
