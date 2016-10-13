function particles = q_mvnrnd(target)
%q_mvnrnd Creating M multivariate Gaussian random vectors from a
%D-dimensional quasi random sequence. 
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% target.mu     = mean
% target.chol   = Cholesky factor of the covariance matrix
% target.M      = desired number of samples
% target.Dims   = dimensionality
% target.qsi    = quasi-random point set
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% 
% particles     = multivariate Gaussian quasi-random vectors
%
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        05-Aug-2016
%   Last edited:    05-Aug-2016
%   Matlab version: 8.6.0.267246 (R2015b)
%	Email:          khoa.tran@unsw.edu.au
%
%   All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met: 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement: This product includes
%    software developed by the organization.
% 4. Neither the name of the organization nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY Khoa T. Tran ''AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL Khoa T. Tran BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% _________________________________________________________________________

%% Generate the D-dimensional quasi random set if none is found
if ~isfield(target, 'qsi') || isempty(target.qsi)
    target.qsi = haltonset(target.Dims,'Skip',randi(1e3),'Leap',randi(1e2));
end
%% Pick a random interval of length M
L       = length(target.qsi);
first   = randi(L); 
last    = first+target.M;
if last>=L; first = first - target.M; end
%% Sampling M points from the quasi random set
particles = target.qsi(first:last-1,:); 
particles = bsxfun(@minus,particles,min(particles));
particles = bsxfun(@rdivide,particles,max(particles));
%% Normalise the range of the set to [0,1]
idx = particles == 0 | particles == 1;
while any(idx(:))
particles(idx) = rand(nnz(idx),1);
idx = particles == 0 | particles == 1;
end
%% Take the Normal inverse cumulative distribution function
particles = norminv(particles,0,1); 
%% Add spatial correlation and shift the center to mu
particles = bsxfun(@plus,particles * target.chol,target.mu);