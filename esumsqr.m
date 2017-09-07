function SSQ=esumsqr(varargin)
% Element-wise sum of squares for arrays
% SSQ = esumsqr(A1,A2,...)
%

nVarargs=nargin;    % number of variable arg inputs (all)

% Input error checking:
%   - input arguments must be N-dim arrays of same size
ii=1;
asize=size(varargin{ii});
while ii<nVarargs
    ii=ii+1;
    if ~isequal(size(varargin{ii}),asize)
        error('Input arrays must be of constant size.');
    end
end

% do the elementary sum of squares!
SSQ=0;
for ii=1:nVarargs
    SSQ=SSQ+varargin{ii}.^2;
end