function M = spblkdiag(varargin)
% SPBLKDIAG  Block diagonal concatenation of matrix input arguments,
%            with sparse output.
%
%
%                                         ( A 0 .. 0 )
% Y = SPBLKDIAG(A,B,...)  produces  SPARSE( 0 B .. 0 )
%                                         ( 0 0 ..   )
%
% but much more efficient than calling sparse after blkdiag
    % Basic checks on input/output consistency
    assert(nargin > 0,...
          [mfilename ':inargcount'],...
          'SPBLKDIAG requires at least 1 input matrix.');
    assert(nargout <= 1,...
          [mfilename  ':outargcount'],...
          'SPBLKDIAG produces at most 1 sparse output matrix.');
    assert(all(cellfun('isclass',varargin,'double')) && all(cellfun(@ndims,varargin)==2),...
          [mfilename  ':invalid_input'],...
          'SPBLKDIAG takes double-valued matrices as input arguments.');
    % MATLAB's own blkdiag() supports sparse output, but it is undocumented.
    % Therefore, force sparse output, and delegate to blkdiag():
    if nargin==1
        M = sparse(varargin{1});
    else
        if ~any(cellfun(@issparse, varargin))
            [~,ind] = min(cellfun('prodofsize', varargin));
            varargin{ind} = sparse(varargin{ind});
        end
        M = blkdiag(varargin{:});
    end
end
