function [msh Q A] = compQ_v3(msh, trancsts, nfrac, vfield)
%    Copyright 2014 Bart Vermeulen
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.
%
%    input
%    msh       : adcpttools mesh
%    trancsts : trancsts to process (default: all)
%
%
% TODO zero value boundary conditions should be given at bottom and
%      side walls in combination with linear interpolation,
%      while interpolating the top still with the nearest neighbour method
%      better: quadratic interpolation with enforced bank angle of say 30 deg at surface
% TODO linear extrapolation of bank position may yield implausible large values
%      if the bottom slope is low or if the GPS coordinates contain an outlier
% automatically select all trancsts
if (nargin() < 2 || isempty(trancsts))
    trancsts = 1:length(msh);
end

% relative size of cross section that is used to determination the bank position
%nfrac=[0.1 0.3 0.3 0.5 0.3 0.1 0.1] ;
if (nargin()<3 || isempty(nfrac))
    nfrac = 0.2*ones(size(trancsts));
end

if (nargin < 4)
    vfield = 'pars'; %estimated velocity field
end
progvelfield = ['prog' vfield];
sfield = 'cs';


end