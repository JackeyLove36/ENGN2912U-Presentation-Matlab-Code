function [x_sliced,t_sliced] = plotHarcColor(t,j,x,L,jstar,resolution)
%PLOTHARCCOLOR   Hybrid arc plot with color (n states and n hybrid time
%domains).
%   [x_sliced,t_sliced] = PLOTHARCCOLOR(t,j,x,L) plots hybrid time vector
%   (matrix) (t,j) versus vector (matrix) x taking into account jumps j.
%   The hybrid arc is plotted with L data (matrix) as color. The input
%   vectors (matrices) t, j, x, L must have the same length (number of
%   rows).
%
%   [x_sliced,t_sliced] = PLOTHARCCOLOR(t,j,x,L,jstar) plots hybrid time
%   vector (matrix) (t,j) versus vector (matrix) x taking into account
%   jumps j. The hybrid arc is plotted with L data (matrix) as color, and
%   the plot is cut regarding the jstar interval (jstar = [j-initial
%   j-final]). The parameter L is NOT optional.
%
%   [x_sliced,t_sliced] = PLOTHARCCOLOR(t,j,x,L,[jstar],resolution) plots
%   hybrid time vector (matrix) (t,j) versus vector (matrix) x taking into
%   account jumps j, and the plot is cut regarding the jstar interval
%   (jstar = [j-initial j-final]). jstar is optional. Also, a maximum
%   resolution in between jumps is given by the input variable resolution.

if ~exist('L','var') || isempty(L)
    L = [];
end
if ~exist('jstar','var') || isempty(jstar)
    jstar = [];
end
if ~exist('resolution','var') || isempty(resolution)
    resolution = [];
end

[x_sliced,t_sliced] = plotarc(t,j,x,L,jstar,[],[],resolution);
