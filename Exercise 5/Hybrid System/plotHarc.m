function [x_sliced,t_sliced] = plotHarc(t,j,x,jstar,modificatorF,modificatorJ,resolution)
%PLOTHARC   Hybrid arc plot (n states and n hybrid time domains).
%   [x_sliced,t_sliced] = PLOTHARC(t,j,x) plots hybrid time vector (matrix)
%   (t,j) versus vector (matrix) x taking into account jumps j. If x is a
%   matrix, then the vector is plotted versus the rows or columns of the
%   matrix, whichever line up. If t and j are a matrices, then each column
%   of x will be plotted according to the hybrid time domain composed for
%   each column of t and j. The function returns an array of cell elements
%   with x and t data indexed by j.
%
%   [x_sliced,t_sliced] = PLOTHARC(t,j,x,jstar) plots hybrid time vector
%   (matrix) (t,j) versus vector (matrix) x taking into account jumps j,
%   and the plot is cut regarding the jstar interval (jstar = [j-initial
%   j-final]).
%
%   PLOTHARC(t,j,x,[jstar],modificatorF,modificatorJ) plots hybrid time
%   vector (matrix) (t,j) versus vector (matrix) x taking into account
%   jumps j, and the plot is cut regarding the jstar interval (jstar =
%   [j-initial j-final]). The inputs modificatorF and modificatorJ modifies
%   the type of line used for flows and jumps, respectively. modificatorF
%   (modificatorJ) must be a cell array that contains the standard matlab
%   ploting modificators (see example). The default values are
%   modificatorF{1} = '', and modificatorJ{1} = '*--'.
%
%   PLOTHARC(t,j,x,[jstar],[modificatorF],[modificatorJ],resolution) plots
%   hybrid time vector (matrix) (t,j) versus vector (matrix) x taking into
%   account jumps j, and the plot is cut regarding the jstar interval
%   (jstar = [j-initial j-final]). Modificators must be cell arrays that
%   contains the standard matlab ploting modificators (see example). Also,
%   a maximum resolution in between jumps is given by the input variable
%   resolution.
%
%   Various line types, plot symbols and colors may be obtained with
%   PLOTHARC(t,j,x,jstar,modificator) where modificator is a cell array
%   created with the following strings:
%
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot
%          c     cyan          +     plus               --    dashed
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram

if ~exist('jstar','var') || isempty(jstar)
    jstar = [];
end
if ~exist('modificatorF','var') || isempty(modificatorF)
    modificatorF = [];
end
if ~exist('modificatorJ','var') || isempty(modificatorJ)
    modificatorJ = [];
end
if ~exist('resolution','var') || isempty(resolution)
    resolution = [];
end

[x_sliced,t_sliced] = plotarc(t,j,x,[],jstar,modificatorF,modificatorJ,resolution);
