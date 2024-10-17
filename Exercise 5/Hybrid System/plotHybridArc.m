function plotHybridArc(t,j,x,jstar,resolution)
%PLOTHYBRIDARC   Hybrid arc plot (n states and n hybrid time domains).
%   PLOTHYBRIDARC(t,j,x) plots (in blue) the trajectory x on hybrid time
%   domains. The intervals [t_j,t_{j+1}] indexed by the corresponding j are
%   depicted in the t-j plane (in red).
%
%   PLOTHYBRIDARC(t,j,x,jstar) plots hybrid time vector (matrix) (t,j) versus
%   vector (matrix) x taking into account jumps j, and the plot is cut
%   regarding the jstar interval (jstar = [j-initial j-final]).
%
%   PLOTHYBRIDARC(t,j,x,[jstar],resolution) plots hybrid time vector (matrix)
%   (t,j) versus vector (matrix) x taking into account jumps j, and the
%   plot is cut regarding the jstar interval (jstar = [j-initial j-final]).
%   Also, a maximum resolution in between jumps is given by the input
%   variable resolution.

if ~exist('jstar','var') || isempty(jstar)
    jstar = [];
end
    modificatorF{1} = 'b';
    modificatorJ{1} = 'LineStyle';
    modificatorJ{2} = 'none';
if ~exist('resolution','var') || isempty(resolution)
    resolution = [];
end
DDD = true;

plotarc(t,j,x,[],jstar,modificatorF,modificatorJ,resolution,DDD);
hold on

modificatorF{1} = 'r';
modificatorJ{1} = 'LineStyle';
modificatorJ{2} = 'none';

plotarc(j,j,t,[],jstar,modificatorF,modificatorJ,resolution);

