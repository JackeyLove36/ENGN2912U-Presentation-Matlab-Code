function [x_sliced,t_sliced] = plotarc(t,j,x,L,jstar,modificatorF,modificatorJ,resolution,DDD,true3D)
%   PLOTARC   Hybrid arc plot (n states and n hybrid time domains).
%   [x_sliced,t_sliced] = PLOTARC(t,j,x) plots hybrid time vector (matrix)
%   (t,j) versus vector (matrix) x taking into account jumps j. If x is a
%   matrix, then the vector is plotted versus the rows or columns of the
%   matrix, whichever line up. If t and j are a matrices, then each column
%   of x will be plotted according to the hybrid time domain composed for
%   each column of t and j. The function returns an array of cell elements
%   with x and t data indexed by j.
%
%   [x_sliced,t_sliced] = PLOTARC(t,j,x,L) plots hybrid time vector
%   (matrix) (t,j) versus vector (matrix) x taking into account jumps j.
%   The hybrid arc is plotted with L data (matrix) as color. The input
%   vectors (matrices) t, j, x, L must have the same length (number of
%   rows).
%
%   [x_sliced,t_sliced] = PLOTARC(t,j,x,[L],jstar) plots hybrid time vector
%   (matrix) (t,j) versus vector (matrix) x taking into account jumps j,
%   and the plot is cut regarding the jstar interval (jstar = [j-initial
%   j-final]). The parameter L is optional.
%
%   [x_sliced,t_sliced] = PLOTARC(t,j,x,[L],[jstar],modificatorF,modificatorJ)
%   plots hybrid time vector (matrix) (t,j) versus vector (matrix) x taking
%   into account jumps j, and the plot is cut regarding the jstar interval
%   (jstar = [j-initial j-final]). The inputs modificatorF and modificatorJ
%   modifies the type of line used for flows and jumps, respectively.
%   modificatorF (modificatorJ) must be a cell array that contains the
%   standard matlab ploting modificators (see example). The default values
%   are modificatorF{1} = '', and modificatorJ{1} = '*--'.
%
%   Various line types, plot symbols and colors may be obtained with
%   PLOTARC(t,j,x,[],[jstar],modificatorF,modificatorJ) where modificator
%   is a cell array created with the following strings:
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
%
%   [x_sliced,t_sliced] = PLOTARC(t,j,x,[L],[jstar],[modificatorF],
%   [modificatorJ],resolution) plots hybrid time vector (matrix) (t,j)
%   versus vector (matrix) x taking into account jumps j, and the plot is
%   cut regarding the jstar interval (jstar = [j-initial j-final]).
%   Modificators must be cell arrays that contains the standard matlab
%   ploting modificators (see example). Also, a maximum resolution in
%   between jumps is given by the input variable resolution.

%% check matlab version
mver = ver('matlab');
mverok = ~verLessThan('matlab', '8.4');
if mverok==0 % pre R2014b plot behaviour
    ColOrd = get(0,'DefaultAxesColorOrder');
end


%  if possible, fast plot
fastpl = false;


%% input management
if ~exist('L','var')
    L = [];
end
% Define type of plot
if isempty(L) % hybrid arc plot
    color = false;
else
    color = true;
end
if ~exist('DDD','var') || isempty(DDD)
    DDD = false;
end
if ~exist('true3D','var') || isempty(true3D)
    true3D = false;
end
    
[rx,cx] = size(x);
[rt,ct] = size(t);
[rj,cj] = size(j);
[rl,cl] = size(L);
% reorientate x, t, or j if necessary
if rt~=rj || ct~=cj || (color==true && (rl~=rj||cl~=cj))
    fprintf('Error, t, j, or L does not have the proper sizes\n')
    return
elseif rx==rt
    nt = rt;
elseif rx==ct && rx~=cx
    nt = rx;
    t = t';
    j = j';
    L = L';
    [rt,ct] = size(t);
    [rj,cj] = size(j);
    [rl,cl] = size(L);
elseif ct==1 && cx==rt
    nt = rt;
    x = x';
    [rx,cx] = size(x);
else
    fprintf('Error, x does not have the proper size\n')
    return
end

if ct>cx
    fprintf('Error, there are more hybrid time domains than arcs\n')
    return
elseif rt == 1
    nt = ct;
    x = x';
    t = t';
    j = j';
    L = L';
    [rx,cx] = size(x);
    [rt,ct] = size(t);
    [rj,cj] = size(j);
    [rl,cl] = size(L);
elseif ct<cx
    if ct==1
        fastpl = true; % fast plot
    end
    ttemp = repmat(t(:,1),1,cx-ct);
    t = [t,ttemp];
    jtemp = repmat(j(:,1),1,cx-ct);
    j = [j,jtemp];
    if color
        Ltemp = repmat(L(:,1),1,cx-ct);
        L = [L,Ltemp];
        [rl,cl] = size(L);
    end
    [rt,ct] = size(t);
    [rj,cj] = size(j);
    % fprintf('Warning, there are some missing hybrid time domains\n')
end

if true3D && cx == 3
    true3D = true;
    DDD = false;
else
    true3D = false;
end

if ~exist('jstar','var') || isempty(jstar) % pick the max and the min!
    jstar = [min(j(1,:)) max(j(end,:))];
    iini = jstar(1);
    ifini = jstar(end);
elseif exist('jstar','var')
    iini = jstar(1);
    try
        ifini = jstar(2);
    catch
        ifini = jstar(1);
    end
    
end

if ~exist('modificatorF','var') || isempty(modificatorF)
    modificatorF{1} = '';
end
if DDD || true3D
    if ~exist('modificatorJ','var') || isempty(modificatorJ)
        modificatorJ{1} = '--';
    end
else
    if ~exist('modificatorJ','var') || isempty(modificatorJ)
        modificatorJ{1} = '*--';
    end
end

if ~exist('resolution','var') || isempty(resolution)
    resolution = nt;
end

%% output management

% Outputs
nout = nargout;

%  if possible, fast plot
if color==false && cj>1
    for ij = 1:rj
        fastpl = sum(j(ij,:)==j(ij,1))==cj;
        if fastpl==0
            break;
        end
    end
elseif cj==1
    fastpl = 1;
end

% ------------- Fast plot
if fastpl==1 && color == false % Fast plot possible
    if nout~=0 % improve speed
        x_sliced{ifini-iini+1,cx} = [];
        t_sliced{ifini-iini+1,cx} = [];
    end
    for ij=1:ifini-iini+1
        indexi = find(j(:,1) == iini+ij-1);
        iij = indexi(1);
        if (length(indexi) > resolution)
            step = round(length(indexi) / resolution);
            indexitemp = indexi(1:step:end);
            if indexitemp(end) ~= indexi(end)
                indexitemp(end+1) = indexi(end);
            end
            indexi = indexitemp;
        end
        if nout~=0 % improve speed
            for ix=1:cx
                x_sliced{ij,ix} = x(indexi,ix);
                t_sliced{ij,ix} = t(indexi,ix);
            end
        end
        if mverok==1 % post R2014b plot behaviour
            set(gca,'ColorOrderIndex',1);
        end
        if length(indexi)>1
            if DDD
                plot3(j(indexi,:),t(indexi,:),x(indexi,1),modificatorF{1:end})
            elseif true3D
                plot3(x(indexi,1),x(indexi,2),x(indexi,3),modificatorF{1:end})                
            else
                plot(t(indexi,:),x(indexi,:),modificatorF{1:end})
            end
            hold on;
        end
        if ij>1
            if mverok==1 && color==false % post R2014b plot behaviour
                set(gca,'ColorOrderIndex',1);
            end
            if DDD
                plot3(j(iij-1:iij,:),t(iij-1:iij,:),x(iij-1:iij,:),...
                    modificatorJ{1:end});
            elseif true3D
                plot3(x(iij-1:iij,1),x(iij-1:iij,2),x(iij-1:iij,3),...
                    modificatorJ{1:end});                
            else
                plot(t(iij-1:iij,:),x(iij-1:iij,:),...
                    modificatorJ{1:end});
            end
        end
    end
    hold off;
    return
    % ------------- arc plot and color plot
elseif color || fastpl==0 % color plot or no fast plot
    if true3D
        eacharc = 1;
    else
        eacharc = cx;
    end
    for ix=1:eacharc % for each arc
        for ij=1:ifini-iini+1
            indexi = find(j(:,ix) == iini+ij-1);
            if ~isempty(indexi)
                iij = indexi(1);
                if (length(indexi) > resolution)
                    step = round(length(indexi) / resolution);
                    indexitemp = indexi(1:step:end);
                    if indexitemp(end) ~= indexi(end)
                        indexitemp(end+1) = indexi(end);
                    end
                    indexi = indexitemp;
                end
                if nout~=0 % improve speed
                    for tempix=1:cx
                        x_sliced{ij,tempix} = x(indexi,tempix);
                        t_sliced{ij,tempix} = t(indexi,tempix);
                    end
                end
                if color % COLOR PLOT
                    col{ij} = L(indexi,ix);
                    if length(indexi)>1
                        if DDD
                            Xdata = [j(indexi,ix),j(indexi,ix)];
                            Ydata = [t(indexi,ix),t(indexi,ix)];
                            Zdata = [x(indexi,ix),x(indexi,ix)];
                            Cdata = [col{ij},col{ij}];
                        elseif true3D
                            Xdata = [x(indexi,1),x(indexi,1)];
                            Ydata = [x(indexi,2),x(indexi,2)];
                            Zdata = [x(indexi,3),x(indexi,3)];
                            Cdata = [col{ij},col{ij}];                            
                        else
                            Xdata = [t(indexi,ix),t(indexi,ix)];
                            Ydata = [x(indexi,ix),x(indexi,ix)];
                            Zdata = [zeros(size(indexi)),zeros(size(indexi))];
                            Cdata = [col{ij},col{ij}];
                        end
                    else
                        if DDD
                            Xdata = [j(indexi,ix)*ones(2)];
                            Ydata = [t(indexi,ix)*ones(2)];
                            Zdata = [x(indexi,ix)*ones(2)];
                            Cdata = [col{ij}*ones(2)];
                        elseif true3D                            
                            Xdata = [x(indexi,1)*ones(2)];
                            Ydata = [x(indexi,2)*ones(2)];
                            Zdata = [x(indexi,3)*ones(2)];
                            Cdata = [col{ij}*ones(2)];
                        else
                            Xdata = [t(indexi,ix)*ones(2)];
                            Ydata = [x(indexi,ix)*ones(2)];
                            Zdata = [zeros(size(indexi))*ones(2)];
                            Cdata = [col{ij}*ones(2)];
                        end                        
                    end
                    surface(...
                        'XData',Xdata,...
                        'YData',Ydata,...
                        'ZData',Zdata,...
                        'CData',Cdata,...
                        'facecol','no','edgecol','flat','linew',2);                    
                    hold on;
                    try
                        if DDD
                            Xdata = [j(iij-1:iij,ix),j(iij-1:iij,ix)];
                            Ydata = [t(iij-1:iij,ix),t(iij-1:iij,ix)];
                            Zdata = [x(iij-1:iij,ix),x(iij-1:iij,ix)];
                            Cdata = [L(iij-1:iij,ix),L(iij-1:iij,ix)];
                            surface(...
                                'XData',Xdata,...
                                'YData',Ydata,...
                                'ZData',Zdata,...
                                'CData',Cdata,...
                                'LineStyle','--',...
                                'facecol','no','edgecol','flat','linew',1);
                        elseif true3D
                            Xdata = [x(iij-1:iij,1),x(iij-1:iij,1)];
                            Ydata = [x(iij-1:iij,2),x(iij-1:iij,2)];
                            Zdata = [x(iij-1:iij,3),x(iij-1:iij,3)];
                            Cdata = [L(iij-1:iij,ix),L(iij-1:iij,ix)];
                            surface(...
                                'XData',Xdata,...
                                'YData',Ydata,...
                                'ZData',Zdata,...
                                'CData',Cdata,...
                                'LineStyle','--',...
                                'facecol','no','edgecol','flat','linew',1);                            
                        else
                            Xdata = [t(iij-1:iij,ix),t(iij-1:iij,ix)];
                            Ydata = [x(iij-1:iij,ix),x(iij-1:iij,ix)];
                            Zdata = [zeros(size(t(iij-1:iij,ix))),zeros(size(t(iij-1:iij,ix)))];
                            Cdata = [L(iij-1:iij,ix),L(iij-1:iij,ix)];
                            surface(...
                                'XData',Xdata,...
                                'YData',Ydata,...
                                'ZData',Zdata,...
                                'CData',Cdata,...
                                'LineStyle','--','Marker','*',...
                                'facecol','no','edgecol','flat','linew',1);
                        end
                    catch
                    end
                elseif fastpl==0 % NO FAST PLOT
                    if  mverok==0 % pre R2014b plot behaviour
                        set(0,'DefaultAxesColorOrder',circshift(ColOrd,1-ix))
                        set(gca,'ColorOrder',circshift(ColOrd,1-ix))
                    end
                    if length(indexi)>1
                        if mverok==1 % post R2014b plot behaviour
                            set(gca,'ColorOrderIndex',ix);
                        end
                        if DDD
                            plot3(j(indexi,ix),t(indexi,ix),x(indexi,ix),...
                                modificatorF{1:end})
                        elseif true3D
                            plot3(x(indexi,1),x(indexi,2),x(indexi,3),...
                                modificatorF{1:end})                            
                        else
                            plot(t(indexi,ix),x(indexi,ix),...
                                modificatorF{1:end})
                        end
                        hold on;
                    end
                    try
                        if mverok==1 % post R2014b plot behaviour
                            set(gca,'ColorOrderIndex',ix);
                        end
                        if DDD
                            plot3(j(iij-1:iij,ix),t(iij-1:iij,ix),x(iij-1:iij,ix),...
                                modificatorJ{1:end});
                        elseif true3D
                            plot3(x(iij-1:iij,1),x(iij-1:iij,2),x(iij-1:iij,3),...
                                modificatorJ{1:end});                            
                        else
                            plot(t(iij-1:iij,ix),x(iij-1:iij,ix),...
                                modificatorJ{1:end});
                        end
                    catch
                    end
                    
                end
            end
        end
    end
    if  mverok==0% pre R2014b plot behaviour
        set(0,'DefaultAxesColorOrder',ColOrd)
    end
end
