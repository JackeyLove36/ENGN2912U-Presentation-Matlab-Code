function [t,j,x] = HyEQsolver(f,g,C,D,x0,TSPAN,JSPAN,rule,options,solver,E)
%HYEQSOLVER solves hybrid equations.
%   Syntax: [t j x] = HYEQSOLVER(f,g,C,D,x0,TSPAN,JSPAN,rule,options,solver,E)
%   computes solutions to the hybrid equations
%
%   \dot{x} = f(x,t,j)  x \in C x^+ = g(x,t,j)  x \in D
%
%   where x is the state, f is the flow map, g is the jump map, C is the
%   flow set, and D is the jump set. It outputs the state trajectory (t,j)
%   -> x(t,j), where t is the flow time parameter and j is the jump
%   parameter.
%
%   x0 defines the initial condition for the state.
%
%   TSPAN = [TSTART TFINAL] is the time interval. JSPAN = [JSTART JSTOP] is
%       the interval for discrete jumps. The algorithm stop when the first
%       stop condition is reached.
%
%   rule (optional parameter) - rule for jumps
%       rule = 1 (default) -> priority for jumps rule = 2 -> priority for
%       flows
%
%   options (optional parameter) - options for the solver see odeset f.ex.
%       options = odeset('RelTol',1e-6);
%       options = odeset('InitialStep',eps);
%
%   solver (optional parameter. String) - selection of the desired ode
%       solver. All ode solvers are suported, exept for ode15i.  See help
%       odeset for detailed information.
%
%   E (optional parameter) - Mass matrix [constant matrix | function_handle]
%       For problems: 
%       E*\dot{x} = f(x) x \in C 
%       x^+ = g(x)  x \in D
%       set this property to the value of the constant mass matrix. For
%       problems with time- or state-dependent mass matrices, set this
%       property to a function that evaluates the mass matrix. See help
%       odeset for detailed information.

if ~exist('rule','var')
    rule = 1;
end

if ~exist('options','var')
    options = odeset();
end
if exist('E','var') && ~exist('solver','var')
    solver = 'ode15s';
end
if ~exist('solver','var')
    solver = 'ode45';
end
if ~exist('E','var')
    E = [];
end
% mass matrix (if existent)
isDAE = false;
if ~isempty(E)
    isDAE = true;
    switch isa(E,'function_handle')
        case true % Function E(x)
            M = E;
            options = odeset(options,'Mass',M,'Stats','off',...
                'MassSingular','maybe','MStateDependence','strong',...
                'InitialSlope',f_hdae(x0,TSPAN(1))); 
        case false % Constant double matrix
            M = double(E);
            options = odeset(options,'Mass',M,'Stats','off',...
                'MassSingular','maybe','MStateDependence','none');
    end
end

odeX = str2func(solver);
nargf = nargin(f);
nargg = nargin(g);
nargC = nargin(C);
nargD = nargin(D);

% simulation horizon
tstart = TSPAN(1);
tfinal = TSPAN(end);
jout = JSPAN(1);
j = jout(end);

% simulate
tout = tstart;
[rx,cx] = size(x0);
if rx == 1
    xout = x0;
elseif cx == 1
    xout = x0.';
else
    error('Error, x0 does not have the proper size')
end

% Jump if jump is prioritized:
if rule == 1
    while (j<JSPAN(end))
        % Check if value it is possible to jump current position
        insideD = fun_wrap(xout(end,:).',tout(end),j,D,nargD);
        if insideD == 1
            [j, tout, jout, xout] = jump(g,j,tout,jout,xout,nargg);
        else
            break;
        end
    end
end
fprintf('Completed: %3.0f%%',0);
while (j < JSPAN(end) && tout(end) < TSPAN(end))
    options = odeset(options,'Events',@(t,x) zeroevents(x,t,j,C,D,...
        rule,nargC,nargD));
    % Check if it is possible to flow from current position
    insideC = fun_wrap(xout(end,:).',tout(end),j,C,nargC);
    if insideC == 1
        if isDAE
            options = odeset(options,'InitialSlope',f(xout(end,:).',tout(end)));
        end
        [t,x] = odeX(@(t,x) fun_wrap(x,t,j,f,nargf),[tout(end) tfinal],...
            xout(end,:).', options);
        nt = length(t);
        tout = [tout; t];
        xout = [xout; x];
        jout = [jout; j*ones(1,nt)'];
    end
    
    %Check if it is possible to jump
    insideD = fun_wrap(xout(end,:).',tout(end),j,D,nargD);
    if insideD == 0
        break;
    else
        if rule == 1
            while (j<JSPAN(end))
                % Check if it is possible to jump from current position
                insideD = fun_wrap(xout(end,:).',tout(end),j,D,nargD);
                if insideD == 1
                    [j, tout, jout, xout] = jump(g,j,tout,jout,xout,nargg);
                else
                    break;
                end
            end
        else
            [j, tout, jout, xout] = jump(g,j,tout,jout,xout,nargg);
        end
    end
    fprintf('\b\b\b\b%3.0f%%',max(100*j/JSPAN(end),100*tout(end)/TSPAN(end)));
end
t = tout;
x = xout;
j = jout;
fprintf('\nDone\n');
end

function [value,isterminal,direction] = zeroevents(x,t,j,C,D,rule,nargC,nargD)
switch rule
    case 1 % -> priority for jumps
        isterminal(1) = 1; % InsideC
        isterminal(2) = 1; % Inside(C \cap D)
        isterminal(3) = 1; % OutsideC
        direction(1) = -1; % InsideC
        direction(2) = -1; % Inside(C \cap D)
        direction(3) =  1; % OutsideC
    case 2 %(default) -> priority for flows
        isterminal(1) = 1; % InsideC
        isterminal(2) = 0; % Inside(C \cap D)
        isterminal(3) = 1; % OutsideC
        direction(1) = -1; % InsideC
        direction(2) = -1; % Inside(C \cap D)
        direction(3) =  1; % OutsideC
end

insideC = fun_wrap(x,t,j,C,nargC);
insideD = fun_wrap(x,t,j,D,nargD);
outsideC = -fun_wrap(x,t,j,C,nargC);


value(1) = 2*insideC;
value(2) = 2-insideC - insideD;
value(3) = 2*outsideC;

end


function [j, tout, jout, xout] = jump(g,j,tout,jout,xout,nargfun)
% Jump
j = j+1;
y = fun_wrap(xout(end,:).',tout(end),jout(end),g,nargfun); 
% Save results
tout = [tout; tout(end)];
xout = [xout; y.'];
jout = [jout; j];
end

function xdelta = fun_wrap(x,t,j,h,nargfun)
%fun_wrap   Variable input arguments function (easy use for users).
%   fun_wrap(x,t,j,h,nargfun) depending on the function h written by the
%   user, this script selects how the HyEQ solver should call that
%   function.
%    x: state
%    t: time
%    j: discrete time
%    h: function handle
%    nargfun: number of input arguments of function h    

switch nargfun
    case 1
        xdelta = h(x);
    case 2
        xdelta = h(x,t);
    case 3
        xdelta = h(x,t,j);        
end
end