%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Region of Attraction Boundary Tracing (RABT) algorithm
%   - uses COSY numerical continuation tool to approximate the boundary for
%   the region of attraction of a stable equilibrium
%
% Inputs:
%   F       - function as handle
%   equil   - stable equilibrium value
%   plimit  - parameter continuation limits, last entry is step size
%   steps   - number of RK steps for each continuation step
%   ztol    - COSY zero tolerance
% 
% Author: Peter Uth
% Created: July 2017
%--------------------------------------------------------------------------
function basin = rabt(F,equil,plimit,steps,ztol)
%% Run test simulations to find initial boundary points

% move in positive direction of x2, starting from equilibrium
Fsim = str2func(strrep(func2str(F),'@(x,u,p)','@(t,x)'));
x0 = [equil(1) equil(2)];
xb0 = [];
simstep = 0.1;
bound = 0;
while x0 <= plimit(2,2)
    x0 = [x0;x0(end,1) x0(end,2)+simstep];
    tf = 20;
    [~,x] = ode15s(Fsim,[0 tf],x0(end,:));
    if abs(x(end,1)-equil(1)) > 0.1 || abs(x(end,2)-equil(2)) > 0.1
        x0 = [x0;x0(end,1) x0(end-1,2)];
        simstep = simstep/10;
        bound = 1;
    end
    
    if simstep < 0.0001
        break
    end  
    
    disp(simstep)
    disp(x0(end,:))
end
if bound == 1
    xb0 = x0(end,:);
end

% move in negative direction of x2, starting from equilibrium
x0 = [equil(1) equil(2)];
simstep = 0.1;
bound = 0;
while x0 >= plimit(2,1)
    x0 = [x0;x0(end,1) x0(end,2)-simstep];
    tf = 20;
    [~,x] = ode15s(Fsim,[0 tf],x0(end,:));
    if abs(x(end,1)-equil(1)) > 0.1 || abs(x(end,2)-equil(2)) > 0.1
        x0 = [x0;x0(end,1) x0(end-1,2)];
        simstep = simstep/10;
        bound = 1;
    end
    
    if simstep < 0.0001
        break
    end
    
    disp(simstep)
    disp(x0(end,:))
end

if bound == 1
    xb0 = [xb0;x0(end,:)];
end
        
%% Apply COSY numerical continuation to each boundary point
% Additional parameter: integrator step size
global fun step type
fun = F;
step = steps;
type = 'stepsize';

nx = size(plimit,1)-1;      % # of states
np = nx+1;                  % # of parameters = nx + 1 (step size h)
ng = nx;                    % # of constraints = # of states

X = [];
for i = 1:size(xb0,1)
    clear b x0 branches p0 j
    
    b = cosy(@rabt_iter,nx,0,np,ng,0);

    b.p_min = plimit(:,1); 
    b.p_max = plimit(:,2);

    b.cm.solution_zero_tol = ztol;
    b.cm.max_steps = 10000;
    b.cm.ds_min = 0.001;

    b.ig_zero = 1:ng;
    b.iy_zero = [];
    b.iu_free = [];
    u0 = 0;

    x0 = equil;
    p0 = xb0(i,:);
    p0 = [p0 plimit(end,1)];
    p0 = p0'; 

    [branches,~] = b.trace_solution_branches(x0,u0,p0);  
    
    for j = 1:size(branches,2)
        X = [X branches(j).p];
    end
    assignin('base','branches',branches);
end

% basin = X;
if ~isempty(X)
    X_idx = boundary(X(1,:)',X(2,:)',0.99);
    basin = X(:,X_idx);
else
    basin = X;
end

%% Check for incompleteness
% if incomplete, rerun COSY with epsilon tolerance as additional parameter
if isempty(basin) || pdist([basin(1,2) basin(2,2);basin(1,end) basin(2,end)],'euclidean') > 0.1 
    type = 'tolerance';
    X = [];
    for i = 1:size(xb0,1)
        clear b x0 branches p0 j

        b = cosy(@rabt_iter,nx,0,np,ng,0);

        b.p_min = plimit(:,1);
        b.p_min(np) = -0.5;
        b.p_max = plimit(:,2);
        b.p_max(np) = 0.5;

        b.cm.solution_zero_tol = ztol;
        b.cm.max_steps = 10000;
        b.cm.ds_min = 0.001;

        b.ig_zero = 1:ng;
        b.iy_zero = [];
        b.iu_free = [];
        u0 = 0;

        x0 = equil;
        p0 = xb0(i,:);
        p0 = [p0 0];
        p0 = p0'; 

        [branches,~] = b.trace_solution_branches(x0,u0,p0);  

        for j = 1:size(branches,2)
            X = [X branches(j).p];
        end
        assignin('base','branches',branches);
    end
    
    basin = [basin X];

end
if ~isempty(X)
    X_idx = boundary(basin(1,:)',basin(2,:)',0.99);
    basin = basin(:,X_idx);
else
    basin = X;
end
   