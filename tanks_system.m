clear all, close all 

%% 5 tanks system params
g   = 9.8;                 % Gravity acceleration
at  = 2.5*2.5*pi;        % Cross sectional area of cylinders
Aij = 0.15*0.15*pi;      % Orifice area between tanks 
Aii = 0.05*0.05*pi;      % Outlet orifice area 
Cii = 0.5;               % Coefficient of discharge (outlet)
Cij = 0.75;              % Coefficient of discharge (between tanks)

kii = Cii*Aii;
kij = Cij*Aij;

% Set-point levels between 5m and 10m
hmin = 5; 
hmax = 10;

ai_min = kii*sqrt(2*g)/2/at/sqrt(hmax);
ai_max = kii*sqrt(2*g)/2/at/sqrt(hmin);

aij_max = kij*sqrt(2*g)/2/at/sqrt(abs(0.1));
aij_min = 0;


% Tank 1 => p1 : a1, p2 : a12
% Tank 2 => p3 : a2, p4 : a23, p5 : a25
% Tank 3 => p6 : a3, p7 : a34
% Tank 4 => p8 : a4
% Tank 5 => p9 : a5 

Ap = @(p)[-p(1)-p(2) p(2)                   0               0           0;
           p(2)      -p(3)-p(2)-p(4)-p(5)   p(4)            0           p(5);
           0         p(4)                  -p(6)-p(4)-p(7)  p(7)        0;
           0         0                      p(7)           -p(8)-p(7)   0;
           0         p(5)                   0               0          -p(5)-p(9)];

B = eye(5,5);

m = 5;
n = 5;
             %p1     %p2       %p3     %p4      %p5       %p6     %p7       %p8       %p9                                                          
pBounds = [ ai_min   aij_min   ai_min  aij_min  aij_min   ai_min  aij_min   ai_min   ai_min      
            ai_max   aij_max   ai_max  aij_max  aij_max   ai_max  aij_max   ai_max   ai_max];


% Create vertices set
vv = {};
for i=1:length(pBounds)
    vv{end+1} = pBounds(:,i);
end
Dv = vv;
[Dv{:}] = ndgrid(vv{:});
paramCombs = cell2mat(cellfun(@(m)m(:),Dv,'uni',0));

%% Variables for the LMI design 
K_struct = [1 1 0 0 0
            1 1 1 0 1
            0 1 1 1 0
            0 0 1 1 0
            0 1 0 0 1];      
% F_0  has the same strcuture of K
F_0 = sdpvar(m,n,'full').*K_struct;

K_struct_1 = [1 1 0 0 0
              1 1 1 0 1
              0 0 0 0 0
              0 0 0 0 0
              0 0 0 0 0];
F_1 = sdpvar(m,n,'full').*K_struct_1;

K_struct_2 = K_struct_1;
F_2 = sdpvar(m,n,'full').*K_struct_2;

K_struct_3 = [1 1 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 0 0 0 0
              0 1 0 0 1];
F_3 = sdpvar(m,n,'full').*K_struct_3;

K_struct_4 = [0 0 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 0 0 0 0
              0 0 0 0 0];
F_4 = sdpvar(m,n,'full').*K_struct_4;

K_struct_5 = [0 0 0 0 0
              1 1 1 0 1
              0 0 0 0 0
              0 0 0 0 0
              0 1 0 0 1];
F_5 = sdpvar(m,n,'full').*K_struct_5;

K_struct_6 = [0 0 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 0 1 1 0
              0 0 0 0 0];
F_6 = sdpvar(m,n,'full').*K_struct_6;

K_struct_7 = [0 0 0 0 0
              0 0 0 0 0
              0 1 1 1 0
              0 0 1 1 0
              0 0 0 0 0];
F_7 = sdpvar(m,n,'full').*K_struct_7;

K_struct_8 = K_struct_7;
F_8 = sdpvar(m,n,'full').*K_struct_8;

K_struct_9 = K_struct_5;
F_9 = sdpvar(m,n,'full').*K_struct_9;

% M (the auxilary matrix) has the structure defined by X_struct
% Define the structure set, X, and Lambda 
L = basisSet(K_struct); % Define L
X = sym('q',n);
Lambda = sym('l',nnz(K_struct));
% Define equation 
eq = L*kron(eye(nnz(K_struct)),X)-L*kron(Lambda,eye(n));
% Find structure of X and Lambda
[X_struct, Lambda_struct] = findStruct(X,Lambda,eq);
Xs = zeros(n,n);
Xs(find(X_struct))=1;

% Xs = [1     1     0     0     0
%       0     1     0     0     0
%       0     0     1     0     0
%       0     0     1     1     0
%       0     1     0     0     1];
X = sdpvar(n,n,'full').*Xs;

% Lyapunov function matrix P
Q = sdpvar(n,n,'symmetric');

%% LMI Problem
eps = 1e-7; % Tolerance value
LMIoptions = sdpsettings('solver','Sedumi','Verbose',1); 

LMI = [];
LMI = [LMI; Q >= eps*eye(n)]; % Positive definiteness of P

alpha = 0.5;

for i = 1:length(paramCombs)
    % Vertex evaluation
    v_i    = paramCombs(i,:);  % Current vertex
    A_i    = Ap(v_i);          % Current A    
    
    BF_i   = B*F_0 + ...
             B*F_1*v_i(1)+ B*F_2*v_i(2)+ B*F_3*v_i(3) + B*F_4*v_i(4) + B*F_5*v_i(5) + ...
             B*F_6*v_i(6)+ B*F_7*v_i(7)+ B*F_8*v_i(8) + B*F_9*v_i(9);

     % Structured feedback LMI
     M_i    = [A_i*X+X'*A_i'+BF_i+BF_i' + Q.*(2*alpha),      A_i*X+BF_i-X'+Q;
              X'*A_i'+BF_i'-X+Q,             -X-X'] <=-eps*eye(2*n); 

    LMI = [LMI; M_i]; 
end

disp('Starting optimization.')
diagnostics = solvesdp(LMI,[],LMIoptions) % Solve problem 
% Find Lyapunov matrix
Q = value(Q);

% Find controller
F_0 = value(F_0); 
F_1 = value(F_1); F_2 = value(F_2); F_3 = value(F_3); 
F_4 = value(F_4); F_5 = value(F_5); F_6 = value(F_6); 
F_7 = value(F_7); F_8 = value(F_8); F_9 = value(F_9); 

X = value(X);
K_0 = F_0/X; 
K_1 = F_1/X; K_2 = F_2/X; K_3 = F_3/X; 
K_4 = F_4/X; K_5 = F_5/X; K_6 = F_6/X; 
K_7 = F_7/X; K_8 = F_8/X; K_9 = F_9/X;




function [X_struct, Lambda_struct] = findStruct(X,Lambda,eq)
% New X and Lambda
X_struct = X;
Lambda_struct = Lambda;

% Temporary Q and Lambda
X_temp = X_struct;
Lambda_temp = Lambda_struct;

while 1

    % Find zero lambda_ij and q_ij and set corresponding variables to zero
    for i=1:numel(eq)
        sv = symvar(eq(i)); % Variables in each equation
        if length(sv) == 1   % If is one variable then set corresponding variables to zero 
            hsl = has(Lambda,sv(1)); 
            hsq = has(X,sv(1)); 
            if nnz(hsl) > 0
                Lambda_struct(hsl)=0;
            end
            if nnz(hsq) > 0
                X_struct(hsq)=0;
            end
        end
    end

    % Now substitute new information in eq
    eq = subs(eq,X,X_struct);
    eq = subs(eq,Lambda,Lambda_struct);
  
    if isequal(X_struct,X_temp) && isequal(Lambda_struct,Lambda_temp)
        break;
    else
        X_temp = X_struct;
        Lambda_temp = Lambda_struct;
    end

end

end


% Find basis set of a matrix with zero elements
function L = basisSet(K)
L=[];
[ii,jj] = find(K);
for i=1:length(ii)
   z = zeros(size(K));
   z(ii(i),jj(i)) = 1;
   L = [L z];
end
end

