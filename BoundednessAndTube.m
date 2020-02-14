% it works  better if we have a longer in time set omega. Anyway it is not
% necessary to have it until the end of the prediction horizon!
%% General parameters
delta = 0.01;
w = 0.08/2;%half of the width of the box

x(1) = 0.0;
xd(1) = 0.0;
zx(1) = 0.0;

y(1) = 0;
yd(1) = 0;
zy(1) = 0;

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

additionalFirstStepDuration = 10;

fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);%x component ofsupport the center of the polygon will be stored in here
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);%y component ofsupport the center of the polygon will be stored in here

for i = 1:12 %number of the steps
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];%unit step range that 
end
%building the box
fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

zx_min = fs_sequence_x - w;
zx_max = fs_sequence_x + w;  

zy_min = fs_sequence_y - w;
zy_max = fs_sequence_y + w;
%beginning support polygon rectangle
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;

%% Compute matrices
p = ones(N,1);
P = delta*tril(ones(N,N));
A = [P;-P];

Vu = [];
Vs = [];

 %% discretization; this part could have been computed outside the loop
 ch = cosh(omega*delta);
 sh = sinh(omega*delta);
 A_upd = [ch, sh/omega, 1-ch; omega*sh, ch, -omega*sh; 0, 0, 1];
 B_upd = [delta-sh/omega; 1-ch; delta];
 C_forComTracking=[1,0,0];
 C_forTubeMPC=[0,0,1];
 p=[0.8 0.81 0.7];
 p=[0.4001 0.51 0.45];%% yeaaah!!
Q_r=0.1*eye(3);
R=1;
N_2=[0;0;0];
[F,S,e] = dlqr(A_upd,B_upd,Q_r,R,N_2);
%  F=place(A_upd,B_upd,p);
 K=(C_forComTracking*(eye(3)-(A_upd-B_upd*(F)))^(-1)*B_upd)^-1; % tracking matrix
 load('variables','zx_model','zy_model', 'x_model','y_model','xd_model','yd_model');
 D=[0;0.0094;0.01]; d=0.009;
%% Tube MPC restriction
ymax=0.04;   % Support Polygon for the ZMP
ymin=-0.04;
P_y=Polyhedron('A',[1  ;-1 ],'B',[ymax;-ymin]);
wmax=0.011;   % Bounded Disturbances Set
wmin=-0.011;
wmax=0.1330012;   % Bounded Disturbances Set which can surely be
                      % rejected in a very conservative way!
wmin=-0.1330012;
% it almost reaches the meddle and a solution is found almost to the center
% of the support polygon!
t=wmax;
% t=0.008000012;
P_w=Polyhedron('A',[1  ;-1 ],'B',[wmax;-wmin]);
% Initializations
array_y=[];  
array_y_err=[];  
array_y_old=[];
Rj=0*C_forTubeMPC*D*P_w;
Rj_err=0*C_forTubeMPC*D*P_w;
Phi= A_upd-B_upd*F; % Stabilized State Transition Matrix
array_Rj=[];
array_Rj_err=[];
y_M=[];
y_m=[];
y_M_err=[];
y_m_err=[];
e=0;

for i=1:N % Loop to Compute the Erosion
Rj_err=Rj_err.minVRep()+C_forTubeMPC*(A_upd)^(i-1)*D*P_w; %.minVRep
array_Rj_err=[array_Rj_err;Rj_err];
scop=Rj_err; % Polygon to Remove from P_y .minVRep()
array_y_err=[array_y_err;P_y-scop]; % Pontryagin Difference (Erosion)
if array_y_err(i).isEmptySet() % 36 for Phi
e=e+1;
y_M_err=[y_M_err; array_y_err(i-e).V(1,1)];%(i)
y_m_err=[y_m_err;array_y_err(i-e).V(2,1)];   
else
y_M_err=[y_M_err; array_y_err(i).V(1,1)];%(i)
y_m_err=[y_m_err;array_y_err(i).V(2,1)];
% end
end
end
limit_h=i-e
%{
array_y=[];  
array_y_err=[];  
array_y_old=[];
Rj=0*C_forTubeMPC*D*P_w;
Rj_err=0*C_forTubeMPC*D*P_w;
Phi= A_upd-B_upd*F; % Stabilized State Transition Matrix
array_Rj=[];
array_Rj_err=[];
y_M=[];
y_m=[];
y_M_err=[];
y_m_err=[];

g=1;
while g
    e=0;
wmax=wmax-0.0002;   % Bounded Disturbances Set
wmin=-0.012+0.0002;
P_w=Polyhedron('A',[1  ;-1 ],'B',[wmax;-wmin]);
t_star=wmax;
for i=1:N % Loop to Compute the Erosion
Rj_err=Rj_err.minVRep()+C_forTubeMPC*(A_upd)^(i-1)*D*P_w; %.minVRep
array_Rj_err=[array_Rj_err;Rj_err];
scop=Rj_err; % Polygon to Remove from P_y .minVRep()
array_y_err=[array_y_err;P_y-scop]; % Pontryagin Difference (Erosion)
if array_y_err(i).isEmptySet() % 36 for Phi
e=e+1; 
end
end
(N-e)
if (N-e) >= N-limit_h
 g=0;   
end
end
t_star
pause
 %}
restriction=ymax*ones(N,1)-y_M_err;
Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) / (1-exp(-omega*delta*N));
%%
input_x=[];
input_y=[];
Q = 0;
P = delta*tril(ones(N,N));
H = Q*(P'*P) + eye(N);
% %% Solve
for i = 1:691

    
    f_x = Q*P'*(- fs_sequence_x(i:i+N-1) + zx(i));
    f_y = Q*P'*(- fs_sequence_y(i:i+N-1) + zy(i));
    beq_x = x(i) + xd(i)/omega - zx(i);
    beq_y = y(i) + yd(i)/omega - zy(i);
     
%   { 
 b_x = [ zx_max(i:i+N-1) - zx(i)-restriction;...       %-restriction
        - zx_min(i:i+N-1) + zx(i)-restriction];
    b_y = [ zy_max(i:i+N-1) - zy(i)-restriction;...
        - zy_min(i:i+N-1) + zy(i)-restriction];
    %}
% % {
%     b_x = [ zx_max(i:i+N-1) - zx(i);...      
%         - zx_min(i:i+N-1) + zx(i)];
%     b_y = [ zy_max(i:i+N-1) - zy(i);...
%         - zy_min(i:i+N-1) + zy(i)];
 %}  
    
    
     zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);
     zd_y = quadprog(H,f_y,A,b_y,Aeq,beq_y);
   
    
    
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    input_x=[input_x,zd_x(1)];
    input_y=[input_y,zd_y(1)];
    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(zd_x(1))-D*t;
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*(zd_y(1))-D*t;
    %updates
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
    %graph from now on
    clf
    hold on

    rect_x = [w,w,-w,-w,w];
    rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w];
    plot(rect_x,rect_y,'m','lineWidth',2,'HandleVisibility','off');
    
    rect_x = [w,w,-w,-w,w];
    rect_y = [w,-w,-w,w,w];
    
    nPlottedFootsteps = 12;
    
    for j = 1:nPlottedFootsteps
        rect_x = [w,w,-w,-w,w];
        rect_y = [w,-w,-w,w,w];
        h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2,'HandleVisibility','off');
    end
    
    h2 = plot(x,y,'r','lineWidth',2);
    h3 = plot(zx,zy,'b','lineWidth',2);
    h4 = plot(z_pred_x(1:end),z_pred_y(1:end)','g','lineWidth',2);
    legend('CoM', 'ZMP', 'ZMPprediction') 
    title('Disturbance action with restricted constraints plus boundedness condition')
    axis equal
    axis([-0.2 2.5 -0.5 0.5])
    legend
    xlabel('x [m]')
    ylabel('y [m]') 
    grid on
    
    drawnow
end
%%

save('In', 'input_x','input_x');
close all
load('In', 'input_x','input_x');


%{
x(1) = 0.0;
xd(1) = 0.0;
zx(1) = 0.0;

y(1) = 0;
yd(1) = 0;
zy(1) = 0;
num_states=3;
C=C_forTubeMPC;
T_bar_y = zeros(N, num_states);
for j = 1:N
    T_bar_y(j, :) = C*A_upd^j;
end
C = [0 0 1];
S_bar_y = zeros(N, N);
for j = 1:N
    for k = 1:j
        S_bar_y(j,k) = C*A_upd^(j-k)*B_upd;
    end
end
H = 2*(eye(N,N)+1*eye(N,N));
A = [S_bar_y; -S_bar_y];
Q = 0;
P = delta*tril(ones(N,N));
H = Q*(P'*P) + eye(N);
%}
