%% tracking parameters
%{
load('tracking_vars','x','xd','zx','y','yd','zy');
ref_x = [x; xd; zx];
ref_y = [y; yd; zy];
clear x xd zx y yd zy
%}

%% General parameters
delta = 0.01;       % time step
w = 0.08/2;         % half foot size = 0.04

%{
x(1) = 0.0;         % xc
xd(1) = 0.0;        % xc_dot
zx(1) = 0.0;

y(1) = 0;           % yc
yd(1) = 0;          % yc_dot
zy(1) = 0;
%}

% %{
shifting_x = 0*100;
shifting_y = 0*100;
x(1) = shifting_x;         % xc
xd(1) = 0.0;               % xc_dot
zx(1) = shifting_x;
x_int(1) = 0;

y(1) = shifting_y;         % yc
yd(1) = 0;                 % yc_dot
zy(1) = shifting_y;
y_int(1) = 0;
%}

warning('off','all');

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

% additionalFirstStepDuration = 10;           % why wait for a while at first???
additionalFirstStepDuration = 0;

% one step duration                         % why zeros at first step???
fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);

% fs_matrix is the preassigned footsteps matrix
for i = 1:12
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    % [old sequence; single support has same value; 
    %  double support has an interpolated value (linear)]
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];
end

% why this deletion???
fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

% ZMP constraints are the boundaries of the support polygon at time step
%{
zx_min = fs_sequence_x - w;
zx_max = fs_sequence_x + w;

zy_min = fs_sequence_y - w;
zy_max = fs_sequence_y + w;
%}

% %{
zx_min = fs_sequence_x - w + shifting_x;
zx_max = fs_sequence_x + w + shifting_x;

zy_min = fs_sequence_y - w + shifting_y;
zy_max = fs_sequence_y + w + shifting_y;
%}

% make support polygon of stand still contain both feet
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;

%% repeated computations

% system discretization
ch = cosh(omega*delta);     % ( e^x + e^-x )/2
sh = sinh(omega*delta);     % ( e^x - e^-x )/2

A_upd = [ch      , sh/omega, 1-ch;
         omega*sh, ch      , -omega*sh;
         0       , 0       , 1];
B_upd = [delta-sh/omega;
         1-ch;
         delta];

Ts=0.01;

% %{
% k_lqr = dlqr(A_upd, B_upd, diag([100,100,100]), 1);     % not working
% k_lqr = dlqr(A_upd, B_upd, diag([1,1,1]), 1);     % not working
% k_lqr = dlqr(A_upd, B_upd, diag([0.1,0.1,0.1]), 1);     % working
% k_lqr = dlqr(A_upd, B_upd, diag([0,10^-10,10^-10]), 10^-10);    % better
k_lqr = dlqr(A_upd, B_upd, diag([0,0.1,0]), 1); 

phi = A_upd - B_upd * k_lqr;
% disp(eig(phi))
%}

%{
K_lqt = dare(A_upd,B_upd,diag([10^2,0,10^4]),1);
%}

%% disturbance

D1 = [sh/omega;ch-1;0];                      % D = [1;0;0]
D2 =[ch/(omega^2)-1/(omega^2);sh/omega;0];    % D = [0;1;0]
D3 = [delta-sh/omega; 1-ch; delta];          % D = [0;0;1]
D_par_unc=[0;0.0094;0.01]; %coming from sys=ss(A,B,[],[0;0.045;0]); and  sys=c2d(sys,delta);
% %{
% D3 = 1.5*D3;
% D2 = 0.2*D2;
% D1 = 0.06*D1*0;
% D3 = 0.5*D3;
% D2 = 0.1*D2;
% D1 = 0.06*D1;
D3 = D3;
D2 = D2;
D1 = 0.03*D1;
D_par_unc=(0.1^2+2*0.1*omega)*D_par_unc;
d = 1;
%}

%{
d = 0;
D3 = 0*D3;
D2 = 0*D2;
D1 = 0*D1;
%}

D = 0*D1 + D2 + D3;
D=D_par_unc;
d_vector = ones(N,1) * d;

% my H,Q,S,T,D
R = 10;
R_bar = diag(ones(N,1)*R);
num_states = size(A_upd, 1);
% Q = eye(num_states);
Q = eye(num_states)*0;
Q_bar = [];
for j = 1:N
    Q_bar = blkdiag(Q_bar, Q);
end

S_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        % num_states*(j-1)+1 = j*num_states-2
        S_bar(num_states*(j-1)+1:j*num_states,k) = A_upd^(j-k)*B_upd;
    end
end

T_bar = zeros(N*num_states, num_states);
for j = 1:N
    T_bar(num_states*(j-1)+1:j*num_states, :) = A_upd^j;
end

D_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        D_bar(num_states*(j-1)+1:j*num_states,k) = A_upd^(j-k)*D;
    end
end

H = 2*(R_bar + S_bar'*Q_bar*S_bar);
F = 2*S_bar'*Q_bar*T_bar;
L = 2*S_bar'*Q_bar*D_bar;

C = [0 0 1];
S_bar_y = zeros(N, N);
for j = 1:N
    for k = 1:j
        S_bar_y(j,k) = C*A_upd^(j-k)*B_upd;
    end
end

% %{
F_bar = zeros(N, N);
for j = 1:N
    for k = 1:j
        if j==k
            F_bar(j,k) = 1;
        else
            F_bar(j,k) = -k_lqr*phi^(j-k-1)*B_upd;
        end
    end
end
%}

T_bar_y = zeros(N, num_states);
for j = 1:N
    T_bar_y(j, :) = C*A_upd^j;
end

% %{
T_F_bar = zeros(N, num_states);
for j = 1:N
    T_F_bar(j, :) = -k_lqr*phi^(j-1);
end
%}

% K = [1 1 1];
K = [0 0 1];
D_bar_y = zeros(N, N);
D_bar_y_abs = zeros(N, N);

%{
for j = 1:N-1
    for k = 1:j
%         D_bar_y_abs(j+1,k) = abs(K*A_upd^(j-k)*D);
        D_bar_y(j+1,k) = K*A_upd^(j-k)*D;
    end
end
%}

% %{
for j = 1:N-1
    for k = 1:j
%         D_bar_y_abs(j+1,k) = abs(K*phi^(j-k)*D);
        D_bar_y(j+1,k) = K*phi^(j-k)*D;
    end
end
%}

constraints_restrict = D_bar_y*d_vector;
% constraints_restrict_abs = D_bar_y_abs*d_vector;
% constraints_restrict = constraints_restrict_abs;
% disp([constraints_restrict constraints_restrict_abs])
 err=0.1;
 ch = cosh((omega+err*omega)*delta);     % ( e^x + e^-x )/2
 sh = sinh((omega+err*omega)*delta);
 A_upd_real = [ch, sh/(omega+err*omega), 1-ch; (omega+err*omega)*sh, ch, -(omega+err*omega)*sh; 0, 0, 1];
 B_upd_real = [delta-sh/(omega+err*omega); 1-ch; delta];
 phi_real = A_upd_real - B_upd_real * k_lqr;
%  D2 =[ch/(omega^2)-1/(omega^2);sh/omega;0];
%% Solve
for i = 1:500

    P = delta*tril(ones(N,N));      % tril: Lower triangular part of matrix

% %{
    f_x = F*[x(i); xd(i); zx(i)] + L*d_vector;
    f_y = F*[y(i); yd(i); zy(i)] + L*d_vector;
%}

%{
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.03
            constraints_restrict(j) = 0.03;
        end
    end

    A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict];

    Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) / (1-exp(-omega*delta*N));

    % x(i) is xc, xd(i) is xc_dot
    beq_x = x(i) + xd(i)/omega - zx(i);
    beq_y = y(i) + yd(i)/omega - zy(i);
    
    zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);
    zd_y = quadprog(H,f_y,A,b_y,Aeq,beq_y);

    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*zd_x(1) + D*d;
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*zd_y(1) + D*d;
    
    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

%{
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.03
            constraints_restrict(j) = 0.03;
        end
    end

    A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict];
    
    zd_x = quadprog(H,f_x,A,b_x);
    zd_y = quadprog(H,f_y,A,b_y);

    u_lqt_x = -1*B_upd'*K_lqt*([x(i); xd(i); zx(i)] - ref_x(:,i));
    u_lqt_y = -1*B_upd'*K_lqt*([y(i); yd(i); zy(i)] - ref_y(:,i));
    
    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(u_lqt_x+zd_x(1)*0) + D*d;
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*(u_lqt_y+zd_y(1)*0) + D*d;
    
    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

% %{
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.032
            constraints_restrict(j) = 0.032;
        end
    end

    A = [S_bar_y*F_bar; -S_bar_y*F_bar];
    b_x = [ zx_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar)*[x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar)*[x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar)*[y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar)*[y(i); yd(i); zy(i)] - constraints_restrict];

    c_u_x = quadprog(H,f_x,A,b_x);
    c_u_y = quadprog(H,f_y,A,b_y);
    
    if err==0
        swap=1;
    else
        swap=0;
    end
    x_updated = phi_real*[x(i); xd(i); zx(i)] + B_upd*c_u_x(1) + 0.12*sin(6.28/100)*D*(d);
    y_updated = phi_real*[y(i); yd(i); zy(i)] + B_upd*c_u_y(1) + 0.12*sin(6.28/100)*D*(d);

    zd_x = F_bar*c_u_x + T_F_bar*[x(i); xd(i); zx(i)];
    zd_y = F_bar*c_u_y + T_F_bar*[y(i); yd(i); zy(i)];

    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

    clf
    hold on

    rect_x = [w,w,-w,-w,w]+shifting_x;
    rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w]+shifting_y;
    plot(rect_x,rect_y,'m','lineWidth',2,'HandleVisibility','off');
    
    rect_x = [w,w,-w,-w,w]+shifting_x;
    rect_y = [w,-w,-w,w,w]+shifting_y;
    
    nPlottedFootsteps = 10;
    
    for j = 1:nPlottedFootsteps
        rect_x = [w,w,-w,-w,w]+shifting_x;
        rect_y = [w,-w,-w,w,w]+shifting_y;
        h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2,'HandleVisibility','off');
    end
    
    h2 = plot(x,y,'r','lineWidth',2);
    h3 = plot(zx,zy,'b','lineWidth',2);
    h4 = plot(z_pred_x(1:end),z_pred_y(1:end)','g','lineWidth',2);
    legend('CoM', 'ZMP', 'ZMPprediction') 
    title('Parametric uncertainty disturbance action with restricted constraints')
    axis equal
    axis([-0.2 2 -0.5 0.5]+shifting_x)
    legend
    xlabel('x [m]')
    ylabel('y [m]') 
    grid on
    
    drawnow
%     break;
end