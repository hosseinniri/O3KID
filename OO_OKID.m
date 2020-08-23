clc
clear all


%  Y = |         |  q: number of output
%      |  l x q  |  l:length of data
%      |         |  
%Y=Y(:,1)

% p  : dim of V 
% n  : DOF of system

% shoud be : l-p>pq
% p >> n :p > 20*n
% l= 10000,nd =2, n=3


p = 100;
%order ==>> n
n = 4;

%% create state space model
Ts = 0.001;
s=tf('s');
sysc = (2000)/((s+20+40i)*(s+20-40i)*(s+40+100i)*(s+40-100i));
%step(sysc);

sysd = c2d(sysc,Ts);
Num = sysd.Numerator{1};
Den = sysd.Denominator{1};

% Num = sysc.Numerator{1};
% Den = sysc.Denominator{1};

[A,B,C,D] = tf2ss(Num,Den);
clear sysd Num sysc Den
%% simulation with state space

T = 0:Ts:10-Ts;    % simulation time = 10 seconds
U = ones(1,length(T));    % no input
X0 = [0.1 0.1 0.1 0.1 ];    % initial conditions of the three states
sys = ss(A,B,C,D,Ts);     % construct a system model
%;0 1 0 1
[Y, Tsim, X] =lsim(sys, U, T, X0);    % simulate and plot the response (the output)
%figure
%plot(Tsim,Y)

%% Identification Auto Regressive

[l,q] = size(Y);
v=zeros(q*p,l-p);

for k=p+1:l
    %disp(k);
    for i=1:p
        for j=1:q
                v((q*(i-1)+j),k)=Y(k-i,j);
        end
    end
end

disp("      constructed V")

clear k i j

v_c = v(:,p+1:l);
y_c = Y(p+1:l,:)';

clear v

phi = y_c*pinv(v_c);
 
E_e = y_c - phi*v_c ;
Y_e = y_c -E_e;

%% convert markov parameter from bar form to inovation form

psi=zeros(size(phi));

psi(:,1:q) = phi(:,1:q);

N = p;

for j=2:N
    psi(:,(j-1)*q+1:j*q) = phi(:,(j-1)*q+1:j*q);
    for h=1:j-1
        psi(:,(j-1)*q+1:j*q) = psi(:,(j-1)*q+1:j*q)+ phi(:,(h-1)*q+1:h*q )* psi(:,(j-h-1)*q+1:(j-h)* q);
    end
end

%% STEP 2-a : Deterministic intersection
E_Y = zeros(2*q,l);
E_Y(:,p+1:l)=[E_e;Y_e];

% H = | H(p)|
%     | H(f)|


i = 2* (n)/(q);
%l = size(E_Y,2); 

H = zeros(4*i*q,(l-2*i-p+1));


for ix = 1:2*i*q
   for iy = 1:(l-2*i-p+1)
      H(2*(ix-1)*q+1:2*ix*q,iy) = E_Y(:,p+ix+iy-1);
      %disp(p+ix+iy-1);
   end
end
clear ix iy E_Y
disp('                   ')
disp('      constructed H')

%%calculate svd of H

[U,S,V_T] = svd(H);

S11 = S(1:2*q*i+n,1:2*q*i+n);
U11 = U(1:2*q*i,1:2*q*i+n);
U12 = U(1:2*q*i,2*q*i+n+1:4*q*i);
U21 = U(2*q*i+1:4*q*i,1:2*q*i+n);
U22 = U(2*q*i+1:4*q*i,2*q*i+n+1:4*q*i);

UU = U12'*U11*S11;
[U_qr,~,V_T_qr] = svd(UU);

Uq=U_qr(:,1:n);

LHS = Uq'*U12'*U(2*q+1:2*q*(i+1),1:2*q*i+n)*S11;
RHS =[Uq'*U12'*U(1:2*q*i,1:2*q*i+n)*S11;U(2*q*i+1:2*q*i+q,1:2*q*i+n)*S11];

sol = LHS/RHS;
 
A_DI=sol(1:n,1:n);
C_DI = (U(2*q*i+q+1:2*q*(i+1), 1:2*q*i+n)*S11)/(Uq'*U12'*U(1:2*q*i,1:2*q*i+n)*S11);

clear U S V_T S11 U11 U12 U21 U22 UU U_qr Sq V_T_qr Uq LHS RHS sol
%zplane(eig(A_DI),eig(A))
disp('                   ')
disp('      DI completed ')
%% STEP 2-b : Deterministic Projection

% construct hankel matrix of Y and E
Ye = zeros(q,l);
Ee = zeros(q,l);

Ye(:,p+1:l) = Y_e;
Ee(:,p+1:l) = E_e;

Y_h = zeros(i,l-i-p+1);
E_h = zeros(i,l-i-p+1);


for ix=1:i
    for iy=1:l-i-p+1
        Y_h((ix-1)*q+1:ix*q,iy)=Ye(:,p+ix+iy-1);
        E_h((ix-1)*q+1:ix*q,iy)=Ee(:,p+ix+iy-1);
    end
end

% project Y_h on E_h
I = eye(l-i-p+1);

Y_p_E =Y_h*(I-E_h'*pinv(E_h*E_h')*E_h);

[U,S,V_T]=svd(Y_p_E);

S1 = S(1:n,1:n);
U1 = U(:,1:n);
U2 = U(:,n:q*i);
V1_T = V_T(1:n,:);
V2_T = V_T(n+1:(l-i-p+1),:);


U1_L=U1(1:size(U1,1)-q,:);
U1_H=U1(q+1:size(U1,1),:);

C_DP = U1(1:q,:);
A_DP = pinv(U1_L)*U1_H;

clear Ye Ee Y_h E_h Y_p_E U S V_T S1 U1 U2 V1_T V2_T U1_L U1_H
disp('                   ')
disp('      DP completed ')
%% STEP 2-c : ERA method on Markove parameter

%Construct H0 and H1
H0 = zeros((N/2)*q,(N/2)*q);
H1 = zeros((N/2)*q,(N/2)*q);

clear i j

for j=1:(N/2)
    for i=1:(N/2)
        k=i+j-1;
        H0((i-1)*q+1:i*q,(j-1)*q+1:j*q) = psi(:,(k-1)*q+1:k*q);
    end
end
for j=1:(N/2)
    for i=1:(N/2)
        k=i+j;
        H1((i-1)*q+1:i*q,(j-1)*q+1:j*q) = psi(:,(k-1)*q+1:k*q);
    end
end

% compute svd of H0
[U,S,V] = svd(H0);

S1 = S(1:n,1:n);
U1 = U(:,1:n);
V1_T = V(1:n,:);

% compute observability and controllability matrices
Ob = U1*(S1^(0.5));
Co = (S1^(0.5))*V1_T;
%clculate A and C and K matrices
C_ERA =Ob(1:q,:);
K_ERA =Co(:,1:q);
A_ERA = (S1^(-0.5))*U1'*H1*V1_T'*(S1^(-0.5));

clear i j k U S V S1 U1 V1_T 
disp('                   ')
disp('      ERA completed ')
%% ERA-DC method
H0_DC = H0*H0';
H1_DC = H1*H1';

[U_DC,S_DC,V_DC] = svd(H0_DC);

S1_DC = S_DC(1:n,1:n);
U1_DC = U_DC(:,1:n);
V1_T_DC = V_DC(1:n,:);

Ob_DC = U1_DC*(S1_DC^(0.5));
Co_DC = (S1_DC^(-0.5))*U1_DC'*H0;

%clculate A and C and K matrices
C_ERA_DC =Ob_DC(1:q,:);
K_ERA_DC =Co_DC(:,1:q);
A_ERA_DC = (S1_DC^(-0.5))*U1_DC'*H1_DC*V1_T_DC'*(S1_DC^(-0.5));

clear H0 H1 H0_DC H1_DC U_DC S_DC V_DC S1_DC U1_DC V1_T_DC Ob_DC Co_DC 
disp('                   ')
disp('      ERA-DC completed ')
%% visualise Results

eig_A = eig(A);
%eig_A = (1/Ts)*log(eig_A);

real_A = real(eig_A);
imag_A = imag(eig_A);

eig_A_DI = eig(A_DI);
%eig_A_DI = (1/Ts)*log(eig_A_DI);

real_eig_A_DI = real(eig_A_DI);
imag_eig_A_DI = imag(eig_A_DI);

eig_A_DP = eig(A_DP);
%eig_A_DP = (1/Ts)*log(eig_A_DP);

real_eig_A_DP = real(eig_A_DP);
imag_eig_A_DP = imag(eig_A_DP);

eig_A_ERA = eig(A_ERA);
%eig_A_ERA = (1/Ts)*log(eig_A_ERA);

real_eig_A_ERA = real(eig_A_ERA);
imag_eig_A_ERA = imag(eig_A_ERA);

eig_A_ERA_DC = eig(A_ERA_DC);
%eig_A_ERA_DC = (1/Ts)*log(eig_A_ERA_DC);

real_eig_A_ERA_DC = real(eig_A_ERA_DC);
imag_eig_A_ERA_DC = imag(eig_A_ERA_DC);

plot(real_A,imag_A,'r*',real_eig_A_DI,imag_eig_A_DI,'bo',real_eig_A_DP,imag_eig_A_DP,'gx',real_eig_A_ERA,imag_eig_A_ERA,'kx',real_eig_A_ERA_DC,imag_eig_A_ERA_DC,'m*')

xlabel("Real");
ylabel("Imag");
legend("original system","OKID-DI","OKID-DP","OKID-ERA","OKID-ERA-DC")

disp('                   ')
disp('      visualize completed ')

