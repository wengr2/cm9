
%% Continuum Mechanics, Exercise 9
% Paul Kulyk
% Raphael Wenger
%
% paul.kulyk@students.unibe.ch
% raphael.wenger@students.unibe.ch
%

% Due May 9, 201
%% Include the predefined tensor math functions
addpath('../../Matlab/');
if exist('imported','var') ~= 1
    def_symbols;
    imported = 1;
end
%% Parameters
%Toogle symbolic values, used for testing
enableSyms = 1;
%% 8.1 Computation of all necessary variables
% define the geometry and build the tetrahedron
Xi = [[1/10 0 0]; [0 1/10 0]; [0 0 1/10]; [0 0 0]; ]';

% density
rho = 1000; % kg/m^3

%Gravity
g = 10; % m/s^2 

% Time symbols
syms T Tmax t


%%
% Rotation around e1,e2,e3 respectively
R1 = @(theta) [ [ 1 0 0 ]; [ 0 cos(theta) -sin(theta) ]; [ 0 sin(theta) cos(theta) ];];
R2 = @(theta) [ [ cos(theta) 0 -sin(theta) ]; [ 0 1 0 ]; [ sin(theta) 0 cos(theta) ];];
R3 = @(theta) [ [ cos(theta) -sin(theta) 0 ]; [ sin(theta) cos(theta) 0 ]; [ 0 0 1 ];];

%%
%Motion of the tetrahedron as in exercice 7
Tmax = 1/2;

if enableSyms == 1
    T = t;                  %Time defined as a symbol for starters
else
    T = Tmax;
end

%Rotation
Rt = R3(2*pi*T/Tmax);

%Translation
bt = [ 0; 0; 3/20*T/Tmax];

%Transformation matrix
y =@(R,x,b) R*x + b;

%%
%Contact force defined with handles
F_contact  =  @(theta) [-pi^2/15*(cos(theta)-sin(theta)); -pi^2/15*(cos(theta)+sin(theta));5/3];
F_contact0 =  @(fact,theta) fact*[sin(theta)+cos(theta);sin(theta)-cos(theta);0];

%Apply the given values
F_con  = F_contact(4*pi*T);
F_con0 = F_contact0((125+pi^2*(4+60*T))/3000,4*pi*T);

%%
% Function handle to clean up all these triple integrals
TripInt =@(fun,v1,l1,u1,v2,l2,u2,v3,l3,u3) int( int( int( fun, v1, l1, u1), v2, l2, u2), v3, l3, u3);

% Function handle for the are of faces
faceArea = @(nodeB,nodeC,nodeD) 1/2*cm.norm(cm.cross_product((nodeC-nodeD),(nodeB-nodeD)));  

%%
%Computation of the new vertices
yt = y(Rt,Xi(1:3,1:3)*b,bt); % Simple transform of b into x into y...

%%
%Center of gravtiy
V = 1/6000; %
M = 6*V*TripInt(rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);    

% Function handle for the centre of gravity
COG     =@(xx) 6*V/M*TripInt(rho*xx,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);

yc = COG(yt);

%%
%Final positions
yi(:,1) = y(Rt,Xi(:,1),bt);
yi(:,2) = y(Rt,Xi(:,2),bt);
yi(:,3) = y(Rt,Xi(:,3),bt);
yi(:,4) = y(Rt,Xi(:,4),bt);

%Get the normals to the faces, centers and area
Ai(1)=faceArea(Xi(:,2),Xi(:,3),Xi(:,4));
Ai(2)=faceArea(Xi(:,3),Xi(:,4),Xi(:,1));
Ai(3)=faceArea(Xi(:,4),Xi(:,1),Xi(:,2));
Ai(4)=faceArea(Xi(:,1),Xi(:,2),Xi(:,3));

Ait(1)=faceArea(yi(:,2),yi(:,3),yi(:,4));
Ait(2)=faceArea(yi(:,3),yi(:,4),yi(:,1));
Ait(3)=faceArea(yi(:,4),yi(:,1),yi(:,2));
Ait(4)=faceArea(yi(:,1),yi(:,2),yi(:,3));

%%
%Use the provided function for the normals to the surface
[xnormi xcenti] = cm.get_tetra_normal(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
[ynormi ycenti] = cm.get_tetra_normal(yi(:,1),yi(:,2),yi(:,3),yi(:,4));

%%
%Computation of the Area-weighted normas
for i = 1:4
    Aini(:,i) = Ai(i)*ynormi(:,i);
end

%%
%Only non-zero scalar product is at i=4
%because all other faces are on the systems plane.

%Only non-singular matrix is @i=4
i=4;
wt = -4*cm.invert(cm.scalar_product(ycenti(:,i),Aini(:,i))*I+cm.dyadic_product11(ycenti(:,i),Aini(:,i)))*(F_con0-cm.cross_product(yc,F_con));

% Antisymmetric stress tensor W
for i = 1:4
    WtiAini(:,i) = -1/2*cm.cross_product(wt,Aini(:,i));
end

%Symmetric stress tensor
if enableSyms == 1
    syms T11 T12 T13 T21 T22 T23 T31 T32 T33
    Tt = [T11, T12, T13; T12, T22, T23; T13, T23, T33];  %Taking advantage of the symmetry              
else
    Tt = zeros(3)
end

% Forces vector on the vertices
for i=1:4
    ift(:,i) = -1/3*(Tt*Aini(:,i)+ WtiAini(:,i))+1/4*F_con;
end  
%% 9.1 (1) Cauchy stress
%Compute the sum of the forces
Tt = zeros(3);
for i=1:4
    fit(:,i) = Tt*ynormi(:,i)-1/2*cm.cross_product(wt,ynormi(:,i))+ 1/(4*Ai(i))*F_con;
end  

%Substitue for the plots
tmpfit=cm.roundDecimals(double(subs(fit,t,Tmax)),2);
tmpyi=subs(yi,t,Tmax);
tmpycenti=subs(ycenti,t,Tmax);

% Plot the tetras
figure(1)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
% Add the forces vectors
scaleFact = 0.001;
for i=1:4
    hold on
    cm.plot_vector(tmpycenti(:,i),tmpycenti(:,i)+tmpfit(:,i)*scaleFact,2,'red')
end


%% 9.2 (2) Nominal stress vectors
% Compute the sum of the moments
for i=1:4
    fiptfiAi(:,i) = fit(:,i)*Ait(i);
    fipt(:,i)=fiptfiAi(:,i)/Ai(i);
end  

%Substitue for the plots
tmpfiptfiAi=cm.roundDecimals(double(subs(fiptfiAi,t,Tmax/2)),2);


% Plot the tetras
figure(2)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
% Add the forces vectors
scaleFact = 0.2;
for i=1:4
    hold on
    cm.plot_vector(xcenti(:,i),xcenti(:,i)+tmpfiptfiAi(:,i)*scaleFact,2,'red')
end

%% 9.3 (3) Material stress vectors
%Composition of the transformation
F = Rt;

for i=1:4
    fist(:,i) = cm.transpose2(F)*fipt(:,i);
end  

%Substitue for the plots
tmpfist=cm.roundDecimals(double(subs(fipt,t,Tmax/2)),2);


% Plot the tetras
figure(3)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
% Add the forces vectors
scaleFact = 0.001;
for i=1:4
    hold on
    cm.plot_vector(xcenti(:,i),xcenti(:,i)+tmpfist(:,i)*scaleFact,2,'green')
end

%% 9.4 (4) Cauchy Stress Tensor
T   =@(t)   1e6*[...
                [ 2*t*(t+1)*(1+2*t)*cos(4*pi*t)^2   t*(t+1)*(1+2*t)*sin(8*pi*t)     0 ];...
                [ t*(t+1)*(1+2*t)*sin(8*pi*t)       2*t*(t+1)*(1+2*t)*sin(4*pi*t)^2 0 ];...
                [ 0                                 0                               0 ];...
            ];

% Defined in the problem
t_max = 3/8;

% Generate the Cauchy stress tensor at tmax/4
T_4 =   T(t_max/4);

%Compute the additionnal deformation
Ut=(t/t_max)*cm.dyadic_product11(e1,e1)+I ;
%Composition of the rotation and dilations
F = Rt*Ut;

%Update the transformation with the dilation. 
%Not sur what tmax, t etc. we should use anymore...
for i=1:4
    yi(:,i) = y(F,Xi(:,i),bt);
end
[ynormi ycenti] = cm.get_tetra_normal(yi(:,1),yi(:,2),yi(:,3),yi(:,4));

%Compute the stess vectors on the faces
for i=1:4
    fit(:,i) = T_4*ynormi(:,i);
end 

%Substitue for the plots
tmpfit=cm.roundDecimals(double(subs(fit,t,t_max/4)),2);
tmpyi=subs(yi,t,t_max/4);
tmpycenti=subs(ycenti,t,t_max/4);

% Plot the tetras
figure(5)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
hold on
% Add the forces vectors. Somehow not working...
scaleFact = 0.000001;
for i=1:4
    cm.plot_vector(tmpycenti(:,i),tmpycenti(:,i)+tmpfit(:,i)*scaleFact,2,'red')
    hold on
end

%% 9.4 (5) Nominal stress vectors
% What in the 7 hells is that star supposed to be?
Fstar=cm.det2(F)*cm.invert(cm.transpose2(F));
P = T_4*Fstar;

% Compute the nominal stress vectors
for i=1:4
    fipt(:,i)=P*xnormi(:,i);
end  

%Substitue for the plots
tmpfipt=cm.roundDecimals(double(subs(fipt,t,0)),2);

% Plot the tetras
figure(2)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
% Add the forces vectors
scaleFact = 0.0000002;
for i=1:4
    hold on
    cm.plot_vector(xcenti(:,i),xcenti(:,i)+tmpfipt(:,i)*scaleFact,2,'red')
end
%% 9.4 (6) Material stress vector
S = P/F;


for i=1:4
    fist(:,i) = S*xnormi(:,i);
end  

%Substitue for the plots
tmpfist=cm.roundDecimals(double(subs(fipt,t,0)),2);


% Plot the tetras
figure(3)
cm.plot_tetra_dual(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4),tmpyi(:,1),tmpyi(:,2),tmpyi(:,3),tmpyi(:,4))
% Add the forces vectors
scaleFact = 0.000001;
for i=1:4
    hold on
    cm.plot_vector(xcenti(:,i),xcenti(:,i)+tmpfist(:,i)*scaleFact,2,'green')
end

%% 
%Finish
%publish('exercice9.m')
