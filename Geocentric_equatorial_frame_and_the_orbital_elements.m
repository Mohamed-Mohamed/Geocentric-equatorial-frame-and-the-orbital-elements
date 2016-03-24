%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg


% this program is used for showing a satalite about the earth
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
hz=15000; % simulation frequancy
step_sim=150; % simulation step
%% Intial condition
M=[5.972*10^24,1000];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[8000e3;0;0]+[3.844e8;0;0]*0;                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0;7.5;0]*1e3+[0;1.023;0]*1e3*0;                 % velocity of M(2)
% image (inter your URL of the image)
image_file = 'D:\4th year of Aerospace\1st\Orbital Mechanics\AER-427, Orbital Mechanics, Mohamed Mohamed Elsayed,SCE 2, BN 13  By MATLAB\week 11\satalite about earth with tracking/earth.jpg';
%% RK45 parameter
tf=9000+30*24*3600*0;   % final time of soution
dt=100+3600*0;            % rough time step
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
sol(1:12,1)=X0;
order=12;
%% RK4
for n=1:length(0:dt:tf)
    b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0; ...
        0,0,0,0,0,0,0,0,0,1,0,0; ...
        0,0,0,0,0,0,0,0,0,0,1,0; ...
        0,0,0,0,0,0,0,0,0,0,0,1;...
        -b,0,0,b,0,0,0,0,0,0,0,0; ...
        0,-b,0,0,b,0,0,0,0,0,0,0; ...
        0,0,-b,0,0,b,0,0,0,0,0,0; ...
        -c,0,0,c,0,0,0,0,0,0,0,0; ...
        0,-c,0,0,c,0,0,0,0,0,0,0; ...
        0,0,-c,0,0,c,0,0,0,0,0,0 ];
    [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
    sol(1:12,n+1)=XX(1:12,2);
end
R1_x=sol(1,:)';
R1_y=sol(2,:)';
R1_z=sol(3,:)';
R2_x=sol(4,:)';
R2_y=sol(5,:)';
R2_z=sol(6,:)';
V1_x=sol(7,:);
V1_y=sol(8,:);
V1_z=sol(9,:);
V2_x=sol(10,:);
V2_y=sol(11,:);
V2_z=sol(12,:);
%% projected path of the satalite on the earth
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
for P=1:length(R1_x)
    points = intersectLineSphere([0,0,0,R2_x(P)-R1_x(P)+(P-1)*erot*dt,R2_y(P)-R1_y(P)+(P-1)*erot*dt,R2_z(P)-R1_z(P)+(P-1)*erot*dt], [0,0,0,6400e3]);
    XX1(P,1)=points(2,1);
    YY1(P,1)=points(2,2);
    ZZ1(P,1)=points(2,3);
end
%% calculating Right ascension and Declination
for k=1:length(R1_x)
    [ RA(k), Dec(k) ] = lmn2RaDec( [R2_x(k), R2_y(k), R2_z(k)]-[R1_x(k), R1_y(k), R1_z(k)] );
end
%% plotting
%--------------------------------------------------------------------------------------------------------------------------------------------------------
fig=figure();
% slider
menu = uicontrol('Parent',fig,'Style','popupmenu','String',{'omega';'i';'w';'Exit'},'Units','centimeters' ,'Position',[17.5,0.25,3,0.5]);
slider = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',360);
% slider_alpha = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',360);
% slider_beta = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0.75,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',180);
% slider_gamma = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,1.5,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',360);
% figure(fig);
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
M1=1;
p1=0;
i=0; w=0; omega=0;
% AA_old1=0;
% AA_old2=0;
% AA_old3=0;
view(3);
while M1 == 1
    figure(fig);
    S = get(menu,'value');
    if S == 1 
        omega = get(slider,'value');
    elseif S == 2
        i = get(slider,'value')/360*180;
    elseif S == 3
        w = get(slider,'value');
    elseif S == 4
        M1=0;
    end
    [ Q ] = RM ( omega, i, w, '313');
    VV=[[R2_x,R2_y,R2_z]-[R1_x,R1_y,R1_z]]*Q;
    VV1=[XX1,YY1,ZZ1]*Q;
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = 4.89496121282306 - (p1-1)*erot*dt; % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    grid on;
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    hold;
    globe = surf(y, x, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    hold all;
    plot3(VV1(:,1),VV1(:,2),VV1(:,3),'r','LineWidth',2);
    % M2 trajectory
    hold all;
    plot3(I(1)*ones(length(R1_x),1)+VV(:,1)-R1_x,I(2)*ones(length(R1_y),1)+VV(:,2)-R1_y,I(3)*ones(length(R1_z),1)+VV(:,3)-R1_z,'g','LineWidth',2);
    % M2 start
    hold all;
    plot3(VV(1,1)-R1_x(1),VV(1,2)-R1_y(1),VV(1,3)-R1_z(1),'*','color',[.9,0.2,0.5],'LineWidth',10);
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title(['Satellite about Earth w = ' num2str(w) ' i = ' num2str(i) ' \Omega = ' num2str(omega)  ],'Fontsize',18);
    % Equatorial plane
    theta_vec=linspace(0,2*pi,30);
    r_vec=0:1e6:1.8e7;
    [theta_mat, r_mat]=meshgrid(theta_vec, r_vec);
    [x_mat, y_mat]=pol2cart(theta_mat, r_mat);
    z_mat=r_mat*0;
    surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.2);
    legend('Projected path of the satalite on the earth','Satellite trajectory','Satellite perigee','Equatorial plane','location','north');
    p1=p1+1;
    pause(1/500);
end
close all;