%% shear flow calculations
% Everything is in MKS units
% Dhairya Dosi,Aug 2023
clearvars
t = 0.0005;  % thickness in m
V = 35; 
Izz = 0;
Iyz = 0;
q{1,1} = 0;
 
%% defination of the airfoil
% Imported as CSV file from airfoiltools.com
% Chord needs to be in MM.
% make sure size of 'm' and 'n' array is the same
% m array is for X coordinates
m = [290    289.55485   288.22201   286.0096    282.93125   279.00581   274.25706   268.71458   262.41201   255.38821   247.6861    239.35295   230.44009   221.00233   211.09767   200.78672   190.13299   179.20231   168.06109   156.77835   145.42311   134.06497   122.77411   111.61926   100.66915   89.99106    79.65024    69.71049    60.233  51.15658    42.41221    34.15446    26.53877    19.70956    13.78689    8.86095 4.98597 2.18254 0.44138 -0.2697 0   1.16348 3.1291  5.83016 9.2075  13.21385    17.82108    23.02484    28.8463 35.3278 42.52705    50.5035 59.30935    68.7648 78.69237    89.03087    99.71592    110.68169   121.86003   133.18163   144.57689   155.97476   167.30477   178.49703   189.48194   200.19135   210.55972   220.52209   230.01756   238.98697   247.37493   255.12953   262.20292   268.55102   274.13468   278.91939   282.87528   285.9777    288.2078    290];
% n array is fot Y coordinates
% last element should be put to zero to ensure closed loop
n = [0  0.000911    0.003625    0.008097    0.014256    0.022002    0.031218    0.041763    0.053488    0.06623 0.079825    0.094105    0.108892    0.124016    0.139293    0.154541    0.169577    0.184205    0.198224    0.211424    0.223599    0.234526    0.243997    0.251807    0.257758    0.261682    0.26343 0.262894    0.259997    0.254576    0.245729    0.232632    0.215006    0.193108    0.167658    0.139699    0.110417    0.080971    0.052339    0.02523 0   -0.022524   -0.041798   -0.058273   -0.0725 -0.085077   -0.096561   -0.107433   -0.11805    -0.128615   -0.139119   -0.149251   -0.158311   -0.165387   -0.170311   -0.17313    -0.173925   -0.172811   -0.169937   -0.165459   -0.159555   -0.152407   -0.1442 -0.135114   -0.125326   -0.115008   -0.104325   -0.093435   -0.082493   -0.071656   -0.061068   -0.050886   -0.041255   -0.032323   -0.024238   -0.017127   -0.011121   -0.006328   -0.002836   0];
%defining the polygon
polygon = polyshape(m,n);
plot(polygon,"LineWidth",0.5)

%% finding centriod
[Cx, Cy] = centroid(polygon);

 
%% idealisation
area1{1,1} = 2*t*sqrt((m(1,1)-m(1,80))^2+((n(1,1)-n(1,80)))^2)/2000;
area1{1,2} =[m(1,1)/1000];
area1{1,3} =[n(1,1)/1000];
for k=2:(length(m)-1)
   le =  t*(((sqrt((m(1,k+1)-m(1,k))^2+((n(1,k+1)-n(1,k)))^2)/6000)*(2+((n(1,k+1)-Cy)/(n(1,k)-Cy))))+(sqrt((m(1,k)-m(1,k-1))^2+((n(1,k)-n(1,k-1)))^2)/6000)*(2+((n(1,k)-Cy)/(n(1,k-1)-Cy))));
   area1{k,1} =le;
   area1{k,2} =[m(1,k)/1000];
   area1{k,3} =[n(1,k)/1000];
end
area1{80,1} = 2*t*sqrt((m(1,80)-m(1,1))^2+((n(1,80)-n(1,1)))^2)/2000;
area1{80,2} =[m(1,80)/1000];
area1{80,3} =[n(1,80)/1000];

%% finding out moment of inertia
for k =1:(length(m))
    Izz = area1{k,1}*(area1{k,3}^2) + Izz;
    Iyz = area1{k,1}*(area1{k,3}*area1{k,2}) + Iyz;
end
%Izz = 5*10^-8;
Forceterm = V/Izz; 
%% find shear flow
for k =2:(length(m))
    q{k,1} = (Forceterm*area1{k,1}*(area1{k,3}))/1000 + q{k-1,1};
end
 
%% angle of twist is zero
AOT = 0;
for k =1:(length(m)-1)
    AOT = q{k+1,1}*sqrt((m(1,k+1)-m(1,k))^2+((n(1,k+1)-n(1,k)))^2)/1000 + AOT;
end

%finding total parameter of the airfoil
param = 0;
for k =1:(length(m)-1)
    param  = sqrt((m(1,k+1)-m(1,k))^2+((n(1,k+1)-n(1,k)))^2)/1000 + param;
end
param  = .29*2.2;
q{1,1} = AOT/param;

% final results
for k =1:(length(m)-1)
     q{k+1,1} = q{k+1,1} - q{1,1};
end

%% now we visyualise the data generated 
% Warning. it`s just a visualisation tool. 
% for exact values rely on 'q' colmn vector define a part above
%scale factor (choose the one giving the best fit
cheat_code = 1;
scale_factor = 1;
% offset line in respect to the shear flo
for k =1:(length(m)-1)
    plot([m(1,k),m(1,k+1)],[n(1,k),n(1,k+1)]) %plotting airfoil
    xlabel('Chord Position (mm)') 
    hold on;
    sl = (n(1,k+1)-n(1,k))/(m(1,k+1)-m(1,k)); %offsetting the line
    y = n(1,k);
    x = m(1,k);
    a = (sl^2+1);
    b = 2*(y-sl*x);
    c = q{k,1}*scale_factor;
    dy = roots([a b c]);
    if real(dy(1,1)) > 0
        d1y = abs(dy(1,1));
    else
        d1y = -abs(dy(1,1));
    end
    d1x = sl*d1y;
    plot([m(1,k)+d1x,m(1,k+1)+d1x],[n(1,k)+d1y,n(1,k+1)+d1y])
    ylabel('Shear flow(N/M)') 
    title('Shear flow distribution on NACA 23015 wing root without spar contact')
end
