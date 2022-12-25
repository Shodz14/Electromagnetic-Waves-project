%clc;
lambda_o=2.6*10^-6;
wo=10^(-5);
zo_free=(pi*(wo)^2)/lambda_o; %1.2*(10^(-4))
ko = 2*pi/lambda_o;
epsr=2.25;
A=10;
%       %       %       %
dx=2*10^(-6);                           % since pi/k = 1.3*10^(-6) therfore we chose  {dx=2*10^(-6); > pi/k }
x=-50*dx:dx:50*dx;
u = A*exp(-(x.^2)/wo^2);                  %input function
U = u.^2;                                 %Squaring the field equation to give the intensity
plot(x,U);                                % plotting input X vs intensity
xlabel('X(m)'); ylabel('intensity');                       
title ('plotting input');

dkx=pi/(50*dx);                        
kx=-(pi/dx):dkx:(pi/dx);                  %Kx domain
kz1=sqrt((ko^2)-(kx.^2));

x_fft=fft(U);                         
x_fft_shift = fftshift(x_fft);            %fourier transform

figure;                               
plot(kx,abs(x_fft_shift));                %Kx domain plotted rad/m
xlabel('Kx (rad/m)')
title ('Kx domain plotted');

result1 = x_fft_shift .* exp(-i*kz1*0);   % Kx * T.F.
figure_1 =(ifft(ifftshift(result1)))      % 
figure;
plot(x,abs(figure_1));                    %intensity of the incident wave at initial position (Fig1)
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the incident wave at the initial position (Fig1)')

result2 = x_fft_shift .* exp(-i*kz1*zo_free);
figure_2a =(ifft(ifftshift(result2)));
figure;
plot(x,abs(figure_2a));                   %intensity of the incident wave at Zo (Fig 2a)
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the incident wave at distance equal to Zo(Fig 2a)')

result3 = x_fft_shift .* exp(-i*kz1*2*zo_free);   %Ei at interface
figure_2b =(ifft(ifftshift(result3)));
figure;
plot(x,abs(figure_2b));                    %intensity of the incident wave at 2 Zo(Fig 2b)
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the incident wave at distance equal to 2 Zo(Fig 2b)')



%---------->   Reflected & Transmitted parallel polarization          <---------
             
%----------> Variables initialization
lambda_glass = lambda_o/sqrt(epsr);
k2  = 2*pi/lambda_glass;
kz2 = sqrt((k2^2)-(kx.^2))
n1  = 377;
n2  = n1/sqrt(epsr);
theta_i   = atan(kx/kz1);
theta_t   = asin(sqrt(1/epsr)*sin(theta_i));
gamma_ref = (n2*cos(theta_t)-n1*cos(theta_i))/(n2*cos(theta_t)+n1*cos(theta_i));
tau_ref   = 1 + gamma_ref;
zo_glass  = (pi*(wo)^2)/lambda_glass;

%----------> Reflection   <-------%
Er=gamma_ref*result3;                    % Er = gamma * Ei
result4=Er .* exp(-i*kz1*zo_free);
figure_3=(ifft(ifftshift(result4)));
figure;
plot(x,abs(figure_3));                   %intensity of the Reflected wave at Zo (Fig3)
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the reflected wave at distance equal to Zo (Fig3)')



%----------> Transmission <-------%

Et=tau_ref*result3;
result5=Et.*exp(-i*kz2*zo_glass);
figure_4=(ifft(ifftshift(result5)));     
figure;
plot(x,abs(figure_4));                   %intensity of the transmitted wave at Zo(Fig4)
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the transmitted wave at distance equal to Zo (Fig4)')


result6=Et.*exp(-i*kz2*2*zo_glass);
figure_5=(ifft(ifftshift(result6)));
figure;
plot(x,abs(figure_5));                    %intensity of the transmitted wave at Zo(Fig5)
xlabel('X(m)'); ylabel('intensity');
title('intensity of the transmitted wave at distance equal to 2 Zo (Fig5)')

%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%_______________________________________________________________________________________________________________________
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


%% EXRTAS for more visualisation of wave behaviour 

%% to observe the behaviour of incident wave
figure;
plot(x,abs(figure_1),"g" ,x,abs(figure_2a),"b", x,abs(figure_2b), "r" );
xlabel('X(m)'); ylabel('intensity'); 
title('intensity of the incident wave at differnt distances initial, Zo, 2 Zo ');
legend ("incident wave at initial pos.",
"incident wave at dist. Zo", 
"incident wave at dist. 2Zo");

%%%%%%%%%%%%%%%%%%%%%%%%
%% INCIDENT VS REFLECTED
figure;
plot(x,abs(figure_1),"g" ,x,abs(figure_2a),"b", x,abs(figure_2b), "r" ,x,abs(figure_3) , "--k" );
xlabel('X(m)'); ylabel('intensity');
title('intensity of the INCIDENT VS REFLECTED waves at distances initial, Zo , 2 Zo ')
legend ("incident wave at initial pos.",
"incident wave at dist. Zo","incident wave at dist. 2 Zo"
,"reflected wave at dist. Zo" );

%%%%%%%%%%%%%%%%%%%%%%%%
%% to observe the behaviour of transmitted wave
figure;
plot(x,abs(figure_4),"-.b" , x,abs(figure_5) , "-.r"); 
xlabel('X(m)'); ylabel('intensity');
title('intensity of the transmitted wave at distances Zo , 2 Zo ')
legend ("transmitted wave at dist. Zo", "transmitted wave at dist. 2 Zo");

%%%%%%%%%%%%%%%%%%%%%%%%
%% %% INCIDENT VS TRANSMITTED
figure;
plot(x,abs(figure_1),"g" ,x,abs(figure_2a),"b", x,abs(figure_2b), "r" ,
 x,abs(figure_4),"-.b" , x,abs(figure_5) , "-.r"); 
xlabel('X(m)'); ylabel('intensity');
title('intensity of the INCIDENT VS TRANSMITTED wave at distances Zo , 2 Zo ')
legend ("incident wave at initial pos.","incident wave at dist. Zo", 
"incident wave at dist. 2 Zo",
"transmitted wave at dist. Zo", "transmitted wave at dist. 2 Zo");
%%%%%%%%%%%%%%%%%%%%%%%%


%%%% <---------------> INCIDENT VS REFLECTED VS TRANSMITTED <--------------->
figure;
plot(x,abs(figure_1),"g" ,x,abs(figure_2a),"b", x,abs(figure_2b), "r" 
,x,abs(figure_3) ,"--k" 
,x,abs(figure_4),"-.b" , x,abs(figure_5) , "-.r"); 
xlabel('X(m)'); ylabel('intensity');
title('intensity of the INCIDENT VS REFLECTED VS TRANSMITTED  waves at different distances initial , Zo , 2 Zo ')
legend ("incident wave at initial pos.","incident wave at dist. Zo", 
"incident wave at dist. 2 Zo", "reflected wave at dist. Zo" ,
"transmitted wave at dist. Zo", "transmitted wave at dist. 2 Zo");

