lambda_o=2.6*10^-6;
wo=10^(-5);
zo_free=(pi*(wo)^2)/lambda_o; %1.2*(10^(-4))
ko = 2*pi/lambda_o;
epsr=2.25;
A=10;
kz1=sqrt((ko^2)-(kx.^2));
%%%%%%%%%%%%%%%%%%%%%%%
dx=2*10^(-6);
x=-50*dx:dx:50*dx;
U = A*exp(-(x.^2)/wo^2);              %function
plot(x,U);                            
dkx=pi/(50*dx);                        
kx=-(pi/dx):dkx:(pi/dx);              
x_fft=fft(U);                         
x_fft_shift = fftshift(x_fft);        %fourier transform
figure;                               
plot(kx,abs(x_fft_shift));            %kx domain plotted
result1 = x_fft_shift .* exp(-i*kz1*0);   
figure_1 =(ifft(ifftshift(result1)))
figure;
plot(x,abs(figure_1));                 %intensity of the incident wave at z= 0(Fig1)
title('intensity incident at the initial position (Fig1)')
result2 = x_fft_shift .* exp(-i*kz1*zo_free);
figure_2a =(ifft(ifftshift(result2)));
figure;
plot(x,abs(figure_2a));              %intensity of the incident wave at zo(Fig2a)
title('intensity of the incident wave at distance equal to zo(Fig2a)')

result3 = x_fft_shift .* exp(-i*kz1*2*zo_free);   %Ei@ interface
figure_2b =(ifft(ifftshift(result3)));
figure;
plot(x,abs(figure_2b));              %intensity of the incident wave at 2zo(Fig2b)
title('intensity of the incident wave at distance equal to 2zo(Fig2b)')

%reflection&transmitted_parallel_polarization
lambda_glass=lambda_o/sqrt(epsr);
k2=2*pi/lambda_glass;
kz2=sqrt((k2^2)-(kx.^2))
n1=377;
n2=n1/sqrt(epsr);
theta_i=atan(kx/kz1);
theta_t=asin(sqrt(1/epsr)*sin(theta_i));
gamma_ref=(n2*cos(theta_t)-n1*cos(theta_i))/(n2*cos(theta_t)+n1*cos(theta_i));
tau_ref=1+gamma_ref;
zo_glass=(pi*(wo)^2)/lambda_glass;
%%%%%%
Er=gamma_ref*result3;
result4=Er.*exp(-i*kz1*zo_free);
figure_3=(ifft(ifftshift(result4)));
figure;
plot(x,abs(figure_3));                   %intensity of the transmitted wave at zo(Fig3)
title('intensity of the reflected wave at distance equal to zo(Fig3)')
%%%%%
Et=tau_ref*result3;
result5=Et.*exp(-i*kz2*zo_glass);
figure_4=(ifft(ifftshift(result5)));     
figure;
plot(x,abs(figure_4));                   %intensity of the transmitted wave at zo(Fig4)
title('intensity of the transmitted wave at distance equal to zo (Fig4)')
%%%%%
result6=Et.*exp(-i*kz2*2*zo_glass);
figure_5=(ifft(ifftshift(result6)));
figure;
plot(x,abs(figure_5));                    %intensity of the transmitted wave at zo(Fig5)
title('intensity of the transmitted wave at distance equal to 2zo (Fig5)')