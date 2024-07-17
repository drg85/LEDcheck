
%Load in solar spectra data!
clear
load('DATASPEC.MAT')


%%Values for LED%%
u = 453; %enter wavelength;;
FWHM = 20; %enter FWHM;
s = FWHM./2.355; %gets standard deviation 
E = 90.0; %joules per cm2, irradaince 

%integral
x = round(u - 3*s):round(u + 3*s);
c = 1/sqrt(2*pi*s*s);
disx = c.*exp((-(x-u).^2)./(2*s*s));
cdfx = 0.5.*(1 + erf((x-u)./(sqrt(2)*s)));

i = 1;

while i < length(x)
    lw  = x(i);
    up = x(i+1);
    
    h = find(wavelength >= lw & wavelength <= up); 
    irrad = mean(irradiance(h)); %wm2 sunlight at this band 
    irrad = 10^-4.*irrad; %wcm2; 
    energyled = abs(E.*(cdfx(i+1)-cdfx(i)));
    %energyledcheck = abs(E.*(disx(i+1)-disx(i)));
    t(i) = energyled./irrad; 
    %tcheck(i) = energyledcheck./irrad;
   i = i +1;  
end

x2 = x(2:length(x)); 

HRS = max(t)./3600


%makeplot
x_fill = [x, fliplr(x)]; % Create a vector for x data including the mirrored part
y_fill = [disx, zeros(size(disx))];
hf = max(disx)/2;
xf = (u - (FWHM)/2):1:(u + (FWHM)/2);
yf = 0.*xf + hf; 


close all
subplot(1,2,1);
%plot(x,disx);
hold on
%fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3);
area(x,disx); alpha(0.3); 
hold on; 
plot(xf,yf)
title('LED Spectra')
xlabel('Wavelength (nm)')
ylabel('Proportion')
%plot spectral region%
xb = find(wavelength >= (u - 3*s) & wavelength <= (u+3*s));
yb = irradiance(xb);
xb = wavelength(xb); 

subplot(1,2,2);
plot(wavelength,irradiance);
hold on
area(xb,yb); alpha(0.3); 
xlim([0,800])
title(['Equivalent solar hours: ' num2str(HRS) ] )
xlabel('Wavelength (nm)')
ylabel ('Spectral irradince')