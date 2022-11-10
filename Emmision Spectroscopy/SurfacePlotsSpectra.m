%% Surface Plots of Spectra

clear;clc;

%% Spectra Surface - Fix Pressure, Vary Voltage
pressure = "1 Torr";
files = "2022-09-15 - Spectra/Calibrated Spectra/" + pressure + "/*.txt";
fileList = dir(files);
power = [5,10,15,20,25,30,35,40,45,50];

for i=1:length(fileList)
    fileName = append(fileList(i).folder,'\',fileList(i).name);
    data = readtable(fileName);
    lambda = data.(1);
    I(:,i) = data.(2);
end
V = [75.8630,91.9614,107.5586,124.9212,141.3685,152.4274,165.7229,182.9509,188.4703,205.9964];

figure
surf(lambda,V,I')
colormap("cool")
title(pressure,'FontSize',44)
ylabel("Voltage (V)",'FontSize',22)
xlabel("Wavelength (nm)",'FontSize',22)
zlabel("Radiance (W/m^3/sr)",'FontSize',22)

figure
hold on
for i = 1:10
    plot(lambda, I(:,i))
end
title(pressure)
ylabel("Radiance (W/m^3/sr)")
xlabel("Wavelength (nm)")
legend("75","90","110","125","140","150","165","180","190","200")

%% Spectra Surface - Fix Power, Vary Pressure

voltage = "100 V";
P = [100E-3,250E-3,500E-3,1,5,10];

I2(:,1) = readtable("2022-09-15 - Spectra\Calibrated Spectra\100 mTorr\100mTorr-10W.txt").(2);
I2(:,2) = readtable("2022-09-15 - Spectra\Calibrated Spectra\250 mTorr\250mTorr-10W.txt").(2);
I2(:,3) = readtable("2022-09-15 - Spectra\Calibrated Spectra\500 mTorr\500mTorr-10W.txt").(2);
I2(:,4) = readtable("2022-09-15 - Spectra\Calibrated Spectra\1 Torr\1Torr-15W.txt").(2);
I2(:,5) = readtable("2022-09-15 - Spectra\Calibrated Spectra\5 Torr\5Torr-20W.txt").(2);
I2(:,6) = readtable("2022-09-15 - Spectra\Calibrated Spectra\10 Torr\10Torr-20W.txt").(2);
lambda = readtable("2022-09-15 - Spectra\Calibrated Spectra\10 Torr\10Torr-20W.txt").(1);

figure
surf(lambda,P,I2')
colormap("cool")
title(voltage)
ylabel("Pressure (Torr)")
xlabel("Wavelength (nm)")
zlabel("Radiance (W/m^3/sr)")

figure
hold on
for i = 1:6
    plot(lambda, I2(:,i))
end
title(voltage)
ylabel("Radiance (W/m^3/sr)")
xlabel("Wavelength (nm)")
legend("100 mTorr","250 mTorr","500 mTorr","1 Torr","5 Torr","10 Torr")

%% Total Population Density 

n = readmatrix("2022-09-15 - Spectra/Total Population Surface Plot.csv");
n = [0.2090    0.2360    0.2880    0.3160    0.3670    0.3990    0.4380    0.5100    0.6080    0.6570
    0.3680    0.5670    0.6520    0.7540    0.8270    0.9270    0.9740    1.0400    1.1500    1.4400
    0.5150    0.7740    0.9500    1.1200    1.2500    1.3500    1.5300    1.8700    1.7400    1.7900
    1.1500    1.3000    1.5200    1.8700    1.7500    1.7800    2.1200    2.0400    2.1000    1.7500
    0.7740    1.2700    1.8500    2.9600    3.5500    4.2700    4.9700    5.5900    6.4400    6.7900
    0.7730    0.9950    1.3700    1.9400    2.4400    3.7100    5.0800    6.4200    7.2900    7.9100]*1E13*3;
W = [5,10,15,20,25,30,35,40,45,50];
V = [83.9887,104.8642,124.6258,134.8502,150.6055,164.8989,174.7557,185.5948,193.8169,207.0062;
    90.4187,110.4760,124.1029,135.8078,150.9806,162.6632,170.5978,185.0311,193.4532,206.0659;
    84.4622,102.1209,121.9359,136.1874,145.1525,157.3846,174.7630,187.4510,197.1321,209.9522;
    75.8630,91.9614,107.5586,124.9212,141.3685,152.4274,165.7229,182.9509,188.4703,205.9964;
    64.6519,73.3336,83.1963,99.5738,111.9207,130.0122,140.3667,146.6210,155.9750,155.4913;
    77.4767,82.4186,86.9436,98.6111,108.6263,119.9059,126.9737,128.7102,133.5577,136.5453];
P = [0.1,0.25,0.5,1,5,10];
figure
surf(V,P,n)
xlabel('Voltage P2P (V)','FontSize',22)
ylabel('Pressure (Torr)')
zlabel('Population Density (#/m^3)','FontSize',22)
title('Ar(4p) Total Population Density','FontSize',44)
