%clc, clear all
a1 = 0.01:0.01:0.4; %Size of Scatterers (sm)
lya1 = 0.2; %Wavelength (sm)
f = 30./lya1; %Frequency (GHz)
X1 = 2.*pi.*a1./lya1;
t1 = 20; %Temperature in degrees Celsius
lyas1 = 1.4662.*exp(-0.0634.*t1)+0.000136.*t1.*t1-0.027296.*t1+1.8735116; %лямбда эстое
e011 = 5.5; %эпсилон
es11 = 0.00081.*t1.*t1-0.40885.*t1+88.2; %эпсилон эстое
e111 = e011+((es11-e011)./(1+((lyas1./lya1).*(lyas1./lya1)))); %эпсилон первое
e211 = ((es11-e011).*(lyas1./lya1))./(1+((lyas1./lya1).*(lyas1./lya1))); %эпсилон второе
m111 = sqrt((e111+sqrt(e111.*e111+e211.*e211))./2);
m211 = sqrt((-e111+sqrt(e111.*e111+e211.*e211))./2);
m1 = m111-1i.*m211;
Y1 = X1.*m1;
k1 = 2.*pi./lya1;

for kk1 = 1:length(X1)
    XX1 = X1(kk1); % берётся значение к-го элемента массива Х, с ним идёт цикл по Nmax
    Nmax1 = round(XX1+4.*(XX1.^(1./3))+2);
    massPsix1(1) = (sin(XX1)./XX1)-cos(XX1);
    massPsix1(2) = (3./XX1).*((sin(XX1)./XX1)-cos(XX1))-sin(XX1);
    massHex1(1) = (cos(XX1)./XX1)+sin(XX1);
    massHex1(2) = (3./XX1).*((cos(XX1)./XX1)+sin(XX1))-cos(XX1);
    massPsiSx1(1) = sin(XX1)-(1/XX1).*((sin(XX1)./XX1)-cos(XX1));
    massPsiSx1(2) = (sin(XX1)./XX1)-cos(XX1)-(2./XX1).*((3./XX1).*((sin(XX1)./XX1)-cos(XX1))-sin(XX1));
    massHeSx1(1) = cos(XX1)-(1/XX1).*((cos(XX1)./XX1)+sin(XX1));
    massHeSx1(2) = (cos(XX1)./XX1)+sin(XX1)-(2./XX1).*((3./XX1).*((cos(XX1)./XX1)+sin(XX1))-cos(XX1));
    for n1 = 2:Nmax1
        massPsix1(n1+1) = ((2.*n1+1)./XX1).*massPsix1(n1)-massPsix1(n1-1);
        massHex1(n1+1) = ((2.*n1+1)./XX1).*massHex1(n1)-massHex1(n1-1);
        massPsiSx1(n1+1) = massPsix1(n1)-((n1+1)./XX1).*massPsix1(n1+1);
        massHeSx1(n1+1) = massHex1(n1)-((n1+1)./XX1).*massHex1(n1+1);
        
            massPsix1_mat(1,kk1) = massPsix1(1);
            massPsix1_mat(n1,kk1) = massPsix1(n1);
            massHex1_mat(1,kk1) = massHex1(1);
            massHex1_mat(n1,kk1) = massHex1(n1);
            massPsiSx1_mat(1,kk1) = massPsiSx1(1);
            massPsiSx1_mat(n1,kk1) = massPsiSx1(n1);
            massHeSx1_mat(1,kk1) = massHeSx1(1);
            massHeSx1_mat(n1,kk1) = massHeSx1(n1);
    end
end

for kk1 = 1:length(Y1)
    YY1 = Y1(kk1); % берётся значение к-го элемента массива Х, с ним идёт цикл по Nmax
    XX1 = X1(kk1);
    Nmax1 = round(XX1+4.*(XX1.^(1./3))+2);
    massPsiy1(1) = (sin(YY1)./YY1)-cos(YY1);
    massPsiy1(2) = (3./YY1).*((sin(YY1)./YY1)-cos(YY1))-sin(YY1);
    massHey1(1) = (cos(YY1)./YY1)+sin(YY1);
    massHey1(2) = (3./YY1).*((cos(YY1)./YY1)+sin(YY1))-cos(YY1);
    massPsiSy1(1) = sin(YY1)-(1./YY1).*((sin(YY1)./YY1)-cos(YY1));
    massPsiSy1(2) = (sin(YY1)./YY1)-cos(YY1)-(2./YY1).*((3./YY1).*((sin(YY1)./YY1)-cos(YY1))-sin(YY1));
    massHeSy1(1) = cos(YY1)-(1./YY1).*((cos(YY1)./YY1)+sin(YY1));
    massHeSy1(2) = (cos(YY1)./YY1)+sin(YY1)-(2./YY1).*((3/YY1).*((cos(YY1)./YY1)+sin(YY1))-cos(YY1));
    for n1 = 2:Nmax1
        massPsiy1(n1+1) = ((2.*n1+1)./YY1).*massPsiy1(n1)-massPsiy1(n1-1);
        massHey1(n1+1) = ((2.*n1+1)./YY1)*massHey1(n1)-massHey1(n1-1);
        massPsiSy1(n1+1) = massPsiy1(n1)-((n1+1)./YY1)*massPsiy1(n1+1);
        massHeSy1(n1+1) = massHey1(n1)-((n1+1)./YY1)*massHey1(n1+1);
        
            massPsiy1_mat(1,kk1) = massPsiy1(1);
            massPsiy1_mat(n1,kk1) = massPsiy1(n1);
            massHey1_mat(1,kk1) = massHey1(1);
            massHey1_mat(n1,kk1) = massHey1(n1);
            massPsiSy1_mat(1,kk1) = massPsiSy1(1);
            massPsiSy1_mat(n1,kk1) = massPsiSy1(n1);
            massHeSy1_mat(1,kk1) = massHeSy1(1);
            massHeSy1_mat(n1,kk1) = massHeSy1(n1);
    end
end

MAT_DZTX1 = massPsix1_mat + 1i.*massHex1_mat;
MAT_DZTSX1 = massPsiSx1_mat + 1i.*massHeSx1_mat;
MAT_DZTY1 = massPsiy1_mat + 1i.*massHey1_mat;
MAT_DZTSY1 = massPsiSy1_mat + 1i.*massHeSy1_mat;

MAT_An1 = (massPsix1_mat.*massPsiSy1_mat-m1.*massPsiy1_mat.*massPsiSx1_mat)./(MAT_DZTX1.*massPsiSy1_mat-m1.*massPsiy1_mat.*MAT_DZTSX1);
MAT_Bn1 = (m1.*massPsix1_mat.*massPsiSy1_mat-massPsiy1_mat.*massPsiSx1_mat)./(m1.*MAT_DZTX1.*massPsiSy1_mat-massPsiy1_mat.*MAT_DZTSX1);
MAT_An1(isnan(MAT_An1))=0;
MAT_Bn1(isnan(MAT_Bn1))=0;

h = zeros(1,6);
h(1,1) = 0.5;
h(1,2) = h(1,1)/2;
for i = 3:length(h)
    h(1,i) = h(1,i-1)/2;
end

theta_1 =  0:h(1,6)*(pi/180):(pi/2); %h(1,6) самое точное ошибка 0.015(1.5%)
for kk = 1:length(theta_1)
     theta1 = theta_1(kk); % берётся значение к-го элемента массива theta, с ним идёт цикл по Nmax
     massPi1(1) = 1;
     massTau1(1) = cos(theta1);
     massPi1(2) = 3.*cos(theta1);
     massTau1(2) = 6.*cos(theta1).*cos(theta1)-3;
    for n1 = 3:Nmax1
            massPi1(n1) = ((1./(n1-1)).*(2.*n1-1)).*cos(theta1).*massPi1(n1-1)-(n1./(n1-1)).*massPi1(n1-2);
            massTau1(n1) = n1.*cos(theta1).*massPi1(n1)-(n1+1).*massPi1(n1-1);
            
            massPi1_mat(1,kk) = massPi1(1);
            massPi1_mat(2,kk) = massPi1(2);
            massPi1_mat(n1,kk) = massPi1(n1);
            massTau1_mat(1,kk) = massTau1(1);
            massTau1_mat(2,kk) = massTau1(2);
            massTau1_mat(n1,kk) = massTau1(n1);
    end
end

for j = 1:length(a1)
        massPi1_mat_s(:,:,j) = massPi1_mat(:,:);
        massTau1_mat_s(:,:,j) = massTau1_mat(:,:);
end
for j = 1:length(theta_1)
        MAT_An1_s(:,j,:) = MAT_An1(:,:);
        MAT_Bn1_s(:,j,:) = MAT_Bn1(:,:);
end

matAP_BT = MAT_An1_s.*massPi1_mat_s + MAT_Bn1_s.*massTau1_mat_s;
matAT_BP = MAT_An1_s.*massTau1_mat_s + MAT_Bn1_s.*massPi1_mat_s;

for j = 1:Nmax1
    matN(j,:) = (2.*j+1)./(j.*(j+1));
end
for j = 1:length(a1)
    for jj = 1:length(theta_1)
        matN_L(:,jj,j) = matN(:);
    end
end

matUSf1 = matAP_BT.*matN_L;
matUSf2 = matAT_BP.*matN_L;


ff1 = squeeze(sum((1i./k1).*matUSf1(:,:,:)));
ff2 = -squeeze(sum((1i./k1).*matUSf2(:,:,:)));
    
ff1_mod = abs(ff1); ff2_mod = abs(ff2);
white1 = ff1_mod.^2 + ff2_mod.^2;
ff2_conj = conj(ff2);
black1 = real(ff1.*ff2_conj);
S_1 = (3/4).*white1 + (1/2).*black1;

I1 = 0:10:100;
%I1 = 10;

    for j=1:length(I1)
        S__1(:,:,j) = S_1(:,:)';
    end
    
    for j=1:length(a1)
        I1_a(j,:) = I1(:);
    end
    
    for j=1:length(I1)
        a1_I(:,j) = a1(:);
    end
    
        %w1 = 5.674.*(10.^-5).*(I1_a.^0.324).*((2.*a1_I).^-1.75).*exp(-98.5.*((2.*a1_I).^2.25).*(I1_a.^-0.522)); %distribution of Best
        %w1 = 7.3656.*(10^-5).*(I1_a.^0.411).*((2.*a1_I).^-1.71).*exp(-116.97.*((2.*a1_I).^2.29).*(I1_a.^-0.456)); %distribution of Lous-Parsons
        %w1 = 0.16.*exp(-41.*(2.*a1_I).*(I1_a.^-0.21)); %distribution of Marshal-Palmer
        w1 = 1.25.*(I1_a.^-0.139).*exp(-((113.*a1_I)./(I1_a.^0.175))); %distribution of Polyakovoy-Shifrin
    
    for j = 1:length(theta_1)
        w1_2(:,j,:) = w1(:,:);
    end
    
    func_und_int = S__1.*w1_2;%alpha/theta/I
    %fff1 = squeeze(sum(func_und_int)).*0.01; %взят интеграл по а1 %+
    fff = trapz(a1,func_und_int,1);
    fg = squeeze(fff);
    %______________________________________________________________
    
    ksi = 0.001:0.001:1;
     
    for j = 1:length(ksi)
         ffff(:,j,:) = fg(:,:);
    end
    fffff = squeeze(ffff);
    
    % константы:
    g1 = 10; % g1 = 1-10 deg
    tettat0 = g1.*(pi/180);
    alphaT = 4.*log(2)./(tettat0.^2);
    g2 = 10; % g2 = 1-10 deg
    tettar0 = g2.*(pi/180);
    alphaR = 4.*log(2)./(tettar0.^2);
    
    kh = 1.5823; % 0.2: 1.5823; 0.36: 1.2101; 0.8: 0.3895;
    alph = 0.6494; % 0.2: 0.6494; 0.36: 0.7052; 0.8: 0.8817;
    GAMMA = kh.*(I1.^(alph));
    L = 10;% 1 km = 10^5 sm; L = 1-10;
    
    % 2D сетка координат
    [K,T,GAMMAA] = meshgrid(ksi,theta_1,GAMMA);% [..GAMMAA]=..GAMMA
    % расчет подынтегральной функции: 
    
    F = fffff.*((10^5).*pi.*L.*(T./(K.^2))).*exp(-(T.*T).*alphaT).*exp(-(T.*T).*alphaR.*(((1-K)./K).^2)).*exp(-(T.*T).*((0.1.*log(10).*L.*GAMMAA)./2).*(((1-K)./K)));
    % метод трапеций
    i1st = trapz(theta_1,F,1);
    i_tr = squeeze(squeeze(trapz(ksi,i1st)));
    
    %loglog(I1,i_tr,'*'); %x,s,o,*,d,.,h,:,--,b--o,
    
    %plot(I1,i_tr,'x'); %x,s,o,*,d,.,h,:,--,b--o,
    %xlabel('Инстенсивность дождя, мм/ч');
    %ylabel('Дисперсия');
    %title('Засивимость дисперсии от интенсивности, разные распределения');
    text(11,1.1,'Длина трассы(L) 10 км, полуширины по 10 град, длина волны 0.2 см');
    %grid on
    %legend('длина волны 0.2 см','длина волны 0.36 см','длина волны 0.8 см');
    %legend('Бест','Лоус-Парсонс','Маршал-Пальмер', 'Полякова-Шифрин');
    hold on
