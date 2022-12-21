Y1 = X1.*m1;
k1 = 2.*pi./lya1;

for kk1 = 1:length(X1)
    XX1 = X1(kk1); % Xk, cycle by Nmax
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
    YY1 = Y1(kk1); % Xk, cycle by Nmax
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
