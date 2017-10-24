% Function to calculate lyapunov exponent of a system defined by input
% parameters. Here f1 and f2 replace d1 and d2 to avoid clashes with
% tangent vectors. Steps per cycle is set to 25, and numver of iterations
% to be skipped or counted can be set below, though current values are
% reccomended. Method is too complicated to explain here, see report.

function [lExponents] = LyapunovExponent(a, b, c, f1, f2, irrationalFreq)
    
    % Initialising values needed
    stepSize = 50;
    h = (2 * pi) / stepSize;
    hTau = h * irrationalFreq;
    
    skipIterations = 10^4;
    numIterations  = 10^5;
    
    % Initial conditons for orbits and tangent vectors
    p    = 1;
    s    = 0;
    z    = 0;
    zTau = 0;
    
    dx1    = 10 ^ (-2);
    dv1    = 0;
    dz1    = 0;
    dzTau1 = 0;
    
    dx2    = 0;
    dv2    = 10 ^ (-2);
    dz2    = 0;
    dzTau2 = 0;
    
    dx3    = 0;
    dv3    = 0;
    dz3    = 10 ^ (-2);
    dzTau3 = 0;
    
    dx4    = 0;
    dv4    = 0;
    dz4    = 0;
    dzTau4 = 10 ^ (-2);
    
    
    % error size we will keep our tangent line distance at.
    initCondSize1 = norm([dx1, dv1, dz1, dzTau1]);
    initCondSize2 = norm([dx2, dv2, dz2, dzTau2]);
    initCondSize3 = norm([dx3, dv3, dz3, dzTau3]);
    initcondSize4 = norm([dx4, dv4, dz4, dzTau4]);
    
    % counters to track progression through driving cycle
    hCount    = 0;
    hTauCount = 0;
    
    % iterations to skip transients
    for i = 1 : skipIterations 
        
        %---------------------------------------------------------------
        K1 = s;
        L1 = a*s + b*p + c*p*p*p + f1*cos(z) + f2*cos(zTau);
  
        M1 = dv1;
        N1 = a*dv1 + b*dx1 + 3*c*p*p*dx1 - f1*sin(z)*dz1 -...
            f2*sin(zTau)*dzTau1;
        
        F1 = dv2;
        G1 = a*dv2 + b*dx2 + 3*c*p*p*dx2 - f1*sin(z)*dz2 -...
            f2*sin(zTau)*dzTau2;
        
        H1 = dv3;
        J1 = a*dv3 + b*dx3 + 3*c*p*p*dx3 - f1*sin(z)*dz3 -...
            f2*sin(zTau)*dzTau3;
        
        U1 = dv4;
        V1 = a*dv4 + b*dx4 + 3*c*p*p*dx4 - f1*sin(z)*dz4 -...
            f2*sin(zTau)*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K2 = s + (0.5*h*L1);
        L2 = a*(s + (0.5*h*L1)) + b*(p + (0.5*h*K1)) +...
            c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(p + (0.5*h*K1)) +...
            f1*cos((z + (0.5*h))) + f2*cos(zTau + (0.5*hTau));
        
        M2 = dv1 + (0.5*h*N1);
        N2 = a*(dv1 + (0.5*h*N1)) + b*(dx1 + (0.5*h*M1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx1 + (0.5*h*M1)) -...
            f1*sin((z + (0.5*h)))*dz1 - f2*sin(zTau + (0.5*hTau))*dzTau1;
        
        F2 = dv2 + (0.5*h*G1);
        G2 = a*(dv2 + (0.5*h*G1)) + b*(dx2 + (0.5*h*F1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx2 + (0.5*h*F1)) -...
            f1*sin((z + (0.5*h)))*dz2 - f2*sin(zTau + (0.5*hTau))*dzTau2;
        
        H2 = dv3 + (0.5*h*J1);
        J2 = a*(dv3 + (0.5*h*J1)) + b*(dx3 + (0.5*h*H1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx3 + (0.5*h*H1)) -...
            f1*sin((z + (0.5*h)))*dz3 - f2*sin(zTau + (0.5*hTau))*dzTau3;
        
        U2 = dv4 + (0.5*h*V1);
        V2 = a*(dv4 + (0.5*h*V1)) + b*(dx4 + (0.5*h*U1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx4 + (0.5*h*U1)) -...
            f1*sin((z + (0.5*h)))*dz4 - f2*sin(zTau + (0.5*hTau))*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K3 = s + (0.5*h*L2);
        L3 = a*(s + (0.5*h*L2)) + b*(p + (0.5*h*K2)) +...
            c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(p + (0.5*h*K2)) +...
            f1*cos((z + (0.5*h))) + f2*cos(zTau + (0.5*hTau));
        
        M3 = dv1 + (0.5*h*N2);
        N3 = a*(dv1 + (0.5*h*N2)) + b*(dx1 + (0.5*h*M2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx1 + (0.5*h*M2)) -...
            f1*sin((z + (0.5*h)))*dz1 - f2*sin(zTau + (0.5*hTau))*dzTau1;
        
        F3 = dv2 + (0.5*h*G2);
        G3 = a*(dv2 + (0.5*h*G2)) + b*(dx2 + (0.5*h*F2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx2 + (0.5*h*F2)) -...
            f1*sin((z + (0.5*h)))*dz2 - f2*sin(zTau + (0.5*hTau))*dzTau2;
        
        H3 = dv3 + (0.5*h*J2);
        J3 = a*(dv3 + (0.5*h*J2)) + b*(dx3 + (0.5*h*H2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx3 + (0.5*h*H2)) -...
            f1*sin((z + (0.5*h)))*dz3 - f2*sin(zTau + (0.5*hTau))*dzTau3;
        
        U3 = dv4 + (0.5*h*V2);
        V3 = a*(dv4 + (0.5*h*V2)) + b*(dx4 + (0.5*h*U2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx4 + (0.5*h*U2)) -...
            f1*sin((z + (0.5*h)))*dz4 - f2*sin(zTau + (0.5*hTau))*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K4 = s + (h*L3);       
        L4 = a*(s + h*L3) + b*(p + h*K3) +...
            c*(p + h*K3)*(p + h*K3)*(p + h*K3) +...
            f1*cos((z + h)) + f2*cos(zTau + hTau);
       
        M4 = dv1 + (h*N3);
        N4 = a*(dv1 + (h*N3)) + b*(dx1 + (h*M3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx1 + (h*M3)) -...
            f1*sin((z + h))*dz1 - f2*sin(zTau + hTau)*dzTau1;
        
        F4 = dv2 + (h*G3);
        G4 = a*(dv2 + (h*G3)) + b*(dx2 + (h*F3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx2 + (h*F3)) -...
            f1*sin((z + h))*dz2 - f2*sin(zTau + hTau)*dzTau2;
        
        H4 = dv3 + (h*J3);
        J4 = a*(dv3 + (h*J3)) + b*(dx3 + (h*H3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx3 + (h*H3)) -...
            f1*sin((z + h))*dz3 - f2*sin(zTau + hTau)*dzTau3;
        
        U4 = dv4 + (h*V3);
        V4 = a*(dv4 + (h*V3)) + b*(dx4 + (h*U3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx4 + (h*U3)) -...
            f1*sin((z + h))*dz4 - f2*sin(zTau + hTau)*dzTau4;
        %---------------------------------------------------------------
        
        p = p + (1/6)*(K1 + (2*K2) + (2*K3) + K4)*h;
        s = s + (1/6)*(L1 + (2*L2) + (2*L3) + L4)*h;
        
        dx1 = dx1 + (1/6)*(M1 + (2*M2) + (2*M3) + M4)*h;
        dv1 = dv1 + (1/6)*(N1 + (2*N2) + (2*N3) + N4)*h;

        dx2 = dx2 + (1/6)*(F1 + (2*F2) + (2*F3) + F4)*h;
        dv2 = dv2 + (1/6)*(G1 + (2*G2) + (2*G3) + G4)*h;

        dx3 = dx3 + (1/6)*(H1 + (2*H2) + (2*H3) + H4)*h;
        dv3 = dv3 + (1/6)*(J1 + (2*J2) + (2*J3) + J4)*h;
        
        dx4 = dx4 + (1/6)*(U1 + (2*U2) + (2*U3) + U4)*h;
        dv4 = dv4 + (1/6)*(V1 + (2*V2) + (2*V3) + V4)*h;
        
        z = z + h;
        hCount = hCount + 1;
        
        zTau = zTau + hTau;
        hTauCount = hTauCount + 1;
        
         if hCount == stepSize
            z = 0;
            hCount = 0;
        end
        
        if hTauCount == stepSize
           zTau = mod(zTau, (2*pi));
           hTauCount = 0; 
        end
        
        % Rescaling step to prevent folding or overlfow etc.
        error1 = [dx1, dv1, dz1, dzTau1]; 
        unitError1 = error1 / sqrt(error1(1)^2 + error1(2)^2 + error1(3)^2 + error1(4)^2);
        
        error2 = [dx2, dv2, dz2, dzTau2];
        error2 = error2 - (error2(1)*unitError1(1) + error2(2)*unitError1(2) + error2(3)*unitError1(3) + error2(4)*unitError1(4))*unitError1;
        unitError2 = error2 / sqrt(error2(1)^2 + error2(2)^2 + error2(3)^2 + error2(4)^2);
        
        error3 = [dx3, dv3, dz3, dzTau3];
        error3 = error3 - (error3(1)*unitError2(1) + error3(2)*unitError2(2) + error3(3)*unitError2(3) + error3(4)*unitError2(4))*unitError2...
            - (error3(1)*unitError1(1) + error3(2)*unitError1(2) + error3(3)*unitError1(3) + error3(4)*unitError1(4))*unitError1;
        unitError3 = error3 / sqrt(error3(1)^2 + error3(2)^2 + error3(3)^2 + error3(4)^2);
        
        error4 = [dx4, dv4, dz4, dzTau4];
        error4 = error4 - (error4(1)*unitError3(1) + error4(2)*unitError3(2) + error4(3)*unitError3(3) + error4(4)*unitError3(4))*unitError3...
            - (error4(1)*unitError2(1) + error4(2)*unitError2(2) + error4(3)*unitError2(3) + error4(4)*unitError2(4))*unitError2...
            - (error4(1)*unitError1(1) + error4(2)*unitError1(2) + error4(3)*unitError1(3) + error4(4)*unitError1(4))*unitError1;
        unitError4 = error4 / sqrt(error4(1)^2 + error4(2)^2 + error4(3)^2 + error4(4)^2);
        
        newValues1 = unitError1 * initCondSize1;
        dx1 = newValues1(1); dv1 = newValues1(2); dz1 = newValues1(3); dzTau1 = newValues1(4);
        
        newValues2 = unitError2 * initCondSize2;
        dx2 = newValues2(1); dv2 = newValues2(2); dz2 = newValues2(3); dzTau2 = newValues2(4);
        
        newValues3 = unitError3 * initCondSize3;
        dx3 = newValues3(1); dv3 = newValues3(2); dz3 = newValues3(3); dzTau3 = newValues3(4);
        
        newValues4 = unitError4 * initcondSize4;
        dx4 = newValues4(1); dv4 = newValues4(2); dz4 = newValues3(3); dzTau4 = newValues4(4);
        
    end
    
    % initialising sums for lyapunov exponents
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    
    % iterations to calculate lyapunov exponent
    for i = 1 : numIterations
        
        %---------------------------------------------------------------
        K1 = s;
        L1 = a*s + b*p + c*p*p*p + f1*cos(z) + f2*cos(zTau);
  
        M1 = dv1;
        N1 = a*dv1 + b*dx1 + 3*c*p*p*dx1 - f1*sin(z)*dz1 -...
            f2*sin(zTau)*dzTau1;
        
        F1 = dv2;
        G1 = a*dv2 + b*dx2 + 3*c*p*p*dx2 - f1*sin(z)*dz2 -...
            f2*sin(zTau)*dzTau2;
        
        H1 = dv3;
        J1 = a*dv3 + b*dx3 + 3*c*p*p*dx3 - f1*sin(z)*dz3 -...
            f2*sin(zTau)*dzTau3;
        
        U1 = dv4;
        V1 = a*dv4 + b*dx4 + 3*c*p*p*dx4 - f1*sin(z)*dz4 -...
            f2*sin(zTau)*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K2 = s + (0.5*h*L1);
        L2 = a*(s + (0.5*h*L1)) + b*(p + (0.5*h*K1)) +...
            c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(p + (0.5*h*K1)) +...
            f1*cos((z + (0.5*h))) + f2*cos(zTau + (0.5*hTau));
        
        M2 = dv1 + (0.5*h*N1);
        N2 = a*(dv1 + (0.5*h*N1)) + b*(dx1 + (0.5*h*M1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx1 + (0.5*h*M1)) -...
            f1*sin((z + (0.5*h)))*dz1 - f2*sin(zTau + (0.5*hTau))*dzTau1;
        
        F2 = dv2 + (0.5*h*G1);
        G2 = a*(dv2 + (0.5*h*G1)) + b*(dx2 + (0.5*h*F1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx2 + (0.5*h*F1)) -...
            f1*sin((z + (0.5*h)))*dz2 - f2*sin(zTau + (0.5*hTau))*dzTau2;
        
        H2 = dv3 + (0.5*h*J1);
        J2 = a*(dv3 + (0.5*h*J1)) + b*(dx3 + (0.5*h*H1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx3 + (0.5*h*H1)) -...
            f1*sin((z + (0.5*h)))*dz3 - f2*sin(zTau + (0.5*hTau))*dzTau3;
        
        U2 = dv4 + (0.5*h*V1);
        V2 = a*(dv4 + (0.5*h*V1)) + b*(dx4 + (0.5*h*U1)) +...
            3*c*(p + (0.5*h*K1))*(p + (0.5*h*K1))*(dx4 + (0.5*h*U1)) -...
            f1*sin((z + (0.5*h)))*dz4 - f2*sin(zTau + (0.5*hTau))*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K3 = s + (0.5*h*L2);
        L3 = a*(s + (0.5*h*L2)) + b*(p + (0.5*h*K2)) +...
            c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(p + (0.5*h*K2)) +...
            f1*cos((z + (0.5*h))) + f2*cos(zTau + (0.5*hTau));
        
        M3 = dv1 + (0.5*h*N2);
        N3 = a*(dv1 + (0.5*h*N2)) + b*(dx1 + (0.5*h*M2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx1 + (0.5*h*M2)) -...
            f1*sin((z + (0.5*h)))*dz1 - f2*sin(zTau + (0.5*hTau))*dzTau1;
        
        F3 = dv2 + (0.5*h*G2);
        G3 = a*(dv2 + (0.5*h*G2)) + b*(dx2 + (0.5*h*F2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx2 + (0.5*h*F2)) -...
            f1*sin((z + (0.5*h)))*dz2 - f2*sin(zTau + (0.5*hTau))*dzTau2;
        
        H3 = dv3 + (0.5*h*J2);
        J3 = a*(dv3 + (0.5*h*J2)) + b*(dx3 + (0.5*h*H2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx3 + (0.5*h*H2)) -...
            f1*sin((z + (0.5*h)))*dz3 - f2*sin(zTau + (0.5*hTau))*dzTau3;
        
        U3 = dv4 + (0.5*h*V2);
        V3 = a*(dv4 + (0.5*h*V2)) + b*(dx4 + (0.5*h*U2)) +...
            3*c*(p + (0.5*h*K2))*(p + (0.5*h*K2))*(dx4 + (0.5*h*U2)) -...
            f1*sin((z + (0.5*h)))*dz4 - f2*sin(zTau + (0.5*hTau))*dzTau4;
        %---------------------------------------------------------------
        
        %---------------------------------------------------------------
        K4 = s + (h*L3);       
        L4 = a*(s + h*L3) + b*(p + h*K3) +...
            c*(p + h*K3)*(p + h*K3)*(p + h*K3) +...
            f1*cos((z + h)) + f2*cos(zTau + hTau);
       
        M4 = dv1 + (h*N3);
        N4 = a*(dv1 + (h*N3)) + b*(dx1 + (h*M3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx1 + (h*M3)) -...
            f1*sin((z + h))*dz1 - f2*sin(zTau + hTau)*dzTau1;
        
        F4 = dv2 + (h*G3);
        G4 = a*(dv2 + (h*G3)) + b*(dx2 + (h*F3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx2 + (h*F3)) -...
            f1*sin((z + h))*dz2 - f2*sin(zTau + hTau)*dzTau2;
        
        H4 = dv3 + (h*J3);
        J4 = a*(dv3 + (h*J3)) + b*(dx3 + (h*H3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx3 + (h*H3)) -...
            f1*sin((z + h))*dz3 - f2*sin(zTau + hTau)*dzTau3;
        
        U4 = dv4 + (h*V3);
        V4 = a*(dv4 + (h*V3)) + b*(dx4 + (h*U3)) +...
            3*c*(p + (h*K3))*(p + (h*K3))*(dx4 + (h*U3)) -...
            f1*sin((z + h))*dz4 - f2*sin(zTau + hTau)*dzTau4;
        %---------------------------------------------------------------
        
        p = p + (1/6)*(K1 + (2*K2) + (2*K3) + K4)*h;
        s = s + (1/6)*(L1 + (2*L2) + (2*L3) + L4)*h;
        
        dx1 = dx1 + (1/6)*(M1 + (2*M2) + (2*M3) + M4)*h;
        dv1 = dv1 + (1/6)*(N1 + (2*N2) + (2*N3) + N4)*h;

        dx2 = dx2 + (1/6)*(F1 + (2*F2) + (2*F3) + F4)*h;
        dv2 = dv2 + (1/6)*(G1 + (2*G2) + (2*G3) + G4)*h;

        dx3 = dx3 + (1/6)*(H1 + (2*H2) + (2*H3) + H4)*h;
        dv3 = dv3 + (1/6)*(J1 + (2*J2) + (2*J3) + J4)*h;
        
        dx4 = dx4 + (1/6)*(U1 + (2*U2) + (2*U3) + U4)*h;
        dv4 = dv4 + (1/6)*(V1 + (2*V2) + (2*V3) + V4)*h;
        
        z = z + h;
        hCount = hCount + 1;
        
        zTau = zTau + hTau;
        hTauCount = hTauCount + 1;
        
         if hCount == stepSize
            z = 0;
            hCount = 0;
        end
        
        if hTauCount == stepSize 
           zTau = mod(zTau, (2*pi));
           hTauCount = 0; 
        end
        
        % Rescaling step to prevent folding or overlfow etc.
        error1 = [dx1, dv1, dz1, dzTau1]; 
        unitError1 = error1 / sqrt(error1(1)^2 + error1(2)^2 + error1(3)^2 + error1(4)^2);
        
        error2 = [dx2, dv2, dz2, dzTau2];
        error2 = error2 - (error2(1)*unitError1(1) + error2(2)*unitError1(2) + error2(3)*unitError1(3) + error2(4)*unitError1(4))*unitError1;
        unitError2 = error2 / sqrt(error2(1)^2 + error2(2)^2 + error2(3)^2 + error2(4)^2);
        
        error3 = [dx3, dv3, dz3, dzTau3];
        error3 = error3 - (error3(1)*unitError2(1) + error3(2)*unitError2(2) + error3(3)*unitError2(3) + error3(4)*unitError2(4))*unitError2...
            - (error3(1)*unitError1(1) + error3(2)*unitError1(2) + error3(3)*unitError1(3) + error3(4)*unitError1(4))*unitError1;
        unitError3 = error3 / sqrt(error3(1)^2 + error3(2)^2 + error3(3)^2 + error3(4)^2);
        
        error4 = [dx4, dv4, dz4, dzTau4];
        error4 = error4 - (error4(1)*unitError3(1) + error4(2)*unitError3(2) + error4(3)*unitError3(3) + error4(4)*unitError3(4))*unitError3...
            - (error4(1)*unitError2(1) + error4(2)*unitError2(2) + error4(3)*unitError2(3) + error4(4)*unitError2(4))*unitError2...
            - (error4(1)*unitError1(1) + error4(2)*unitError1(2) + error4(3)*unitError1(3) + error4(4)*unitError1(4))*unitError1;
        unitError4 = error4 / sqrt(error4(1)^2 + error4(2)^2 + error4(3)^2 + error4(4)^2);
        
        newValues1 = unitError1 * initCondSize1;
        dx1 = newValues1(1); dv1 = newValues1(2); dz1 = newValues1(3); dzTau1 = newValues1(4);
        
        newValues2 = unitError2 * initCondSize2;
        dx2 = newValues2(1); dv2 = newValues2(2); dz2 = newValues2(3); dzTau2 = newValues2(4);
        
        newValues3 = unitError3 * initCondSize3;
        dx3 = newValues3(1); dv3 = newValues3(2); dz3 = newValues3(3); dzTau3 = newValues3(4);
        
        newValues4 = unitError4 * initcondSize4;
        dx4 = newValues4(1); dv4 = newValues4(2); dz4 = newValues3(3); dzTau4 = newValues4(4);
        
        % Adding to sum to calculate lyapunov exponents.
        sum1 = sum1 + log((sqrt(error1(1)^2 + error1(2)^2 + error1(3)^2 + error1(4)^2)) / initCondSize1);
        sum2 = sum2 + log((sqrt(error2(1)^2 + error2(2)^2 + error2(3)^2 + error2(4)^2)) / initCondSize2);
        sum3 = sum3 + log((sqrt(error3(1)^2 + error3(2)^2 + error3(3)^2 + error3(4)^2)) / initCondSize3);
        sum4 = sum4 + log((sqrt(error4(1)^2 + error4(2)^2 + error4(3)^2 + error4(4)^2)) / initcondSize4);
        
    end 
    
    % taking the average to get final answer for LE and representing them
    % as a vector.
    lExponent1 = sum1 / (numIterations * h);
    lExponent2 = sum2 / (numIterations * h);
    lExponent3 = sum3 / (numIterations * h);
    lExponent4 = sum4 / (numIterations * h);
    
    lExponents = [lExponent1, lExponent2, lExponent3, lExponent4];
    
end