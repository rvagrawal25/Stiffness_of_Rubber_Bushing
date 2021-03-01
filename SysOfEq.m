%% Obtaining the k expression in symbolic form

syms A11 A12 A13 A14 mu kappa X a b delta F ;

bc1 = A11*a*a*(kappa-2) + (A13/2)*(-1+(kappa+1)*log(a)) + (A12/2)*(-1+(kappa-1)*log(a)) + A14/(a*a) + 2*mu*X == delta*2*mu;

bc2 = A11*a*a*(kappa+2) + (A13/2)*(-1-(kappa+1)*log(a)) + (A12/2)*(-1-(kappa-1)*log(a)) + A14/(a*a) - 2*mu*X == -delta*2*mu ;

bc3 = A11*b*b*(kappa-2) + (A13/2)*(-1+(kappa+1)*log(b)) + (A12/2)*(-1+(kappa-1)*log(b)) + A14/(b*b) + 2*mu*X == 0 ;

bc4 = A11*b*b*(kappa+2) + (A13/2)*(-1-(kappa+1)*log(b)) + (A12/2)*(-1-(kappa-1)*log(b)) + A14/(b*b) - 2*mu*X == 0 ;

cd1 = A13*(kappa-1)+A12*(kappa+1) == 0;

sol = solve([bc1,bc2,bc3,bc4,cd1],[A11, A13, A12, A14, X]);

A11sol = simplify(sol.A11);
A13sol = simplify(sol.A13);
A12sol  = simplify(sol.A12);
A14sol = simplify(sol.A14);
Xsol = simplify(sol.X);

F = -2*pi*A13sol;
F = simplify(F);

k = F/(mu*delta);
k = simplify(k);


%% Plot for (k) vs (a/b)

% Let k1 be dimensionless stiffness for plane stress and let k2 be dimensionless stiffness for plane strain.
% Let the ratio of radii [a/b = r]. 
% So now we plot the graph for stiffness vs radius ratio {k1,k2 v/s r }
% For making the process simpler, we represent the expression for stiffness k as a
% function of r, instead of radii a,b.

nu = 0.5;               % We assume the rubber material to be incompressible.
K1= (3-nu)/(1+nu);      % Kolosov's Constant for plane stress.
K2= (3-4*nu);           % Kolosov's constant for plane strain.


r1 = 0:0.0001:1;         % We find the values of k1 for different values of r, at an interval of 0.01
for index = 1:10001
    k1(index) = (2*K1*(K1+1)*pi*(r1(index)*r1(index)+1))/(-K1*K1*log(r1(index))*(r1(index)*r1(index)+1)+(r1(index)*r1(index)-1));
end

r2 = 0:0.0001:1;         % We find the values of k2 for different values of r, at an interval of 0.01
for index = 1:10001
    k2(index) = (2*K2*(K2+1)*pi*(r2(index)*r2(index)+1))/(-K2*K2*log(r2(index))*(r2(index)*r2(index)+1)+(r2(index)*r2(index)-1)) ;
end


plot(r1,k1,r2,k2)
axis([0 1 0 100])
legend('Plane Stress','Plane Strain','Location','northwest')
title('Dimensionless Stiffness(F/𝜇𝛿) Vs (a/b) plot')
xlabel('a/b')
ylabel('Dimensionless Stiffness(F/𝜇𝛿)')