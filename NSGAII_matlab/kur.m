%Algorithm that aplyes the test function KUR with n = 3

%Parameters:
%x      -       Array with the values for x
%aptitude   -       Array with the two values for the aptitude
function  aptitude = kur(x)    
    aptitude = zeros(1,2);
    n = length(x);
    f1 = 0;
    for i=1:n-1
       f1 = f1 + (-10 * exp(-0.2 * sqrt(x(i)^2 + x(i+1)^2)));
    end
    aptitude(1) = f1;
    f2 = 0;
    for i=1:n
       f2 = f2 + ((abs(x(i))^0.8) + (5 * sin(x(i)^3)));
    end
    aptitude(2) = f2;