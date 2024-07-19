%Function that implement the test function FON

%Parameters:
%x      -       values of x values

%Return:
%fx     -       Values of Aptitude
function fx = fon(x)
    fx = zeros(1,2);
    sum1 = 0;
    sum2 = 0;
    for i = 1:3
       sum1 = sum1 + (x(i) - (3^(-0.5)))^2;
       sum2 = sum2 + (x(i) + (3^(-0.5)))^2;
    end
    fx(1) = 1 - exp(-sum1);
    fx(2) = 1 - exp(-sum2);