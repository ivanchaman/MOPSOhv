%Function that calculates the delta value used in the no uniform mutation

%Parameters
%t              -       First Parameter
%y              -       Secod parameter
%generation     -       Actual Generation
%maxGen         -       Max number of generations
function res = delta(t, y, generation, maxGen)
    r = rand;
    b = 5;
    res = y * (1 - (r^((1 - (generation/maxGen))^b)));