%Function that simulates the toss of a coin with a certain probability
%Parameters:
%prob       -       Probability of calls true

%Returns:
%res        -       Result of the toss
function res = flip(prob)
    %Generates a random value and if its lesser or equal that the
    %probability it calls true
    if rand(1) <= prob
        res = 1;
    else
        res = 0;
    end