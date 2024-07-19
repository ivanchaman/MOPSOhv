%Function that creates a two childrens using the SBX Crossover

%Parameters:
%parent1        -       Individual 1
%parent2        -       Individual 2
%pc             -       Probability of crossover

%Returns:
%children1
%children2

function [children1 children2] = sbxCrossover(parent1, parent2, pc)
    %'Flip a coin' using the probability of crossover and check if its true
    if flip(pc)
        %Create empty individuals
        children1 = createIndividual();
        children2 = createIndividual();        
        n = length(parent1.cromosome);
        %Create a random value in [0,1]
        u = rand();
        nc = 20;
        %Calculates the value of b
        if u <= 0.5
           b = (2 * u)^(1 / (nc + 1));
        else
           b = (1 / (2 - (2 * u)))^(1 / (nc + 1)); 
        end
        %Create and ampty array of cromosome for the childrens
        children1.cromosome = zeros(1,n);
        children2.cromosome = zeros(1,n);
        %Crossover according SBX
        children1.cromosome(1:n) = 0.5 *((parent1.cromosome + parent2.cromosome) - b*abs(parent2.cromosome - parent1.cromosome));
        children2.cromosome(1:n) = 0.5 *((parent1.cromosome + parent2.cromosome) + b*abs(parent2.cromosome - parent1.cromosome));
    else
        %The parents are preserved
        children1 = parent1;
        children2 = parent2;
    end