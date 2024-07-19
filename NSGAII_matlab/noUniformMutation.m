%Function that applies the No Uniform Mutation Operator

%Parameters:
%individual
%pm             -       Probability of mutation
%generation     -       Actual generation
%maxGen         -       Maimum number of the Generations
%limVar         -       Matrix that contains the limits for the values
function ind = noUniformMutation(individual, pm, generation, maxGen, limVar)
    if  flip(pm)
        k = randint(1, 1, [1 length(individual.cromosome)]);
        if flip(0.5)
            individual.cromosome(k) = individual.cromosome(k) + delta(generation,limVar(k, 2) - individual.cromosome(k), generation, maxGen);
        else
            individual.cromosome(k) = individual.cromosome(k) - delta(generation,individual.cromosome(k) - limVar(k, 1), generation, maxGen);
        end    
    end
    ind = individual;