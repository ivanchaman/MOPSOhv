%Function that apply the uniform mutation over a individual
%Parameters:
%individual     -   Individual to be mutated
%pm             -   Probability of Mutation
%generation     -   Not used in the algoritm, only used to standarize the mutation function
%maxGen         -   Not used in the algoritm, only used to standarize the mutation function
%limVar         -   Matrix that contains the limits of the variables

function ind = uniformMutation(individual, pm, generation, maxGen, limVar)
    %'Toss a Coin' with probability equal of pm, if true apply mutation
    if  flip(pm)
        %Select a random position over the cromosome
        k = randint(1, 1, [1 length(individual.cromosome)]);
        %Generate a new value for the k position of the cromosome that is
        %in the interval defined in limVar
        individual.cromosome(k) = unifrnd(limVar(k,1), limVar(k,2));        
    end
    ind = individual;