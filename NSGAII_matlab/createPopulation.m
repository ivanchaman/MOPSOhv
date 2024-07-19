%Function that generates a new population, applying selection, crossover
%and selection

%Parameters:
%population     -       Parents population
%pc             -       Probabily of Crossover
%pm             -       Probability of Mutation
%generation     -       Actual Generation (Used just if needed in crossover and/or mutation)
%genMax         -       Max generation(Used just if needed in crossover and/or mutation)
%limVar         -       Matrix with the limits of variables
%cost           -       Function to calculate the cost

%Returns:
%pop            -       New Population of n individuals
function pop = createPopulation(population, pc, pm, generation, genMax, limVar, cost)
    %Assign the Croosover and mutation that is used after
    crossover = @intermediateCrossover;
    mutation = @uniformMutation;    
    n = length(population);
    %Creates a list of the parents
    [parentsA parentsB] = tournamentSelection(population); 
    A = population(parentsA);
    B = population(parentsB);
    
    %Create an array of n empty individuals
    pop=repmat(createIndividual(),1,n);    
    for i = 1:floor(n/2)
        %Apply the crossover to create the new childrens
        [children1 children2] = crossover(population(parentsA(i)), population(parentsB(i)), pc);
        
        %Apply the mutation operator in the childrens
        children1 = mutation(children1, pm, generation, genMax, limVar);
        children2 = mutation(children2, pm, generation, genMax, limVar);
        
        %Calculate the Aptitude of the childrens
        children1.Aptitude = cost(children1.cromosome);
        children2.Aptitude = cost(children2.cromosome);
        
        %Save the new individuals
        pop(2*i - 1) = children1;
        pop(2*i) = children2;
    end        