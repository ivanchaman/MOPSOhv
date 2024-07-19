%Function that creates two children using the intermediate Crossover

%Parameters
%parent1        -       Individual 1
%parent2        -       Individual 2
%pc             -       Probability of Crossover

%Returns
%children1
%children2
function [children1 children2] = intermediateCrossover(parent1, parent2, pc)
    %Check if the crossover can take place
    if flip(pc)
        %Create childrens 
        children1 = createIndividual();
        children2 = createIndividual();        
        n = length(parent1.cromosome);
        %Random select a position in the cromosome
        k = randint(1,1,[0 n]);
        %Select a random number in [0,1]
        alpha = rand();
        %Copy the cromosomes from the parents up to k-position
        if(k == 0)
            children1.cromosome = zeros(1,n);
            children2.cromosome = zeros(1,n);
            children1.cromosome(1:n) = ((alpha)*parent2.cromosome)+((1 - alpha)*parent1.cromosome);
            children2.cromosome(1:n) = ((alpha)*parent1.cromosome)+((1 - alpha)*parent2.cromosome);           
        else
            %Croosover the rest of the values of the cromosome
            children1.cromosome = zeros(1,n);
            children2.cromosome = zeros(1,n);
            children1.cromosome(1:k) = parent1.cromosome(1:k);
            children2.cromosome(1:k) = parent1.cromosome(1:k);
            children1.cromosome(k + 1:n) = ((alpha)*parent2.cromosome(k+1:n))+((1 - alpha)*parent1.cromosome(k+1:n));
            children2.cromosome(k + 1:n) = ((alpha)*parent1.cromosome(k+1:n))+((1 - alpha)*parent2.cromosome(k+1:n));            
        end
    else        
        %Preserve the parents
        children1 = parent1;
        children2 = parent2;
    end