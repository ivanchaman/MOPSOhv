%File that Applies the NSGA-II Algorithm
clear all;
clc;

%Parameters of the algorithm
nIndividuals = 50;      %Size of the population
nobj = 2;               %Number of objectives functions
numCrom = 2;            %Length of the cromosome
maxGen = 200;             %Number of Generation
Pm = 0.05;              %Probability of Mutation
Pc = 0.95;              %Probability of Crossover
%limVar = [-5 5; -5 5; -5 5];       %Matrix that defines the range of the variables
limVar = [-pi pi; -pi pi];       %Matrix that defines the range of the variables


fAptitude = @pol;      %Especify the function that calculates the Aptitude

%Create an empty population
pop=repmat(createIndividual(),nIndividuals,1);
for i= 1:length(pop)
    %Randomly assign the values of the cromossome
    pop(i).cromosome = unifrnd(limVar(:,1),limVar(:,2))';
    %Calculate the Aptitude
    pop(i).Aptitude = fAptitude(pop(i).cromosome);
end

generation = 0;

%Create the Q0 from the first individuals
[Fronts pop] = fastNonDominatedSort(pop);
pop = crowdingDistance(Fronts, pop);
Q = createPopulation(pop, Pc, Pm, generation, maxGen, limVar,fAptitude);

generation = 1;
for j = 1:maxGen
    
    %Create the full individuals set combining the Pi with Qi
    R = [pop Q];    
    for i=1:length(R)
       %Assign an auxiliar index to help to choose the individuals in
       %selection
       R(i).index = i; 
    end
    
    %Apply the Fast Non Dominated Sort over R
    [Fronts R] = fastNonDominatedSort(R);
    %Calculate the crowding Population on the individuals of R
    R = crowdingDistance(Fronts, R);
    
    indexPop = [];
    i = 1;
    %Save the indexes values from the individuals in the Fi front
    while (length(indexPop) + length(Fronts{i})) < nIndividuals
       indexPop = [indexPop Fronts{i}];
       i = i + 1; 
    end
    %Calculate the number of Individuals that needs for complete the
    %population Popi+1
    leftIndividuals = nIndividuals - length(indexPop);
    if(leftIndividuals > 0)
        IA = R(Fronts{i});
        vDist = [IA.distance];
        %Sort the individuals in the Fi front according to it's crowding
        %distance
        [val ind] = sort(vDist,'descend');
        %Complete the indexes for the Popi+1 population
        for k = 1: leftIndividuals
            indexPop(end + 1) = IA(ind(k)).index;        
        end
    end
    %Select the individuals according its indexes
    pop = R(indexPop);
    disp(['Generation No.' num2str(generation)]);
    %Apply the fast non dominated Sort for the new individuals, then
    %calculates its crowding distance
    [F pop] = fastNonDominatedSort(pop);
    pop = crowdingDistance(F, pop);
    Q = createPopulation(pop, Pc, Pm, generation, maxGen, limVar,fAptitude);
    %Incremenet the index of the generation
    generation = generation + 1;
end
pop = pop(F{1});
vAptitude = [pop.Aptitude];
vAptitude = reshape(vAptitude ,nobj, length(pop));
plot(vAptitude(1,:),vAptitude(2,:),'*')