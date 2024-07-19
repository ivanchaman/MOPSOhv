%Function that creates the structure of an individual 

%Returns
%individual     -       Empty Individual
function individual=createIndividual()
    individual.cromosome = [];
    individual.Aptitude = [];
    individual.rank = -1;
    individual.Sp = [];
    individual.np = 0;
    individual.distance = 0;
    individual.vDistance = [];
    individual.index = 0;
    
    