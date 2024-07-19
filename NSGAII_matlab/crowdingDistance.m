%Function that calculates the crowding distance int a poblation

%Parameters:
%Fronts         -       Array that returns the fast nondominated sort function
%population     -       Array of Individuals

%Returns
%pop            _       Array of Individuals with the value of the distance
function pop = crowdingDistance(Fronts, population)
    %Calculate the number of Objective Functions
    m = length(population(1).Aptitude);
    
    for k = 1:length(Fronts)
       %Select the individuals in the Front Fi
       I = population(Fronts{k});
       %Calculate the length of the Front Fi
       l = length(I);
       %Initialize the values that gonna be used
       for j = 1 : length(I)
          I(j).distance = 0;
          I(j).vDistance = zeros(1,m);
       end
       %Obtain the list of index sorting by each objective
       listSort = sortList(I);
       for j = 1:m
            %Assign to the individuals that have the minimum and maximum
            %aptitude in the objective j, the value inf to its distance
            I(listSort(j,1)).vDistance(j) = inf;
            I(listSort(j,end)).vDistance(j) = inf;
            %Calculate the distance for the other individuals
            if (l - 2) >= 1
                maxV = I(listSort(j,end)).Aptitude(j);
                minV = I(listSort(j,1)).Aptitude(j);
                %Calculate distance
                for i = 2 : l - 1
                    I(listSort(j,i)).vDistance(j) = (I(listSort(j,i + 1)).Aptitude(j) - I(listSort(j,i - 1)).Aptitude(j))/(maxV - minV);                    
                end
            end
       end       
       for i = 1 : length(I)          
          I(i).distance = sum(I(i).vDistance);       
       end
       population(Fronts{k}) = I;       
    end    
    pop = population;