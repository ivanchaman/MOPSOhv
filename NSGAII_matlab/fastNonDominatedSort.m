%Function that calculates the rank of the individuals

%Parameters:
%population     -      Array that contains the individuals

%Returns:
%Fronts         -       Cells with the indexes of the individuals sorted by the i-Fronts
%pop            -       Population that contains its rank value
function [Fronts pop] = fastNonDominatedSort(population)
    if length(population) <= 0
        pop = population;
        Fronts = {};
        return;
    end
    %Create the container of the Fronts
    Fronts = {};
    %Initialize the Firts Front
    Fronts{1} = [];
    for i = 1 : length(population)
       %Initialize the values for Sp and np
       population(i).Sp = [];
       population(i).np = 0;
       p = population(i);
       for j = 1 : length(population)
          if i == j
              continue;
          end
          q = population(j);
          %if p dominates q agregate the index j to Sp
          if paretoDominance(p, q)
              population(i).Sp(end + 1) = j;
          elseif paretoDominance(q, p)
              %Increment the value for np
              population(i).np = population(i).np + 1;              
          end
       end
       %Save the index for the individuals that are in the first front
       if population(i).np == 0
          population(i).rank = 1;
          Fronts{1}(end + 1) = i;
       end
       pop(i) = population(i);
    end
    
    i = 1;
    while 1    
        Q = [];    
        %Loop trought the fron Fi
        for j = 1 : length(Fronts{i})
            %Select the individuals from Fi
            p = pop(Fronts{i}(j));     
            for k = 1 : length(p.Sp)
                %Decrement by 1 the value np, for the individuals in the
                %set Sp
                q = pop(p.Sp(k));
                np = q.np;                        
                np = np - 1;
                pop(p.Sp(k)).np = np; 
                %If the counter np is equal to 0 save into Q
                if np == 0
                   pop(p.Sp(k)).rank = i + 1;
                   Q(end + 1) = p.Sp(k);               
                end
            end
        end
        if isempty(Q)
           break; 
        end    
        i = i + 1;
        %Create the front Fi
        Fronts{i} = Q;
    end    
    F = Fronts;