%Function that creates two lists of n/2 dimension with the indexes of the parents selected by
%binary torunament

%Parameters:
%pop        -   Population

%Returns:
%parents1   -   Indexes of the first parents
%parents2   -   Indexes of the second parents
function [parents1 parents2] = tournamentSelection(pop)
    %Calculate the size of the arrays
    n = floor(length(pop) / 2);
    
    %Create a random list of indexes to compare
    %Will compare the 1 - n/2 versus n/2+1 - n individuals
    participants = randperm(length(pop));
    parents1 = [];
    for i = 1 : n
        %Compares acording the rank 
        if pop(participants(i)).rank < pop(participants(n + i)).rank
            parents1(end + 1) = participants(i);
        elseif pop(participants(i)).rank > pop(participants(n + i)).rank
            parents1(end + 1) = participants(n + i);
        %If is the same rank compares acording the crowding distance
        elseif pop(participants(i)).distance > pop(participants(n + i)).distance
            parents1(end + 1) = participants(i);
        else
            parents1(end + 1) = participants(n + i);
        end
    end
    
    %Create a second random list of indexes to compare
    %Will compare the 1 - n/2 versus n/2+1 - n individuals
    participants = randperm(length(pop));
    parents2 = [];
    for i = 1 : n
        %Compares acording the rank 
        if pop(participants(i)).rank < pop(participants(n + i)).rank
            parents2(end + 1) = participants(i);
        elseif pop(participants(i)).rank > pop(participants(n + i)).rank
            parents2(end + 1) = participants(n + i);
        %If is the same rank compares acording the crowding distance
        elseif pop(participants(i)).distance > pop(participants(n + i)).distance
            parents2(end + 1) = participants(i);
        else
            parents2(end + 1) = participants(n + i);
        end
    end