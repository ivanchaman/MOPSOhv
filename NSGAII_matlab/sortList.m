%Function that return a m-matrix with the indexes of the individuals sorted
%acccording the m objective functions

%Parameters:
%population     -       Array containing the individuals

%Returns:
%list           -       List of the index sorted for the individuals according the m
%objectives

function list = sortList(population)
    m = length(population(1).Aptitude);
    %Obtain the m Aptitude values
    vAptitude = [population.Aptitude];
    %Reshape to obtain m arrays of aptitude
    vAptitude = reshape(vAptitude ,m, length(population))';    
    list = [;];
    for i = 1 : m        
        [val index] = sort(vAptitude(:,i),'ascend');
        list(end +1,:) = index;
    end
    