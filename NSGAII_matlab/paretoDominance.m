%Function that compares if A is Pareto Dominate by B

%Parameters:
%a          -       Individual A
%b          -       Individual B

%Returns:
%res        -       Result of the comparision
function res = paretoDominance(a , b)
    res = 0;
    %Check if all the Aptitudes are less or equal and at least some
    %aptitude lesser
    if all(a.Aptitude <= b.Aptitude) && any(a.Aptitude < b.Aptitude)
        res = 1;
    end