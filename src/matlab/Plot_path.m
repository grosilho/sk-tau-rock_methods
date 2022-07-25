clc;
clear;

problem_num = 2;
filename = 'sol';

path = './';
problem_names = {'ReversibleIsomerization', 'NonlinearReversibleReaction', 'GeneticPositiveFeedbackLoop', 'MichaelisMenten', 'SchloglReaction','DecayingDimerizing','Ecoli'};
N = [2,1,5,4,1,3,28]; % number of species, per problem

fileID = fopen([path problem_names{problem_num} '/' filename '.bin']);
X = fread(fileID,'double')';
fclose(fileID);

N = N(problem_num);
X = reshape(X,[N+1,numel(X)/(N+1)])';
t = X(:,1);
X = X(:,2:end);

% figure;
for i=1:N
    %subplot(N,1,i);
    figure;
    plot(t,X(:,i));
    hold on;
end

%T = table(t,X(:,1),X(:,2),X(:,3),X(:,4),'variablenames',{'t','x1','x2','x3','x4'});
%writetable(T,[problem_names{problem_num} '_path.csv']);
