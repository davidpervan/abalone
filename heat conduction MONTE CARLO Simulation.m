%clear all
%clc

T1=90000; T0=29500; l=0.45;

NumberOfNodes = 50; % User is able to choose the number of nodes along the wall. 
NumberOfCells = NumberOfNodes - 1; % Resulting number of cells, according to the chosen number of nodes. 
thetax = l/NumberOfCells; % Resulting distance between cells. 

T1=T1-29500; T0=T0-29500;

T_MC = zeros(1,NumberOfNodes);
T_MC(1:NumberOfNodes) = T0;

HeatSource = 1;
T_MC(HeatSource) = T1;
%T_MC = T_MC*1e+2;

%MC_runs = 1e+4;
MC_runs = 1000;
GraphsToBePlotted = 5;

Ball_Count=zeros(MC_runs,1);
Ball_Count(1)=sum(T_MC);
counter=0;
balls=0;

for i = 1:GraphsToBePlotted
    for j = 1:(MC_runs/GraphsToBePlotted)
        for x = 1:NumberOfNodes;
            balls=T_MC(x);
            balls=round(balls);
            if balls > 0;
                for b = 1:balls
                    r=rand;
                    if r < 0.5 && x > 1
                        T_MC(x)=T_MC(x)-1;
                        T_MC(x-1)=T_MC(x-1)+1;
                    elseif r >= 0.5 && x < NumberOfNodes
                        T_MC(x)=T_MC(x)-1;
                        T_MC(x+1)=T_MC(x+1)+1;
                    end
                end
            end
            T_MC(1)=T1;
            if T_MC(x) > T1
                T_MC(x) = T1;
            elseif T_MC(x) < T0
                Tn(x) = T0;
            end
        end
        counter=counter+1;
        Ball_Count(counter)=sum(T_MC);
    end
    T_MC=(T_MC/100)+295;
    figure(1);
    plot(T_MC)
    hold on   
    T_MC=(T_MC-295)*100;
    %T_MC
end
T_MC=(T_MC/100)+295;
figure(2);
plot(Ball_Count);