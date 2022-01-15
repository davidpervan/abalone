clear all
%clc

T1=900; T0=295; Dh=4.3e-7; l=0.45; 
thetaT=T1-T0;

MaxNumberOfNodes = 200;
MinNumberOfNodes = 10;
NodesStepSize = 10;
NodesAmountOfSteps = (MaxNumberOfNodes-MinNumberOfNodes)/NodesStepSize;


RMSD = zeros(NodesAmountOfSteps,1);
Nodes = zeros(NodesAmountOfSteps,1);
r = 0;

for NumberOfNodes = MinNumberOfNodes:NodesStepSize:MaxNumberOfNodes;
    NumberOfNodes = 150;
    NumberOfCells = NumberOfNodes - 1;
    thetax = l/NumberOfCells;

    tend = 70.840e+3;
    thetat = ((l/NumberOfCells).^2)/(2*Dh); 
    NumberOfROUNDEDTimeSteps = tend/thetat;
    NumberOfROUNDEDTimeSteps = round(NumberOfROUNDEDTimeSteps); %use integer values rather than floating points for NumberOfPointInTime
    thetat = tend/NumberOfROUNDEDTimeSteps;
    NumberOfPointsInTime = NumberOfROUNDEDTimeSteps + 1;

    Tanalytical = zeros(NumberOfPointsInTime,NumberOfNodes);
    Tanalytical(1,2:NumberOfNodes) = T0;
    Tanalytical(1,1) = T1;

    i=1; 
    x=l;  
    N=1;

    for t=thetat:thetat:tend;
        i=i+1;
        j=0;
        for x=0:thetax:l;
            j=j+1;
            thetasum=0;
            for n=0:N;
                err1=erfc((2*l*n+x)/(2*sqrt(Dh*t)));
                err2=erfc((2*(n+1)*l-x)/(2*sqrt(Dh*t)));
                thetasum=thetasum+((-1.)^n)*(err1+err2);
            end
            Tanalytical(i,j) = thetasum*thetaT+T0;
        end
    end

    % Numerical part

    Tnumerical = zeros(NumberOfPointsInTime,NumberOfNodes);
    Tnumerical(1,2:NumberOfNodes) = T0;
    Tnumerical(1:NumberOfPointsInTime,1) = T1;

    for t = 2:1:NumberOfPointsInTime;
        for x = 2:1:NumberOfCells;
            %T(t,1) = T(t,1+1);
        %Tnumerical(t,1) = Tnumerical(t-1,NumberOfNodes);
        Tnumerical(t,x) = Dh*(thetat/((thetax).^2))*(Tnumerical(t-1,x-1)-2*Tnumerical(t-1,x)+Tnumerical(t-1,x+1))+Tnumerical(t-1,x);
    end
    x = NumberOfNodes;
    Tnumerical(t,NumberOfNodes) = Tnumerical(t,NumberOfNodes-1); 
    end

    % Square-Root-Mean-Difference

    %if size(Tanalytical) == size(Tnumerical) 
    Tdifference = zeros(NumberOfPointsInTime,NumberOfNodes);
    SigmaError = 0; ArithmeticMean = 0;
    SizeofTmatrices = NumberOfPointsInTime*NumberOfNodes;
    for SRT = 1:1:SizeofTmatrices
        Tdifference(SRT) = (Tanalytical(SRT)-Tnumerical(SRT)).^2;
        SigmaError = SigmaError + Tdifference(SRT);
    end
    ArithmeticMean = SigmaError/SizeofTmatrices;
    RootMeanSquareDifference = sqrt(ArithmeticMean);
    r = r + 1;
    RMSD(r) = RootMeanSquareDifference;
    Nodes(r) = NumberOfNodes;
end