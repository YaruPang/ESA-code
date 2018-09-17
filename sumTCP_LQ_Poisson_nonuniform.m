function [] = sumTCP_LQ_Poisson_nonuniform()

clc

format long %data is represented in a long precision.

% for k = 7:1:8

% read DVH data
data_file_name1 = '18297-D-%_PTV.txt';
data_file_name2 = '18359-D-%_PTV.txt';
%     data_file_name = sprintf('1829%d-C_PTV.txt', k);

fileID1 = fopen(data_file_name1);
data = textscan(fileID1, '%f %f %f');
D1 = data{2}/100;
V1 = data{3};
V1i=V1/sum(V1);

fileID2 = fopen(data_file_name2);
data = textscan(fileID2, '%f %f %f');
D2 = data{2}/100;
V2 = data{3};
V2i=V2/sum(V2);
%     delta_V = zeros(201,1);
%     for n = 1 : 1 : size(V, 1) - 1
%         delta_V(n, 1) = V(n, 1) - V(n+1, 1);
%     end


%     plot(D, v, 'r')

gamma50 = 3;
N = log(2) * exp(gamma50 * 2/log(2));
alpha = 0.59;
beta = 0.071;
n = 90;
D50 = 60;



TCP1 = zeros(201, 1);
TCP2 = zeros(201, 1);
size(V2)
size(D2)
for i = 2:1:200
    
    TCP1(i, 1) = TCP1(i-1, 1)+V1i(i, 1)*exp(-N * exp(-(alpha + beta * D1(i, 1)/n)*D1(i, 1)));
    TCP2(i, 1) = TCP2(i-1, 1)+V2i(i, 1)*exp(-N * exp(-(alpha + beta * D2(i, 1)/n)*D2(i, 1)));
    
end

%     TCP_LQ_Poisson_nonuniform = sum(TCP);

plot (D1(1:199,:), TCP1(2:200, :), 'b')
hold on
plot (D2(1:199,:), TCP2(2:200, :), 'b')

%     TCP_LQ_Poisson_nonuniform
%     EUD = D50/((1/TCP_LQ_Poisson_nonuniform - 1)^(1/(4*gamma50)))
%     C=-log(-1/N*log(TCP_LQ_Poisson_nonuniform));
%     B=alpha;
%     A=beta/30;
%
%     syms Deff
%     eqn = A*x^2 + B*x - C == 0;
%     Deff = solve(eqn)

end





