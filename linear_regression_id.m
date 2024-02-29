load ('proj_fit_16.mat')

% identification input data
x1_id = id.X{1};
x2_id = id.X{2};

% identification output data
y_id = id.Y;

% validation input data
x1_val = val.X{1};
x2_val = val.X{2};

% validation output data
y_val = val.Y;

figure
mesh(x1_id, x2_id, y_id);
title('Identification Data','FontSize',16);

figure
mesh(x1_val, x2_val, y_val);
title('Validation Data','FontSize',16);

% vector containig the number of terms for each polynomial degree
nr_terms_poly = ones(100,1);
nr_terms_poly(1) = 3;

increment = 3;
for i = 2:100
    nr_terms_poly(i) = nr_terms_poly(i-1) + increment;
    increment = increment + 1;
end

% the degree to which we check
length_check = 40;
t_MSE = 1:length_check;

MSE_id_list = zeros(length_check,1);
MSE_val_list = zeros(length_check,1);

for i = 1:length_check
   [x,y] = LR_Function(i,nr_terms_poly);
   MSE_id_list(i) = x;
   MSE_val_list(i) = y;
end

figure
plot(t_MSE, MSE_id_list,'b');
title('Mean Squared Error','FontSize',16);
grid; hold
plot(t_MSE, MSE_val_list,'r');
legend('for identification','for validation')

% taking the optimum degree
min_MSE = (min(MSE_val_list));
min_MSE_index = 0;

% finding the index of the optimum degree
for i = 1:length_check
    if MSE_val_list(i) == min_MSE
        min_MSE_index = i;
    end
end

% displaying the optimum degree's index
fprintf("Optimum degree's index: %.f\n",min_MSE_index);

% displaying the optimum degree's value
fprintf("Optimum degree's value: %f\n",min_MSE);

% finding the optimum degree's value for identification
for i = 1:length_check
    if min_MSE_index == 5
        min_MSE_id = MSE_id_list(i);
    end
end

% displaying the optimum degree's value for identification
fprintf("Optimum degree's value for identification: %f\n",min_MSE_id);

% function for finding the approximation 
function [MSE_id,MSE_val] = LR_Function(poly_degree,nr_terms_poly)
load ('proj_fit_16.mat')

% identification input data
x1_id = id.X{1};
x2_id = id.X{2};

% identification output data
y_id = id.Y();

% validation input data
x1_val = val.X{1};
x2_val = val.X{2};

% validation output data
y_val = val.Y();

N = length(y_id);

y_id_vector = ones(N*N,1);

ind = 1;
for i = 1:N
    for j = 1:N
        y_id_vector(ind) = y_id(i,j);
        ind = ind + 1;       
    end 
end

phi = zeros(N*N,nr_terms_poly(poly_degree));

% vector containing the regressors
phi_line = ones(nr_terms_poly(poly_degree), 1);

ind2 = 1;
for i = 1:N
    for j = 1:N
        phi_line_index = 1;
        phi_line_value = 0;

        for u = 0:poly_degree
            for v = 0:poly_degree
                if (u+v) <= poly_degree
                    phi_line_value = phi_line_value + ( (x1_id(i)^u) * (x2_id(j)^v) );
                    phi_line(phi_line_index) = phi_line_value;
                    phi_line_index = phi_line_index + 1;
                    phi_line_value = 0;
                end
            end
        end
     
        % adding each line to the matrix
        for f = 1:nr_terms_poly(poly_degree)
            phi(ind2,f) = phi_line(f);
        end
        
        ind2 =  ind2 + 1;
    end
end

% computing theta using liniar regression
theta = phi \ y_id_vector;

M = length(y_val);

y_hat = ones(M);

% computing the approximation using the obtained theta and validation data
for i = 1:M
    for j = 1:M
        y_hat_value = 0;
        theta_index = 1;

        for u = 0:poly_degree
            for v = 0:poly_degree
                if (u+v) <= poly_degree
                    y_hat_value = y_hat_value + (theta(theta_index) * (x1_val(i)^u) * (x2_val(j)^v) );
                    theta_index = theta_index + 1;
                end
            end
        end

        y_hat(i,j) = y_hat_value;
    end   
end

y_hat_id = ones(N);

for i = 1:N
    for j = 1:N
        y_hat_id_value = 0;
        theta_index = 1;

        for u = 0:poly_degree
            for v = 0:poly_degree
                if (u+v) <= poly_degree
                    y_hat_id_value = y_hat_id_value + (theta(theta_index) * (x1_id(i)^u) * (x2_id(j)^v) );
                    theta_index = theta_index + 1;
                end
            end
        end

        y_hat_id(i,j) = y_hat_id_value;
    end   
end

MSE_id = 0;

for i = 1:41
    for j = 1:41
        MSE_id = MSE_id + ((y_id(i,j)-y_hat_id(i,j))^2);
    end
end

MSE_id = ((1/(N*N)) * MSE_id);

MSE_val = 0;

for i = 1:M
    for j = 1:M
        MSE_val = MSE_val + ((y_val(i,j)-y_hat(i,j))^2);
    end
end

MSE_val = ((1/(M*M)) * MSE_val);

if poly_degree == 5

    figure
    subplot(1,2,1)
    mesh(x1_val, x2_val, y_val);
    title('Validation Data','FontSize',16);

    subplot(1,2,2)
    mesh(x1_val,x2_val,y_hat)
    title(['Polynomial degree: ',num2str(poly_degree)],'FontSize',16);

    figure
    subplot(1,2,1)
    mesh(x1_id, x2_id, y_id);
    title('Identification Data','FontSize',16);

    subplot(1,2,2)
    mesh(x1_id,x2_id,y_hat_id)
    title(['Polynomial degree: ',num2str(poly_degree)],'FontSize',16);
end
end
