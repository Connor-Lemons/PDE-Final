function [y_vals, x_vals, X, delta_x, delta_y] = Radiator_fin_MOL(t_iter, x_iter)    
    y_0 = 0;
    y_f = pi/16;
    y_vals = linspace(y_0, y_f, t_iter);
    beta=10000;


    L = pi/8;
    X = zeros(t_iter, x_iter);
    x_vals = linspace(0,L,x_iter);
    X(end, :) = 0;
    X(:,1) = 0;
    X(:,end) = 0;
    X(1,:)= beta*x_vals.*(L-x_vals);
    
    K = zeros(4, x_iter);
    delta_y = y_vals(2) - y_vals(1);
    delta_x = x_vals(2) - x_vals(1);

    for i = 2 : t_iter
        K(1,2:x_iter-1) = (X(i-1,3:x_iter) - 2*X(i-1,2:x_iter-1) + X(i-1,1:x_iter-2))/delta_x^2;
        K(2,2:x_iter-1) = ((X(i-1,3:x_iter) + (delta_y/2)*K(1,3:x_iter)) - 2*(X(i-1,2:x_iter-1) + (delta_y/2)*K(1,2:x_iter-1)) + (X(i-1,1:x_iter-2) + (delta_y/2)*K(1,1:x_iter-2)))/delta_x^2;
        K(3,2:x_iter-1) = ((X(i-1,3:x_iter) + (delta_y/2)*K(2,3:x_iter)) - 2*(X(i-1,2:x_iter-1) + (delta_y/2)*K(2,2:x_iter-1)) + (X(i-1,1:x_iter-2) + (delta_y/2)*K(2,1:x_iter-2)))/delta_x^2;
        K(4,2:x_iter-1) = ((X(i-1,3:x_iter) + (delta_y)*K(3,3:x_iter)) - 2*(X(i-1,2:x_iter-1) + (delta_y)*K(3,2:x_iter-1)) + (X(i-1,1:x_iter-2) + (delta_y)*K(3,1:x_iter-2)))/delta_x^2;
    
        X(i,:) = X(i-1,:) + (delta_y/6)*(K(1,:) + K(2,:) + K(3,:) + K(4,:));
        
    end

 str = "";
    decide = input(str + "Press enter to exit. Type 1 for Finite Difference Solution", "s");
    switch decide
        case "1"

            % Create meshgrid for surface plot
            [X_grid, Y_grid] = meshgrid(x_vals, y_vals(1:80:end));  % Assuming x and y are 1D vectors
            
            figure;
            hold on
            h = surf(Y_grid, X_grid, X(1:80:end,:), "FaceAlpha", 0.25, 'FaceColor', 'b');
            xlabel('y');
            ylabel('x');
            zlabel('Temperature');
            title('2D Heat MOL Solution');
            view(3); 
        case isempty(decide)
            return
    end


end
