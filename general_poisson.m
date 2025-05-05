function [x_vals, y_vals, X_analytical, X, delta_xy] = general_poisson(x_iter, y_iter, TOL, max_iter)
    
    L = 1;
    H = 1;
    X = zeros(x_iter, y_iter, 1);
    x_0 = @(x,y) 0;
    x_vals = linspace(0,L,x_iter);
    y_vals = linspace(0,H,x_iter);
    for i = 1:length(x_vals)
        for j = 1:length(y_vals)
            X(i,j,1) = x_0(x_vals(i),y_vals(j));
        end
    end
    
    X(1,:,1) = 0;
    X(end,:,1) = 0;
    X(:,1,1) = 0;
    X(:,end,1) = sin(pi*x_vals);
    
    delta_xy = x_vals(2) - x_vals(1);

    for i = 2:max_iter
            X(1,:,i) = 0;
            X(end,:,i) = 0;
            X(:,1,i) = 0;
            X(:,end,i) = sin(pi*x_vals);
        X(2:x_iter-1,2:x_iter-1,i) = 0.25*(X(3:x_iter,2:x_iter-1,i-1) + X(1:x_iter-2,2:x_iter-1,i-1) + X(2:x_iter-1,3:x_iter,i-1) + X(2:x_iter-1,1:x_iter-2,i-1) - 0);
        if max(abs(X(:,:,i) - X(:,:,i-1)),[],"all") <= TOL
            break
        end
    end

    X_analytical = zeros(x_iter, y_iter);
    
    f = @(x,y) sinh(pi*y)/sinh(pi)*sin(pi*x);
    
    for j = 2:length(x_vals) - 1
        for k = 2:length(y_vals)-1
            X_analytical(j,k) = f(x_vals(j),y_vals(k));
        end
    end

    str = "";
    decide = input(str + "Press enter to exit. Type 1 for Finite Difference vs Analytical", "s");
    switch decide
        case "1"

            % Create meshgrid for surface plot
            [X_grid, Y_grid] = meshgrid(x_vals, y_vals);  % Assuming x and y are 1D vectors
            
            figure;
            hold on
            h = surf(X_grid, Y_grid, X(:,:,end), "FaceAlpha", 0.25, 'FaceColor', 'b');
            g = surf(X_grid, Y_grid, X_analytical(:,:), "FaceAlpha", 0.25, 'FaceColor', 'r');
            xlabel('x');
            ylabel('y');
            zlabel('Temperature');
            title('2D Heat Equation Solution');
            % Fix z-axis limits to avoid recalculating each frame
            zlim([min(X(:)), max(X(:))]);
            view(3); 
        case isempty(decide)
            return
    end

end