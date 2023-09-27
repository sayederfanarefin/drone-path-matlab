function[]=pro3_SayedErfan_Arefin(typeOfRun)


startPos = [0,0];
[pois] = scale_free();
r = 100;

[num_rows, num_cols] = size(pois);
poisTrack = zeros(num_rows, 2); % Initialize the matrix with zeros, col1 for count, col2 for visited or not
for row = 1:num_rows
    
    poisTrack(row, 1) = poisTrack(row, 1) + 1;
    for row2 = 1:num_rows
        
        if poisTrack(row2, 2)==0
           % meaning not included in any cluster
            x1 = pois(row, 1);
            y1 = pois(row, 2);
    
            x2 = pois(row2, 1);
            y2 = pois(row2, 2);
    
            distance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
            if distance <= r
                poisTrack(row, 1) = poisTrack(row, 1) + 1;
                poisTrack(row2, 2) = 1;
                %display(distance)
            end
        end
    end

    % disp(poisTrack(row, 3));
end
poisFull = [pois, poisTrack];
sorted_poisFull = sortrows(poisFull, -3);
clusters = sorted_poisFull(sorted_poisFull(:, 3) > 1, :);
clusters = clusters(:, 1:3);


% disp(sorted_poisFull);
% disp(clusters);
clusters_size = size(clusters);
disp("------clusters_size-----");
disp(clusters_size);

for row = 1:clusters_size 
    scatter(clusters(row, 1), clusters(row, 2), 20, 'b', '^', 'filled');
    rectangle('Position', [clusters(row, 1) - r, clusters(row, 2) - r, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
    hold on;
end
hold off;

grid on;
set(gca, 'FontSize', 8);
set(gca, 'xTick', [-200:200:1200]);
set(gca, 'yTick', [-200:200:1200]);
xlabel('X', 'FontSize', 8);
ylabel('Y', 'FontSize', 8); 
axis([-200 1200 -200 1200]);
print -depsc sf_topo;
savefig('clusters.fig');

figure;
for row = 1:clusters_size 
    scatter(clusters(row, 1), clusters(row, 2), 30, 'b', '^', 'filled');
    
    hold on;
end


scatter(startPos(1, 1), startPos(1, 2), 80, 'b', '^', 'filled');
hold on;

totalDist = 0;
    if typeOfRun == 0
        
        
        % Clustering and RAND 
        
        for rowx = 1:clusters_size 
            if rowx == 1
                result = drawLine(startPos(1,1), startPos(1,2), clusters(rowx, 1), clusters(rowx, 2));
                totalDist = totalDist  + result;
            
            else
               result =  drawLine(clusters(rowx-1, 1), clusters(rowx-1, 2), clusters(rowx, 1), clusters(rowx, 2));
                totalDist = totalDist  + result;
            end
        end
        hold off;


    elseif typeOfRun == 1

        %  Clustering and NNF
        clustersJustPoints = clusters(:, 1:2);

        tracker = zeros(clusters_size(1), 2);
        [nearest_neighbor, nearest_neighbor_index] = findNearestNeighbor(clustersJustPoints, startPos);

       
        result =  drawLine(startPos(1,1), startPos(1,2), nearest_neighbor(1, 1), nearest_neighbor(1, 2));
        totalDist = totalDist  + result;
        new_matrix = clustersJustPoints;
        new_matrix(nearest_neighbor_index, :) = [];
        backUp = nearest_neighbor;
        for rowx = 1:clusters_size 
            [nearest_neighbor, nearest_neighbor_index] = findNearestNeighbor(new_matrix, backUp);

       
            result =   drawLine(backUp(1,1), backUp(1,2), nearest_neighbor(1, 1), nearest_neighbor(1, 2));
            totalDist = totalDist  + result;
            backUp = nearest_neighbor;
            new_matrix(nearest_neighbor_index, :) = [];
            
        end
        
    else
        % Clustering and DF 
        disp(clusters);

        for rowx = 1:clusters_size 
            if rowx == 1
                drawLine(startPos(1,1), startPos(1,2), clusters(rowx, 1), clusters(rowx, 2));
            else
                drawLine(clusters(rowx-1, 1), clusters(rowx-1, 2), clusters(rowx, 1), clusters(rowx, 2));
            
            end
        end

    end
    disp(totalDist );
end










function [nearest_neighbor, nearest_neighbor_index] = findNearestNeighbor(data_matrix, query_point)

    num_data_points = size(data_matrix, 1);

    % Initialize variables to store the minimum distance and index of the nearest neighbor
    min_distance = Inf;
    nearest_neighbor_index = 0;

    % Iterate over all data points
    for i = 1:num_data_points
        % Compute the Euclidean distance between the current data point and the query point
        distance = norm(data_matrix(i, :) - query_point);

        % Check if the current distance is smaller than the minimum distance found so far
        if distance < min_distance
            min_distance = distance;
            nearest_neighbor_index = i;
        end
    end

    nearest_neighbor = data_matrix(nearest_neighbor_index, :);
end





function [lineDistance] = drawLine(x1, y1, x2, y2)
        
    lineDistance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    plot([x1, x2], [y1, y2], 'r', 'LineWidth', 1);
    hold on;

end



function[loc]=scale_free()
 
    w_begin= 0;
    w_end = 1000;
    h_begin = 0;
    h_end = 1000;
    
    % number of cells (subareas): 5-by-5, by default
    n_cell = 5;
    tot_cell = n_cell * n_cell;
    size_cell = w_end / n_cell;
    
    % n number of target points: 
    n = 100;
    
    % a set of rectangular subareas: 5-by-5
    x_ = linspace(w_begin, w_end - size_cell, n_cell);
    ux = [];
    for i = 1:n_cell
        ux = [ux, x_]; 
    end 
    ux = ux'
    
    y_ = ones(1, n_cell);  
    uy = [];
    for i = 1:n_cell
        uy = [uy, y_ .* (size_cell * (i - 1))];
    end 
    uy = uy'
    
    % n number of weights: w, n-by-1, uniform
    uw = ones(n, 1);
    
    % n number of weights: w, n-by-1, uniform
    % -- between the interval (w_begin, w_end) 
    w_begin = 0;
    w_end = 10;
    w = w_begin + (w_end - w_begin) .* rand(n, 1);
    
    % coverage area with radius, r (m), by default 100
    r = 100;
    
    % -----------------------
    % scale-free distribution 
    % -----------------------
    
    % clustering exponent, alpha
    alpha = 1.4;
    
    % population, pop  
    % -- initialize to zero
    pop = ones(tot_cell, 1) - 1;
    
    % probability, prob
    % -- initialize to zero
    prob = ones(tot_cell, 1) - 1;
    
    % a set of rectangular subareas, 25-by-5
    subarea_loc = [ux, uy]
    
    % the first target point is randomly assigned to one of cells
    pos_subarea = randi(tot_cell)
    
    pos_x = randi(size_cell) + ux(pos_subarea)
    pos_y = randi(size_cell) + uy(pos_subarea)
    pop(pos_subarea) = pop(pos_subarea) + 1
    
    % the first target point - randomly assigned
    loc(1, 1) = pos_x
    loc(1, 2) = pos_y
    
    % generate all scale-free target points (x, y)
    for i = 2:n
        % calculate probabilities
        % -- sigma_pop = sum(pop, "all")
        sigma_pop = 0;
        for j = 1: tot_cell
            sigma_pop = sigma_pop + power(pop(j) + 1, alpha);
        end
        for j = 1: tot_cell
            prob(j) = power(pop(j) + 1, alpha) / sigma_pop; %power(sigma_pop, alpha);
            %prob(j) = power(pop(j), alpha) / power(sigma_pop, alpha)
        end
        % sanity check: if total probabilities are one
        %tot_prob = sum(prob, "all")
    
        % randomly choose one of subareas
        % -- pos_subarea = randi(tot_cell);
        
        % choose one of subareas based on the probability
        % -- generate a random and compare with cumulative probabilities 
        rand_prob = rand(1, 1); % generate between 0 to 1
        cumu_prob = 0; 
        for j = 1: tot_cell
            cumu_prob = cumu_prob + prob(j);
            if (cumu_prob >= rand_prob)
                pos_subarea = j;
                break
            end
        end
    
        % generate a position within the chosen subarea
        pos_x = randi(size_cell) + ux(pos_subarea);
        pos_y = randi(size_cell) + uy(pos_subarea);
        % increment the population of subarea
        pop(pos_subarea) = pop(pos_subarea) + 1;
    
        % add a new target point's (x, y) into a row
        loc = [loc; [pos_x, pos_y]];
    end    
    
    % draw target points
    plot(loc(:, 1), loc(:, 2), "rx")
    hold on
    % 
    
end