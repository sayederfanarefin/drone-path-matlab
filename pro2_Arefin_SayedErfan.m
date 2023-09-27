function pro2_Arefin_SayedErfan(mobilityModel, velocity, pausingTime, numberOfDummy)

    height = 1000;
    width = 1000;
    timesToRepeat = 5;
    startX=0;
    startY=0;

    rectplot = @(x1,x2) rectangle('Position', [0, 0, width, height], 'EdgeColor', 'k', 'LineWidth', 2);

    totalDistance = 0;
    totalTime = 0;
    figure;
    

    if mobilityModel == 0
        %rwp
        oldX=0;
        oldY=0;
        for i = 1:timesToRepeat
            srcX=0;
            srcY=0;

            dstX=0;
            dstY=0;

            if i == 1
                [randomX, randomY]  = randomPoint(startX,startY, height, width);
                dstX=randomX;
                dstY=randomY;
                srcX=startX;
                srcY=startY;
                oldX=randomX;
                oldY=randomY;
            else
                srcX=oldX;
                srcY=oldY;
                [randomX, randomY]  = randomPoint(startX,startY, height, width);
                dstX=randomX;
                dstY=randomY;
                oldX = randomX;
                oldY = randomY;
            end

            % disp([i, srcX,srcY,dstX,dstY, oldX, oldY]);
            [dis, timeTaken] = drawLine(srcX,srcY,dstX,dstY);
            totalDistance = dis + totalDistance;
            totalTime = totalTime + timeTaken + pausingTime;
            
            % pause(pausingTime)
            hold on;

            

        end
        
        rectplot(height,width);
        title('Drone path using RWP');
        xlabel('X');
        ylabel('Y');
        axis equal;
        grid on;
        hold off;

        fprintf('Total Flying Distance: %.2f meters\n', totalDistance);
        fprintf('Total Flying Time: %.2f seconds\n\n', totalTime);

    else
        % ppr


        counter = 1;
        srcX=startX;
        srcY=startY;

        for i = 1:timesToRepeat
            
            [randomX, randomY]  = randomPoint(srcX,srcY, height, width);
            dstX=0;
            dstY=0;
            

            [x1, y1, x2, y2, x3, y3, x4, y4] = getRectanglePoints(srcX, srcY, randomX, randomY);

            dummyLocations = zeros(numberOfDummy, 2);
            % disp("Calculating dummies");
            for ii = 1:numberOfDummy

                [dumX, dumY] = randomPoint(x1, y1, x4, y3);
                dummyLocations(ii, :) = [dumX, dumY];

                % disp ([dumX, dumY]);
            end

            oldX=0;
            oldY=0;

            for ii = 1:numberOfDummy
                if counter ==1
                    dstX=dummyLocations(ii,1);
                    dstY=dummyLocations(ii,2);
                    oldX=dstX;
                    oldY=dstY;

                else
                    dstX=dummyLocations(ii,1);
                    dstY=dummyLocations(ii,2);
                    srcX=oldX;
                    srcY=oldY;
                    oldX=dstX;
                    oldY=dstY;
                end
                % disp([srcX, srcY, dstX, dstY]);
                [dis, timeTaken] = drawLine(srcX,srcY,dstX,dstY);
                totalDistance = dis + totalDistance;
                totalTime = totalTime + timeTaken + pausingTime;
                
                % pause(pausingTime)
                counter = counter + 1;
                
            end
            srcX = dstX;
            srcY=dstY;



                % disp("----");
                % disp(dummyLocations());
                % disp(dummyLocations(1,2));
                % oldX=randomX;
                % oldY=randomY;
           

        end
        rectplot(height,width);
        title('Drone path using PPR');
        xlabel('X');
        ylabel('Y');
        axis equal;
        grid on;
        hold off;

        fprintf('Total Flying Distance: %.2f meters\n', totalDistance);
        fprintf('Total Flying Time: %.2f seconds\n\n', totalTime);


    end






    


    

    function [lineDistance, timeTaken] = drawLine(x1, y1, x2, y2)

        lineDistance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        timeTaken = lineDistance / velocity;

        t = 0:0.1:timeTaken;
        x = x1 + (x2 - x1) * t / timeTaken;
        y = y1 + (y2 - y1) * t / timeTaken;

        % disp(length(t));
        % disp(length(x));
        % disp(length(y));
        
        for j = 1:length(t)
            
            plot(x(j), y(j), '-*', 'MarkerSize', 4 , 'LineWidth', 1, 'Color', 'r');
            hold on;

            plot(x(end), y(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
            hold on;

            % pause(t(j));
            % totalTime = totalTime + t(j);
        end
        



        

        % Plot the source and destination points
        % plot(x1, y1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        % plot(x2, y2, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
        % Draw the straight line between the source and destination points
        % plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 2);

        

    end








    function [randomX, randomY] = randomPoint(xmin, ymin, xmax, ymax)
        randomValue1 = rand;
        randomValue2 = rand;
        
        randomX = ceil(xmin + randomValue1 * (xmax - xmin));
        randomY = ceil(ymin + randomValue2 * (ymax - ymin));
      
    end


    function [x1, y1, x2, y2, x3, y3, x4, y4] = getRectanglePoints(x1, y1, x2, y2)

        minX = min(x1, x2);
        maxX = max(x1, x2);
        minY = min(y1, y2);
        maxY = max(y1, y2);
    

        x3 = minX;
        y3 = maxY;
        x4 = maxX;
        y4 = minY;

        
    end
end
