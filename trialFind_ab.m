function trialFind_ab

value_eigMin = linspace(-2,-5,3);
value_eigMax = linspace(-9,-50,3);
value_gapMax = linspace(-1.5,7,3);
value_gapMin = linspace(-8,-50,3);

paramData  = NaN(length(value_eigMin),length(value_eigMax),length(value_gapMax),length(value_gapMin));

stepFin_eigMin = length(value_eigMin);
stepFin_eigMax = length(value_eigMax);
stepFin_gapMin = length(value_gapMin);
stepFin_gapMax = length(value_gapMax);

parfor step_eigMin = 1:stepFin_eigMin
    for step_eigMax = 1:stepFin_eigMax
        for step_gapMin = 1:stepFin_gapMin
            for step_gapMax = 1:stepFin_gapMax
                if  value_eigMin(step_eigMin) > value_eigMax(1,step_eigMax)
                    if value_gapMin(1,step_gapMin) < value_gapMax(1,step_gapMax)
                        [A, B] = find_ab(value_eigMin(step_eigMin), value_eigMax(step_eigMax), value_gapMax(step_gapMax), value_gapMin(step_gapMin));
                        if sum(sum(isnan(A))) == 0
                            paramData(step_eigMin, step_eigMax, step_gapMax, step_gapMin) = 1;
                        end
                    end
                end
            end
        end
    end
end


save('paramData','-v7.3');
end
