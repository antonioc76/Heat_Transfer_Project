clear;
radius_total = 1.5 * 0.0254;
delta_x = 0.1 * 0.0254;
total_time = 1 * 60 * 60;
delta_t = 0.01;

% fluid properties
radius_fluid = 1 * 0.0254;
fluid_initial_temperature = 96;
k_fluid = 0.571;
rho_fluid = 1000;
C_p_fluid = 4205;

% insulation layer 1 properties
radius_insulation = 0.1 * 0.0254;

insulation_initial_temperature = 23; % metal
k_insulation = 15.6; % metal
rho_insulation = 7913; % metal
C_p_insulation = 456; % metal
epsilon_insulation = 0.075; % metal

% vacuum properties
radius_vacuum = 0.3 * 0.0254;
vacuum_temperature = 23;

% insulation layer 2 properties
radius_insulation_2 = 0 * 0.0254;

insulation_2_initial_temperature = 10; %glass
k_insulation_2 = 1; %glass
rho_insulation_2 = 2500; %glass
C_p_insulation_2 = 840; %glass
epsilon_insulation_2 = 0.89; %glass

% metal properties
radius_metal = 0.1 * 0.0254;
metal_initial_temperature = 23;
k_metal = 15.6;
rho_metal = 7913;
C_p_metal = 456;
epsilon_metal = 0.075;

% ambient properties
ambient_temperature = -5;
h_air = 24;
rho_air = 1.184;
C_p_air = 1007;

% bottle dimensions
height = 6.64 * 0.0254;
diameter = 3 * 0.0254;
A_surface = 2 * pi * diameter/2 * height + 2 * pi * (diameter/2)^2;

% compute the number of nodes
number_distance_nodes = int16(radius_total / delta_x + 1);
number_time_nodes = floor(total_time / delta_t + 1);
disp(number_time_nodes)
disp(number_time_nodes * delta_t)

% compute the rounded down number of each node in the temperature matrix
number_distance_nodes = double(number_distance_nodes);
proportion_fluid_nodes = radius_fluid / radius_total;
number_fluid_nodes = floor(proportion_fluid_nodes * number_distance_nodes);
proportion_insulation_nodes = radius_insulation / radius_total;
number_insulation_nodes = floor(proportion_insulation_nodes * number_distance_nodes);
proportion_vacuum_nodes = radius_vacuum / radius_total;
number_vacuum_nodes = floor(proportion_vacuum_nodes * number_distance_nodes);
proportion_insulation_2_nodes = radius_insulation_2 / radius_total;
number_insulation_2_nodes = floor(proportion_insulation_2_nodes * number_distance_nodes);
proportion_metal_nodes = radius_metal / radius_total;
number_metal_nodes = floor(proportion_metal_nodes * number_distance_nodes);
node_sum = number_fluid_nodes + number_insulation_nodes + number_vacuum_nodes + number_insulation_2_nodes + number_metal_nodes;

% steady state matrix size
temperature_matrix = zeros(number_time_nodes, int16(node_sum));

% compute the radius & area between every node
current_radius = delta_x / 2;

% there are M-1 cross-sections
cross_areas = zeros(int16(node_sum)-1, 1);
for m = 1:number_distance_nodes-1
    
    width = 2 * sqrt(-1*current_radius^2 + radius_total^2);
    cross_areas(m) = width * height;
    %cross_areas(m) = 0.0128; % just for testing
    current_radius = current_radius + delta_x;

end

cross_areas(end) = cross_areas(end-1);

% fill in the steady state fluid nodes
%number_distance_nodes = double(number_distance_nodes);
%proportion_fluid_nodes = radius_fluid / radius_total;
%number_fluid_nodes = floor(proportion_fluid_nodes * number_distance_nodes);

% initial conditions
for i = 1:number_fluid_nodes
    temperature_matrix(1, i) = fluid_initial_temperature;
end

%proportion_insulation_nodes = radius_insulation / radius_total;
%number_insulation_nodes = floor(proportion_insulation_nodes * number_distance_nodes);

for i = 1:number_insulation_nodes
    temperature_matrix(1, number_fluid_nodes + i) = insulation_initial_temperature;
end

%proportion_vacuum_nodes = radius_vacuum / radius_total;
%number_vacuum_nodes = floor(proportion_vacuum_nodes * number_distance_nodes);

for i = 1:number_vacuum_nodes
    temperature_matrix(1, number_fluid_nodes + number_insulation_nodes + i) = vacuum_temperature;
end

%proportion_insulation_2_nodes = radius_insulation_2 / radius_total;
%number_insulation_2_nodes = floor(proportion_insulation_2_nodes * number_distance_nodes);

for i = 1:number_insulation_2_nodes
    temperature_matrix(1, number_fluid_nodes + number_insulation_nodes + number_vacuum_nodes + i) = insulation_2_initial_temperature;
end

%proportion_metal_nodes = radius_metal / radius_total;
%number_metal_nodes = floor(proportion_metal_nodes * number_distance_nodes);

for i = 1:number_metal_nodes
    temperature_matrix(1, number_fluid_nodes + number_metal_nodes + number_vacuum_nodes + number_insulation_2_nodes + i) = metal_initial_temperature;
end

%node_sum = number_fluid_nodes + number_insulation_nodes + number_vacuum_nodes + number_insulation_2_nodes + number_metal_nodes;

% step forward i-1 time steps
for i = 1:number_time_nodes-1
    current_distance_node = 1; % just for sanity
    % 0th node
    temperature_matrix(i+1, 1) = backward_conduction(temperature_matrix(i, 1), temperature_matrix(i, 2), k_fluid, cross_areas(1), ...
        rho_fluid, C_p_fluid, delta_x, delta_t);
    current_distance_node = 2;

    % fluid interior nodes
    for m = current_distance_node:number_fluid_nodes-1
        temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), ...
            temperature_matrix(i, m), temperature_matrix(i, m+1), k_fluid, ...
            cross_areas(m-1), k_fluid, cross_areas(m), rho_fluid, C_p_fluid, ...
            delta_x, delta_t);

        nodes_added = m-current_distance_node+1;
    end
    current_distance_node = current_distance_node + nodes_added;
    nodes_added = 0;

    % fluid-fluid-insulation boundary nodes (last fluid node)
    for m = current_distance_node:current_distance_node
        temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
            k_fluid, cross_areas(m-1), k_insulation, cross_areas(m), rho_fluid, C_p_fluid, ...
            delta_x, delta_t);

        nodes_added = m-current_distance_node+1;
    end
    current_distance_node = current_distance_node + nodes_added;
    nodes_added = 0;

    %% Inner insulation nodes section
    if number_insulation_nodes > 1 % changed 
        % fluid-insulation-insulation boundary nodes
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_fluid, cross_areas(m-1), k_insulation, cross_areas(m), rho_insulation, C_p_insulation, ...
                delta_x, delta_t);

            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % inner insulation nodes
        for m = current_distance_node:current_distance_node+number_insulation_nodes-3
            temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_insulation, cross_areas(m-1), k_insulation, cross_areas(m), rho_insulation, C_p_insulation, ...
                delta_x, delta_t);
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % metal-insulation-vacuum boundary node
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = conduction_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_metal, cross_areas(m-1), epsilon_insulation, cross_areas(m), rho_insulation, C_p_insulation, ...
                delta_x, delta_t);
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

    elseif number_insulation_nodes == 1 % changed
        % fluid-insulation-radiation boundary node
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = conduction_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_fluid, cross_areas(m-1), epsilon_insulation, cross_areas(m), rho_insulation, C_p_insulation, ...
                delta_x, delta_t);

            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;
    else
        %disp("The number of nodes you have is too few");
    end

    %% Vacuum nodes section
    if number_vacuum_nodes > 1
        % insulation-vacuum-vacuum boundary nodes
        for m = current_distance_node:current_distance_node
            % reza uses only emissivity of metal, so we're also doing that
            % here
            temperature_matrix(i+1, m) = dual_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), epsilon_metal, cross_areas(m), rho_air, C_p_air, ...
                delta_x, delta_t);

            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % inner vacuum nodes
        for m = current_distance_node:current_distance_node+number_vacuum_nodes-3
            temperature_matrix(i+1, m) = dual_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), epsilon_metal, cross_areas(m), rho_air, C_p_air, ...
                delta_x, delta_t);
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % vacuum-vacuum-metal boundary node
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = dual_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), epsilon_metal, cross_areas(m), rho_air, C_p_air, ...
                delta_x, delta_t);

            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

    elseif number_vacuum_nodes == 1
        % insulation1-vacuum-insulation2
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = dual_radiation(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), epsilon_metal, cross_areas(m-1), rho_air, C_p_air, ...
                delta_x, delta_t);
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;
    else
        % metal-metal-metal
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = 23;
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;
    end

    %% insulation layer 2 nodes section
    if number_insulation_2_nodes > 1
        % vacuum-insulation2-insulation2 boundary nodes
        for m = current_distance_node:current_distance_node
            % reza uses only emissivity of metal, so we're also doing that
            % here
            temperature_matrix(i+1, m) = radiation_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), k_insulation_2, cross_areas(m), rho_insulation_2, C_p_insulation_2, ...
                delta_x, delta_t);
            %temperature_matrix(i+1, m) = -5;
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % inner insulation2 nodes
        for m = current_distance_node:current_distance_node+number_insulation_2_nodes-3
            temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_insulation_2, cross_areas(m-1), k_insulation_2, cross_areas(m), rho_insulation_2, C_p_insulation_2, ...
                delta_x, delta_t);
            %temperature_matrix(i+1, m) = 10;
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % insulation2-insulation2-metal boundary node
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_insulation_2, cross_areas(m-1), k_insulation_2, cross_areas(m), rho_insulation_2, C_p_insulation_2, ...
                delta_x, delta_t);
            %temperature_matrix(i+1, m) = 5;
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

    elseif number_insulation_2_nodes == 1
        % vacuum-insulation2-metal
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = radiation_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), k_metal, cross_areas(m), rho_insulation_2, C_p_insulation_2, ...
                delta_x, delta_t);
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;
    else
        %disp("increase the proportion of insulation2 nodes, you don't have any");
    end

    %% Outer metal nodes section
    if number_metal_nodes > 1
        % vacuum-metal-metal boundary nodes
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = radiation_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                epsilon_metal, cross_areas(m-1), k_metal, cross_areas(m), rho_metal, C_p_metal, ...
                delta_x, delta_t);
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % inner metal nodes
        for m = current_distance_node:current_distance_node+number_metal_nodes-3
            temperature_matrix(i+1, m) = dual_conduction(temperature_matrix(i, m-1), temperature_matrix(i, m), temperature_matrix(i, m+1), ...
                k_metal, cross_areas(m-1), k_metal, cross_areas(m), rho_metal, C_p_metal, ...
                delta_x, delta_t);
            nodes_added = m-current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

        % metal-metal-convection boundary node
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = conduction_convection(temperature_matrix(i, m-1), temperature_matrix(i, m), ambient_temperature, ...
                k_metal, cross_areas(m-1), h_air, A_surface, rho_metal, C_p_metal, ...
                delta_x, delta_t); 
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;

    elseif number_metal_nodes == 1
        % vacuum_metal_convection
        for m = current_distance_node:current_distance_node
            temperature_matrix(i+1, m) = radiation_convection(temperature_matrix(i, m-1), temperature_matrix(i, m), ambient_temperature, ...
                epsilon_metal, cross_areas(m-1), h_air, A_surface, rho_metal, C_p_metal, ...
                delta_x, delta_t);
            nodes_added = m - current_distance_node+1;
        end
        current_distance_node = current_distance_node + nodes_added;
        nodes_added = 0;
    else
        disp("Too few metal nodes");
    end
end

mirror_matrix = flip(temperature_matrix, 2);
visual_matrix = [mirror_matrix, temperature_matrix];

%% animation
nodes = linspace(-double(node_sum), double(node_sum), double(node_sum*2));

% for i = 1:1000:size(visual_matrix, 1)
%     current_node = number_fluid_nodes;
%     xline(current_node, '-k');
%     hold on
%     xline(-current_node, '-k');
%     current_node = current_node + number_insulation_nodes;
%     xline(current_node, '-y');
%     xline(-current_node, '-y');
%     current_node = current_node + number_vacuum_nodes;
%     xline(current_node, '-c');
%     xline(-current_node, '-c');
%     current_node = current_node + number_insulation_2_nodes;
%     xline(current_node, '-k');
%     xline(-current_node, '-k');
%     current_node = current_node + number_metal_nodes;
%     xline(current_node, '-y');
%     xline(-current_node, 'y');
% 
%     plot(nodes, visual_matrix(i, :), '-b');
%     xlim([-(node_sum+5), node_sum+5]);
%     L = get(gca, 'XLim');
%     number_xticks = L(2)*2+1;
%     set(gca, 'XTick', linspace(L(1), L(2), number_xticks));
%     ylim([-20, 100])
%     xlabel('Node');
%     ylabel('Temperature (C)');
%     title_text = 'temperature distribution @ t = ' + string((i-1) * delta_t) + ' seconds';
%     title(title_text);
%     drawnow
% 
%     if i < size(visual_matrix, 1)
%         clf
%     end
% end

% final distribution
current_node = number_fluid_nodes;
xline(current_node, '-k');
hold on
xline(-current_node, '-k');
current_node = current_node + number_insulation_nodes;
xline(current_node, '-y');
xline(-current_node, '-y');
current_node = current_node + number_vacuum_nodes;
xline(current_node, '-c');
xline(-current_node, '-c');
current_node = current_node + number_insulation_2_nodes;
xline(current_node, '-k');
xline(-current_node, '-k');
current_node = current_node + number_metal_nodes;
xline(current_node, '-y');
xline(-current_node, 'y');


plot(nodes, tail(visual_matrix, 1), '-b');

xlim([-(node_sum+5), node_sum+5]);
L = get(gca, 'XLim');
number_xticks = L(2)*2+1;
set(gca, 'XTick', linspace(L(1), L(2), number_xticks));
ylim([-20, 100])
xlabel('Node');
ylabel('Temperature (C)');
title_text = 'temperature distribution @ t = ' + string((size(visual_matrix, 1)-1) * delta_t) + ' seconds';
title(title_text);