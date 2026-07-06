% LEACH protocol code Starts from here
%% part 1 
clear
stream1 = RandStream('mt19937ar', 'Seed', 5);% random seed for keeping the same results
stream2 = RandStream('mt19937ar', 'Seed', 6);
num_CHs = zeros(1,1200); % number of CHs selected in each round 
alive_round = zeros(1,1200); % variable for keeping the number of alive node for each round
rounds = 800; % number of rounds
died_round =[]; %variable to know at which round each sensor died
dead = 0;%the number of currently dead sensors
p = 0.1; % percentage of CH
Eo = 0.5; % energy for each node
E_elec = 50e-9;  % Energy for tx/rx one bit
packet_size = 4000; %bits
Erx = E_elec*packet_size; %% Eneregy for recieving one packet
Eag = 5e-9; % Energy for aggregation
n = 250; % number of nodes
s = struct('xd',zeros(1,n),'yd',zeros(1,n),'E',zeros(1,n),'type',repmat('N',n),'alive',zeros(1,n),'nearestCH',zeros(1,n)); %preallocate a Struct for the nodes and its properties
for i = 1:n % 
    s(i).xd = 250*rand(stream1); %% nodes locations
    s(i).yd = 250*rand(stream1);
    s(i).E = Eo; %% energy for each node
    s(i).alive = true; %% the condition of the node(true = alive, false =dead)
    s(i).type = 'N'; % the type of the node (S = sink, N= Non-CH sensor, CH = cluster head)
    
end
s(n+1).xd = 250*0.5;%SINK node postion
s(n+1).yd = 250*0.5;
s(n+1).type = 'S';
CH_list = []; % variable for saving the CHs used in the last 1/p rounds
CHs = [];% variable for saving the CHs used in teh current round
avgE = zeros(1,rounds); %save the avg energy of each round
%% part 2 Protocol Operation
for r = 1:rounds %Total time of the network 
    r
    for i = 1:n %% return all sensors to be nodes
        if s(i).E <= 0 && s(i).alive % check if the node died
            dead = dead+1; % number of currently dead nodes
            s(i).alive = false; % change the condition for the node to dead
            s(i).E = 0; % avoid negative values 
            died_round = [died_round; r dead];% at what round each sensor died (imprtant for ploting)
        end
        alive_round(r) = sum([s.alive]); %number of alive nodes in each round (important for ploting)
        s(i).type = 'N'; % return all nodes to be non-CH nods
    end
    Tn = p/(1-p*(mod(r,1/p))); % Threshold for selecting CHs
    % keep the CH from the last 1/p rounds in a list
    if (rem(r,10) ~= 0) 
        CH_list = [CH_list,CHs];
    else 
        CH_list = [];
    end
    CHs = []; 
    k = 0;
    again = 0; %number of times to retry the CH selection process (in case no CHs were selected)
    while isempty(CHs)  % retry if no CHs in this round
        again=again+1;
        if again >= 5 
            CH_list =[]; % if after five selection rounds no nodes has selected it self as CH,forget about CH_list 
        end
        for j = 1:n %setup phase
            if (s(j).alive && (rand(stream1) <= Tn) && not(ismember(j,CH_list)) )%% check if it was CHs in the last 10 rounds
                s(j).type='CH';
                k = k+1;
                CHs(k) = j; % adds the node to the array of CHs
            end
        end
    end
    colors = lines(numel(CHs));   % K clusters  K different colors
    num_CHs(r) = numel(CHs);
    members = ones(1,k);%the members of each cluster
    totalE = 0; % total residual energy of sensors
    for i = 1:n % nodes
        s(i).nearestCH = 0; %variable for saving the nearest CH for the sensor
        dist_to_CH = []; % variable for saving the distance to the nearest CH  
        if strcmp(s(i).type,'N') && s(i).alive % alive non-CH nodes
            for ii = 1:numel(CHs), j = CHs(ii);% find the closest CHs to the node/ iterate through CHs and their indicies
                if j ~= 0 
                    dist_to_CH(ii) = dist([s(i).xd s(i).yd],[s(j).xd s(j).yd]);
                end
            end
            [d,argmin] = min(dist_to_CH);% the closest CH and the distance to it
            CH_id = CHs(argmin);% the id of the closest CH
            members(argmin) = members(argmin)+1; %members in each cluster%
            s(i).nearestCH = CH_id;
            s(i).E=s(i).E - txEnergy(packet_size,d);% Residual Energy of each node
            s(CH_id).E = s(CH_id).E - Erx;%for recieving the packet from the non-ch node
        end
        totalE = totalE+s(i).E; % Total residual E in this round
    end
    avgE(r) = totalE/n;% avarge residual E in all rounds  
    for j = 1:numel(CHs) % the energy dissipated by CHs for transmitting and aggregating packets
        s(CHs(j)).E = s(CHs(j)).E - (Eag * packet_size * members(j)); %residual energy after transmitting
        d2 = dist([s(CHs(j)).xd s(CHs(j)).yd],[s(n+1).xd s(n+1).yd]);       % Distance from CH to BS
        s(CHs(j)).E = s(CHs(j)).E-txEnergy(members(j)*packet_size,d2);    % residual energy after transmitting to the BS \the packet length depends on the number of cluster members   
    end
end
%% part 3 ploting
figure(1)
grid on 
title('Network Lifetime for LEACH Protocol')
data_bar = [died_round(1) died_round(30) died_round(50) died_round(70)];       % to plot the bar for dead nodes
bar(data_bar)
xlabel('Network Life Time Events');
ylabel('Round Number');
set(gca,'Xticklabels',{'First node die', '30% nodes die','50% nodes die','70% nodes die'})
figure(2)
grid on
hold on
title('Average Residual Energy of Nodes in LEACH Protocol')
plot(avgE)
axis([0 rounds 0 Eo])
xlabel('Round Number');
ylabel('Average Residual Energy');
figure(3)
grid on
hold on
plot(alive_round)
title('Number of Alive nodes over Simulation')
axis([0 rounds 0 n])
xlabel('Round Number');
ylabel('Number of Alive Nodes');
%%%LEACH protocol code ends here