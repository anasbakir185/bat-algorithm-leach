% BA-LEACH single-hop protocol code Starts from here
%% part 1
clear
alive_round = zeros(1,1000);
stream1 = RandStream('mt19937ar', 'Seed', 5);
stream2 = RandStream('mt19937ar', 'Seed', 6);
Eag = 5e-9;
numbats=2; % number of bats
group = []; % the cluster members without the CH
fmin = 1; % for bat algorithm
fmax = 20;
rounds = 10;
num_CHs = zeros(1,rounds);
died_round = [];
dead = 0;
p=0.1; 
Eo = 0.5;
E_elec = 50e-9;
packet_size = 4000;
Erx = E_elec*packet_size;
n = 250;
s = struct('xd',zeros(1,n),'yd',zeros(1,n),'E',zeros(1,n),'type',repmat('N',n),'alive',zeros(1,n),'nearestCH',zeros(1,n),'dist_to_CH',zeros(1,n)); %preallocate a Struct for the nodes and its properties
for i = 1:n
    s(i).xd = 250*rand(stream1); 
    s(i).yd = 250*rand(stream1);
    s(i).E = Eo;
    s(i).alive = true;
    s(i).type = 'N';
    s(i).dist_to_CH = 0; % the distance between the node and the CH
    
end
s(n+1).xd = 250*0.5;
s(n+1).yd = 250*0.5;
s(n+1).type = 'S';
CH_list = [];
CHs = [];
CHs2 = []; % CHs reselcted by BAT algorithm
avgE = zeros(1,rounds); %save the avg energy of each round
%% part 2
for r = 1:rounds %Total time of the network 
    for i = 1:n
        if s(i).E <= 0 && s(i).alive 
            dead = dead+1;
            s(i).alive = false;
            s(i).E = 0;
            died_round = [died_round; r dead];% at what round each sensor died 
        end
        alive_round(r) = sum([s.alive]);
        s(i).type = 'N';
    end
    Tn = p/(1-p*(mod(r,1/p))); % Threshold
    if (rem(r,10) ~= 0) 
        CH_list = [CH_list,CHs2];
    else
        CH_list = [];
    end
    CHs = []; 
    k = 0;
    again = 0; 
    while isempty(CHs)   
        again=again+1;
        if again >= 5 
            CH_list = [];
        end
        for j = 1:n %setup phase
            if (s(j).alive && (rand(stream1) <= Tn)&& not(ismember(j,CH_list)))%% check if it was CHs in the last 10 rounds/ make elections when CHs<10
                s(j).type = 'CH';
                k = k+1;
                CHs(k) = j;
            end
        end
    end
    num_CHs(r) = numel(CHs);
    CHs2 = CHs;
    members = zeros(1,k);
    totalE = 0; % Total Energy
    cluster = cell(1,numel(CHs)) ; % the node indicies(index) of each cluster
    for i = 1:n % nodes
        s(i).nearestCH = 0; 
        dist_to_CH = [];
        if strcmp(s(i).type,'N') && s(i).alive 
            for ii = 1:numel(CHs), j = CHs(ii);
                dist_to_CH(ii) = dist([s(i).xd s(i).yd],[s(j).xd s(j).yd]); % save the distance to all CHs to cose the min
            end
            [d,argmin] = min(dist_to_CH);
            CH = CHs(argmin);
            members(argmin) = members(argmin)+1;
            s(i).nearestCH = CH;
            s(i).dist_to_CH = d;
            cluster{argmin} = [cluster{argmin} i]; % the node indices of each cluster (add members to each cluster)
            s(i).E = s(i).E - txEnergy(packet_size,d);
            s(CH).E = s(CH).E - Erx;
        end
    end
    bestCH = CHs; %the best CHs in this round/ to be used at the end of the iterations
    bestFitF = repmat(inf,[1 numel(CHs)]);%variable. best fittness function in this round
    for b = 1:numbats % generate bats
        group = cluster;
        v = 10*[rand(stream2) rand(stream2)]; %velocity (BA)
        for j = 1:numel(cluster) % nodes in the cluster
            group2{j} = [cluster{j},CHs(j)]; %for plotting the the old cluster head 
            if numel(group{j}) <= 2 % if the cluster members is less than 3 not bat algorithm is needed
                continue
            end
            newCH = group{j}(randi(stream2,numel(group{j})));% choose a node from the cluster that hasnt been used by the bat
            newPos=[s(newCH).xd s(newCH).yd];
            bestPos = [s(bestCH(j)).xd s(bestCH(j)).yd]; % the postion of the node with the best fitness function
            ra = 0.1; %pulse rate
            A = 1;%loudness
            for search = 1:10 %(iterations) for each bat
                ra = 1.26*ra;
                A = 0.8*A;
                for i = cluster{j} %set the new info for each cluster
                    s(i).nearestCH = newCH;
                    s(i).type = 'N';
                    s(i).dist_to_CH = dist([s(i).xd s(i).yd],[s(newCH).xd s(newCH).yd]);
                end
                s(newCH).type = 'CH';
                s(newCH).nearestCH = 0;
                idx=[s.nearestCH] == newCH;   % indexes for distances to CH j for each cluster
                total_dis_to_CH = [s.dist_to_CH]; %distances to CHs from all nodes 
                total_dis_to_CH = total_dis_to_CH(idx);
                avg_dis(j) = sum(total_dis_to_CH)/members(j); % sum of distances to CH(i) divided by number of members=avg distance CH for each cluster
                d2(j) = dist([s(newCH).xd s(newCH).yd],[s(n+1).xd s(n+1).yd]);       % Distance from CH to BS
                CHe(j) = s(newCH).E ;%CHs energy
                fitF(j) = 0.1*(avg_dis(j)/70)+0.8*(d2(j)/353.5)+0.1*(1-(CHe(j)/0.5));  %fitness fuction of the node/ for the normalized values ()
                if fitF(j) < bestFitF(j)  % the best fit function over the iterations (one round)
                    bestFitF(j) = fitF(j);
                    bestCH(j) = newCH;
                    bestPos = [s(bestCH(j)).xd s(bestCH(j)).yd];
                end
                group{j}(group{j} == newCH) = []; % prevent the the chosen CH from being a future option by the bat
                if isempty(group{j})
                    continue
                end
                f = fmin+(fmax-fmin)*rand(stream2);
                v = v+(newPos-bestPos)*f; %global search
                newPos = newPos-v; 
                if rand(stream2) > ra && ~isequal(newPos,bestPos)
                    newPos = newPos+(2*rand(stream2)-1)*A*40; % local search
                end
                [~,argmin] = min(sqrt(([s(group{j}).xd] - (newPos(1))).^2+ ([s(group{j}).yd] -(newPos(2))).^2)); % the node that is closest to the next movement
                newCH = group{j}(argmin); % the new node that we moved to 
                newPos = [s(newCH).xd s(newCH).yd];
            end
        end 
    end
    CHs = bestCH; % the chs cosen by the bat are the new chs
    for j = 1:numel(CHs)
        s(CHs(j)).E = s(CHs(j)).E -( Eag * packet_size * (members(j)));
        d2(j) = dist([s(CHs(j)).xd s(CHs(j)).yd],[s(n+1).xd s(n+1).yd]);       % Distance from CH to BS
        s(CHs(j)).E = s(CHs(j)).E-txEnergy((members(j)+1)*packet_size,d2(j));    % residual energy after aggregating to the BS \the packet length depends on the number of cluster members
    end
    avgE(r) = sum([s.E])/n;% avarge residual E in each round (for plotting)
end
%% part 3 ploting
figure(4)
grid on
title('Network Lifetime for BA-LEACH Sinle-hop')
data_bar = [died_round(1) died_round(30) died_round(50) died_round(70)];       % to plot the bar for dead nodes
bar(data_bar)
xlabel('Network Life Time Events');
ylabel('Round Number');
set(gca,'Xticklabels',{'First node die', '30% nodes die','50% nodes die','70% nodes die'})
figure(2)
grid on
title('Average Residual Energy of Nodes in LEACH Protocol')
plot(avgE)
axis([0 rounds 0 Eo])
xlabel('Round Number');
ylabel('Average Residual Energy');
figure(3)
grid on
plot(alive_round)
title('Number of Alive nodes over Simulation')
axis([0 rounds 0 n])
xlabel('Round Number');
ylabel('Number of Alive Nodes');
%BA-LEACH single-hop protocol code ends here   