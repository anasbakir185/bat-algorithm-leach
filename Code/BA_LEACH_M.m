% BA-LEACH multi-hop protocol code Starts from here
%% part 1
clear
alive_round = zeros(1,1000);
stream1 = RandStream('mt19937ar', 'Seed', 5);
stream2 = RandStream('mt19937ar', 'Seed', 6);
num_chs = zeros(1000);
numbats = 2;
group = []; 
fmin = 1;
fmax = 20;
rounds = 1000;
died_round = [];
dead = 0;
p = 0.1; 
Eo = 0.5;
E_elec = 50e-9;  
packet_size = 4000; 
Erx = E_elec*packet_size;
Eag = 5e-9;
n = 250;
s = struct('xd',zeros(1,n),'yd',zeros(1,n),'E',zeros(1,n),'type',repmat('N',n),'alive',zeros(1,n),'nearestCH',zeros(1,n),'dist_to_CH',zeros(1,n)); %preallocate a Struct for the nodes and its properties
s(n+1).xd = 250*0.5;
s(n+1).yd = 250*0.5;
s(n+1).type = 'S';
for i = 1:n
    s(i).xd = 250*rand(stream1); %% nodes locations
    s(i).yd = 250*rand(stream1);
    s(i).E = Eo; %% energy for each node
    s(i).alive = true; %% assign true using for loop
    s(i).type = 'N';
    s(i).dist_to_CH = 0;
    s(i).dist_to_sink = dist([s(i).xd s(i).yd],[125 125]); % distance to the sink
end 
CH_list = [];
CHs = [];
CHs2 = [];
avgE = zeros(1,rounds); %save the avg energy of each round
%% part 2
for r = 1:rounds %Total time of the network    for i=1:n %% return all sensors to be nodes
    for i = 1:n
        if s(i).E <= 0. && s(i).alive
            dead = dead+1;
            s(i).alive = false;
            s(i).E = 0;
            died_round = [died_round; r dead];
        end
        alive_round(r) = sum([s.alive]);
        s(i).type = 'N';
    end
    Tn = p/(1-p*(mod(r,1/p)));
    if (rem(r,10) ~= 0) 
        CH_list = [CH_list,CHs2];
    else
        CH_list = [];
    end
    CHs = []; 
    k = 0;
    again = 0; 
    while isempty(CHs)
        again = again+1;
        if again >= 5 
            CH_list = [];
        end
        for j = 1:n %setup phase
            if (s(j).alive && ((rand(stream1)) <= Tn)&& not(ismember(j,CH_list)))%% check if it was CHs in the last 10 rounds/ make elections when CHs<10
                s(j).type = 'CH';
                k = k+1;
                CHs(k) = j;
            end
        end
    end
    colors = lines(numel(CHs));
    num_chs(r) = numel(CHs);
    CHs2 = CHs;
    members = zeros(1,k);
    totalE = 0; % Total Energy
    cluster = cell(1,numel(CHs)) ; 
    for i = 1:n
        s(i).nearestCH = 0;
        s(i).agg = 0;
        dist_to_CH = []; 
        if strcmp(s(i).type,'N') && s(i).alive 
            for ii = 1:numel(CHs), j = CHs(ii);
                dist_to_CH(ii) = dist([s(i).xd s(i).yd],[s(j).xd s(j).yd]);
            end
            [d,argmin] = min(dist_to_CH);
            CH_id = CHs(argmin);
            members(argmin) = members(argmin)+1; %how many members in each cluster%
            s(i).nearestCH = CH_id;
            s(i).dist_to_CH = d;
            cluster{argmin} = [cluster{argmin} i];
            s(i).E = s(i).E - txEnergy(packet_size,d);
            s(CH_id).E = s(CH_id).E - Erx;
        end
    end
    bestCH = CHs;
    bestFitF = repmat(inf,[1 numel(CHs)]);
    for b = 1:numbats
        group = cluster;
        v = 50*[rand(stream2) rand(stream2)];
        for j = 1:numel(cluster)
            group2{j} = [cluster{j},CHs(j)];
            if numel(group{j}) <= 2
                continue
            end
            newCH = group{j}(randi(stream2,numel(group{j})));
            newPos = [s(newCH).xd s(newCH).yd];
            bestPos = [s(bestCH(j)).xd s(bestCH(j)).yd];
            ra = 0.1; %pulse rate
            A = 1;%loudness
            for search = 1:10
                ra = 1.26*ra;
                A = 0.8*A;
                for i = cluster{j} 
                    s(i).nearestCH=newCH;
                    s(i).type='N';
                    s(i).dist_to_CH=dist([s(i).xd s(i).yd],[s(newCH).xd s(newCH).yd]);
                end
                s(newCH).type = 'CH';
                s(newCH).nearestCH = 0;
                idx=[s.nearestCH] == newCH;
                total_dis_to_CH = [s.dist_to_CH];
                total_dis_to_CH = total_dis_to_CH(idx);
                avg_dis(j) = sum(total_dis_to_CH)/members(j); 
                d2(j) = dist([s(newCH).xd s(newCH).yd],[s(n+1).xd s(n+1).yd]);
                CHe(j) = s(newCH).E ;%CHs energy
                fitF(j) = 0.1*(avg_dis(j)/150)+0.8*(d2(j)/353.5)+0.1*(1-(CHe(j)/0.5));  %fitness fuction/ for the normalized values ()
                if fitF(j) < bestFitF(j)
                    bestFitF(j) = fitF(j);
                    bestCH(j) = newCH;
                    bestPos = [s(bestCH(j)).xd s(bestCH(j)).yd];
                end
                group{j}(group{j} == newCH) = [];
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
                newCH = group{j}(argmin);
                newPos = [s(newCH).xd s(newCH).yd]; 
            end
        end
    end
    CHs = bestCH;
    avgE(r) = 0.5;% avarge residual E in each round
    dist_to_N_CH = cell(1,numel(CHs));
    dist_to_S = [];
    for j = 1:numel(CHs)
        s(CHs(j)).packets = members(j)+1; %pakets to be transmitted by one CHs (includes its cluster and other clusters) 
        for i = 1:numel(CHs)
            dist_to_N_CH{j}(i) = dist([s(CHs(j)).xd s(CHs(j)).yd],[s(CHs(i)).xd s(CHs(i)).yd]);
        end
        dist_to_S(j) = dist([s(CHs(j)).xd s(CHs(j)).yd],[s(n+1).xd s(n+1).yd]); %Distance to sink;
        dist_to_N_CH{j}(j) = dist_to_S(j);
    end
    for j = 1:numel(CHs) %fitness function for Multihop transmission
        cost = (dist_to_N_CH{j}/dist_to_S(j)<=0.7) .* (dist_to_S/dist_to_S(j)<=0.7).*([s(CHs).E]/s(CHs(j)).E>=0.8).*(dist_to_S(j)./(dist_to_N_CH{j}));
        [~, argmin] = max(cost);
        if argmin ~= 1  && dist_to_S(j)>60 % forward to relay nodes for nodes thats are bit far from the sink
            d2(j) = dist_to_N_CH{j}(argmin); % The distance between the ch and the intermediate cluster head
            s(CHs(argmin)).agg = s(CHs(argmin)).agg+members(j);
            s(CHs(argmin)).packets = s(CHs(argmin)).packets+s(CHs(j)).packets;% adding the number of packets to the intermediate ch
            d_id(j) = CHs(argmin); % the id of distanation CH to forward to
        else % else, forward to the sink 
            d2(j) = dist_to_S(j);
            d_id(j) = n+1; 
        end
    end
    for j = 1:numel(CHs) 
        s(CHs(j)).E = s(CHs(j)).E - Eag * packet_size * (s(CHs(j)).agg);
        s(CHs(j)).E = s(CHs(j)).E-txEnergy(s(CHs(j)).packets*packet_size,d2(j)); 
    end
    avgE(r) = sum([s.E])/n;
end
%% part 3 ploting
figure(5)
grid on
title('Network Lifetime for BA-LEACH Multi-hop')
data_bar = [died_round(1) died_round(30) died_round(50) died_round(70)];       % to plot the bar for dead nodes
bar(data_bar)
xlabel('Network Life Time Events');
ylabel('Round Number');
set(gca,'Xticklabels',{'First node die', '30% nodes die','50% nodes die','70% nodes die'})
figure(2)
grid on
title('Average Residual Energy of Nodes in LEACH Protocol')
plot(avgE)
legend('LEACH','BA-LEACH single-hop','BA-LEACH Multi-hop','Location','northeast')
%%%% use the legend after running the three simulations
axis([0 rounds 0 Eo])
xlabel('Round Number');
ylabel('Average Residual Energy');
figure(3)
grid on
legend('LEACH','BA-LEACH single-hop','BA-LEACH
Multi-hop','Location','northeast')
%% use the legend after running the three simulations
plot(alive_round)
title('Number of Alive nodes over Simulation')
axis([0 rounds 0 n])
xlabel('Round Number');
ylabel('Number of Alive Nodes');
%BA-LEACH multihop protocol code ends here

