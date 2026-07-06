% The function for transmission energy
function Etx= txEnergy(packet_length,d) %% Transmission  model(mp/fs)    
    E_elec=50e-9;  %% Energy for tx/rx one bit
    Efs=10e-12; %% amplification energy for free space transmisshion
    Emp=0.0013e-12; %% amplification energy for multipath jule/b/m4
    do=sqrt(Efs/Emp); %% Threshold distance
    if d < do
        Etx= packet_length*(E_elec + Efs*d^2);
    else
        Etx= packet_length*(E_elec + Emp*d^4);
    end
end