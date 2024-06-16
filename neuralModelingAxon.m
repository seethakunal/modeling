%{
[t,Vm,AP]=MRGaxon('Nodes', Nodes,...
 'FiberDiameter', Diameter,...
 'PW', PulseWidth,...
 'TimeSpan', TimeSpan,...
 'Extracellular', Vext,...
 'PWonset', PWonset); 

 Usage:
   [t,VM]=MRGaxon(InputField,InputValue,...)
       t: vector of time points
       VM: matrix of membrane voltages at each time point (rows) and for
                        each node (columns)
   [t,VM,AP]=MRGaxon(InputField,InputValue,...)
       AP: indicates whether an action potential occured with a one or zero
                        if over-stimulation occured, AP will be -1
 
INPUT FIELDS:
---------------------------------------------------------------------------
 Required Inputs:
   'Nodes',[nodes]: Where nodes is a scalar (1x1) number of Nodes of
                        Ranvier in the axon
   'FiberDiameter',[diameter]: Where diameter is a scalar (1x1) fiber
                        diameter falling between 4 and 20 microns.
                        REQUIRED.
   'PW',[pw]: Where pw is a scalar (1x1) duriation of stimulation, in
                        microseconds. REQUIRED.
                        NOTE: if using multi electrode stimulation, this
                        must be a vector of the individual pulse widths
   'TimeSpan',[t_start t_stop]: Sets the time span, in microseconds.
                        Optionally, [t_start:dt:t_stop]. REQUIRED.
---------------------------------------------------------------------------
 Required Stimulus Inputs (choose 1):
   'Extracellular',[vext]: Where vext is a vector of extracellular
                        voltages at sequential Nodes of Ranvier, starting 
                        at Node 1 and ending at Node N, in millivolts.
                        For multi electrode stimulation, Vext must be a
                        cell array containing the individual Vext vectors.
                        NOTE: at least 1 stimulus type is REQUIRED.
   'Intracellular',[Iinj]: Where Iinj is a vector of the injected currents
                        at sequential Nodes of Ranvier, starting at
                        Node 1 and ending at Node N, in microamps.
                        NOTE: at least 1 stimulus type is REQUIRED.
---------------------------------------------------------------------------
 Optional Stimulus Inputs:
   'PWonset',[delay]: Where delay is a scalar (1x1) delay before stimulus
                        onset, in microseconds.
                        NOTE: if using multi electrode stimulation, this
                        must be a vector of the individual pulse width onsets
   'Frequency',[Frequency]: Where the input is a scalar (1x1) in Hertz
                        to set the frequency of stimulus
                        NOTE: can only be set with biphasic stimulus
                        NOTE: if using multi electrode stimulation, this
                        must be a vector of the individual frequencies
   'Biphasic': Can be set to 1 to indicate stimulus is biphasic
                        NOTE: If using multi electrode stimulation, this must
                        be a vector containing the values for each electrode
   'PWRatio',[Ratio]: Where ratio is a scalar (1x1) indicating the ratio
                        AnodePW/CathodePW. This will scale the anode PW
                        with respect to the Cathode PW and adjust the
                        anode voltage apropriately
                        NOTE: can only be set with biphasic stimulus
                        NOTE: if using multi electrode stimulation, this must
                        be a vector of the individual pulse width ratios
---------------------------------------------------------------------------
 Optional Stimulus Modulation Inputs:
   'PWM',@yourPWfunction: Where the input is an @ symbol followed by the name
                        of your function, defined in its own .m file. Your function
                        must take a time as input. Adjusted pulse width will be
                        multiplied by the output of your function.
                        NOTE: can not be used with PWRatio or multiple electrode
                        stimulation.
   'PAM',@yourPAfunction: Where the input is an @ symbol followed by the name
                        of your function, defined in its own .m file. Your function
                        must take a time as input. Stimulation voltages will be
                        multiplied by the output of your function.
                        NOTE: can not be used with PWRatio or multiple electrode
                        stimulation.
   'FM',@yourFMfunction: Where the input is an @ symbol followed by the name
                        of your function, defined in its own .m file. Your function
                        must take a time as input. The stimulation frequency will be
                        multiplied by the output of your function.
                        NOTE: user must set a base frequency.
                        NOTES: can not be used with multiple electrode
                        stimulation.
---------------------------------------------------------------------------
 Additional Optional Inputs:
   'Axons',[axons]: Where axons is a scalar (1x1) number of axons
                        to be simulated
   'Sensory': Can be set to 1 to indicate the axons simulated are sensory.
                If not set, the axon defaults to motor axon dynamics.
   'Plot': can be set to 'On'. Otherwise plot does not show.

%}

%{
Node Length = 1 µm
• Myelin (or “Internode”) Length = 2.501
0.001593549+0.032342282𝑒𝑒−0..42(𝐹𝐹𝐹𝐹𝐹𝐹𝐹𝐹𝐹𝐹 𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷𝐷) + 0.8
• The conductivity of the extracellular media, σ = 2E-6 S/µm
• Total number of nodes = 41
• The electrode is centered over the center node (node 21)
%}

nodelength = 1; % um
fiberdiameter = 10; % um
internodalspacing = 2.501/(0.001593549 + 0.032342282*exp(-.42*fiberdiameter)) +.8; % um
conductivity = 2*10^-6; % S/um
numNodes = 41;
centerNode = 21; % node where stimulation occurs

%% q1
heightelectrode = 100:50:1000; % um
PulseWidth = 40;
TimeSpan = [0 10000];

Ithresh = [];
its = [];

for i = 1:length(heightelectrode) % um
    AP = 0;
    % cathodic, Istim varies between 0 and -5.0 mA
    Istim = .1;
    change = -.1;
    iter = 0;
    while (AP ~= 1)
        
        iter = iter + 1;

        if (AP ~= 0)
            Istim = Istim - change;
        end

        if (AP == 0)
            Istim = Istim + change;
        end

        [t,Vm,AP]=MRGaxon('Nodes', numNodes,...
         'FiberDiameter', fiberdiameter,...
         'PW', PulseWidth,...
         'TimeSpan', TimeSpan,...
         'Extracellular', Vextracellular(heightelectrode(i), centerNode, numNodes, Istim, conductivity, internodalspacing));
    end

   % after loop our Istim is our Ithreshold
   Ithresh(i) = Istim;
   
end
figure(1)
subplot(1,2,1)
plot(heightelectrode, abs(Ithresh))
title('Cathodic stimulation, over central node of ranvier')
xlabel('height of electrode (um)')
ylabel('Absolute value of Current threshold stimulation (mA)')
%{
 [t,VM,AP]=MRGaxon(InputField,InputValue,...)
       AP: indicates whether an action potential occured with a one or zero
                        if over-stimulation occured, AP will be -1
%}

%% q2
% now the stimulation is centered over the myelin sheath of internode 21,
% this change is accounted for in the function Vextracellular2

heightelectrode = 100:50:1000; % um
PulseWidth = 40;
TimeSpan = [0 10000];

Ithresh = [];
its = [];

for i = 1:length(heightelectrode) % um
    AP = 0;
    % cathodic, Istim varies between 0 and -5.0 mA
    Istim = .1;
    change = -.1;
    iter = 0;
    while (AP ~= 1)
        
        iter = iter + 1;

        if (AP ~= 0)
            Istim = Istim - change;
        end

        if (AP == 0)
            Istim = Istim + change;
        end

        [t,Vm,AP]=MRGaxon('Nodes', numNodes,...
         'FiberDiameter', fiberdiameter,...
         'PW', PulseWidth,...
         'TimeSpan', TimeSpan,...
         'Extracellular', Vextracellular2(heightelectrode(i), centerNode, numNodes, Istim, conductivity, internodalspacing));
    end

   Ithresh(i) = Istim;
   
end
figure(1)
subplot(1,2,2)
plot(heightelectrode, abs(Ithresh))
title('Cathodic stimulation, over central internode (myelated portion)')
xlabel('height of electrode (um)')
ylabel('Absolute value of Current threshold stimulation (mA)')

%% q3, repeat procedures in q1 and q2, except anodic stimulation
% anodic stimulation varies between 0 and +20 mA
% copy and paste the procedures

heightelectrode = 100:50:1000; % um
PulseWidth = 40;
TimeSpan = [0 10000];

Ithresh = [];
its = [];

for i = 1:length(heightelectrode) % um
    AP = 0;
    % cathodic, Istim varies between 0 and -5.0 mA
    Istim = -.1;
    change = .1;
    iter = 0;
    while (AP ~= 1)
        
        iter = iter + 1;

        if (AP ~= 0)
            Istim = Istim - change;
        end

        if (AP == 0)
            Istim = Istim + change;
        end

        [t,Vm,AP]=MRGaxon('Nodes', numNodes,...
         'FiberDiameter', fiberdiameter,...
         'PW', PulseWidth,...
         'TimeSpan', TimeSpan,...
         'Extracellular', Vextracellular(heightelectrode(i), centerNode, numNodes, Istim, conductivity, internodalspacing));
    end

   % after loop our Istim is our Ithreshold
   Ithresh(i) = Istim;
   
end
figure(2)
subplot(1,2,1)
plot(heightelectrode, abs(Ithresh))
title('Anodic stimulation, over central node of ranvier')
xlabel('height of electrode (um)')
ylabel('Absolute value of Current threshold stimulation (mA)')

heightelectrode = 100:50:1000; % um
PulseWidth = 40;
TimeSpan = [0 10000];

Ithresh = [];
its = [];

for i = 1:length(heightelectrode) % um
    AP = 0;
    % cathodic, Istim varies between 0 and -5.0 mA
    Istim = -.1;
    change = .1;
    iter = 0;
    while (AP ~= 1)
        
        iter = iter + 1;

        if (AP ~= 0)
            Istim = Istim - change;
        end

        if (AP == 0)
            Istim = Istim + change;
        end

        [t,Vm,AP]=MRGaxon('Nodes', numNodes,...
         'FiberDiameter', fiberdiameter,...
         'PW', PulseWidth,...
         'TimeSpan', TimeSpan,...
         'Extracellular', Vextracellular2(heightelectrode(i), centerNode, numNodes, Istim, conductivity, internodalspacing));
    end

   Ithresh(i) = Istim;
   
end
figure(2)
subplot(1,2,2)
plot(heightelectrode, abs(Ithresh))
title('Cathodic stimulation')
xlabel('height of electrode (um)')
ylabel('Absolute value of Current threshold stimulation (mA)')

%%
function [Vext] = Vextracellular(Dcenter, centernode, node, Istim, conductivity, internodalspacing)
% returns the Vextracellular at the given node 'node'
for n=1:node
numNodesBetween = abs(centernode - n);
Dnode = sqrt((numNodesBetween*internodalspacing)^2 + Dcenter^2);
Vext(n) = Istim/(4*pi*conductivity*Dnode);
end
end

%%
function [Vext] = Vextracellular2(Dcenter, centernode, node, Istim, conductivity, internodalspacing)
% returns the Vextracellular at the given node 'node'
% in this function, the stimulation is centered above the myelin sheath of
% the central node
for n=1:node
numNodesBetween = abs(centernode - n);
Dnode = sqrt((numNodesBetween*internodalspacing + internodalspacing/2)^2 + Dcenter^2);
Vext(n) = Istim/(4*pi*conductivity*Dnode);
end
end