function [ elementList,NodeList ] = ElementsUnderInfluenceRadius( SensorLocation,influence_Radius )
load data.mat data;
load elementNodes.mat elementNodes;
numberNodesX   =data(9)+1; %along X-direction(width)
% Finding all the nodes in the influence radius of the sensor
NodeList = zeros((2*influence_Radius+1)*(2*influence_Radius+1),1);
ListCounter =1;
for n=-2:2
    for radius_counter = 0:4
        NodeList(ListCounter,:) = (SensorLocation+(n*numberNodesX)...
                                    -(influence_Radius)+radius_counter);
        ListCounter = ListCounter+1;
    end
end
% with the nodes in the influence radius finding the elements in the range
elementList = zeros((2*influence_Radius)*(2*influence_Radius),1);
elementList_counter =1;
for node1_counter = [NodeList(1:4,1)'...
                        NodeList(6:9,1)'...
                        NodeList(11:14,1)'...
                        NodeList(16:19,1)']
    node = zeros(1,4); %4 nodes for one element
    node(1,1) = node1_counter;
    node(1,2) = node1_counter+1;
    node(1,3) = node1_counter+ numberNodesX+1;
    node(1,4) = node1_counter+ numberNodesX;
    for check_counter =1:length(elementNodes)
    if isequal(node,elementNodes(check_counter,:))
        elementList(elementList_counter,1) = check_counter;
        elementList_counter = elementList_counter+1;
    end
    end
end
end

