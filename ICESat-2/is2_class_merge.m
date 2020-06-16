function [class_consol] = is2_class_merge(class)

% A script designed to simplify the ICESat-2 photon classifications
% to be independent of surface type.
% class is a 5 x N array directly sourced from the IS2 variable "signal_conf_ph"
% class_consol is a 1 x N array with the surface type dependence removed

varsize = size(class);
if varsize(1) ~= 5 | any(class>4) | any(class<-1)
    error('Invalid variable input. Array needs to be of size 5 x N with identifiers [-1, 4].')
end


for i = 1:length(class)
    class_consol(i) = max([class(1,i) class(2,i) class(3, i) class(4,i) class(5,i)]);
end
