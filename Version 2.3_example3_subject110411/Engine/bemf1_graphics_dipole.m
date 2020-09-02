function [ ] = bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, type) 
%   Dipole(s) plot with several options
%   type  = 1 for a planar plot in the coronal plane (XZ)
%   type  = 2 for a planar plot in the sagittal plane(YZ)
%   type  = 3 for a planar plot in the transverse plane(XY)
%   type  = 4 for a 3D plot
%
%   Copyright SNM 2017-2018

% strdipolePplus      = strdipolePplus(1:neuron_ind:end);
% strdipolePminus     = strdipolePminus(1:neuron_ind:end);
% strdipoleCurrent    = strdipoleCurrent(1:neuron_ind:end);

    sz = 25; %   circle size in points squared
    
    if type == 0        % volume
        hold on;
        scatter3(strdipolePplus (:, 1), strdipolePplus (:, 2), strdipolePplus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[1 0 0],...
                                                                    'LineWidth',1.5)    
        scatter3(strdipolePminus (:, 1), strdipolePminus (:, 2), strdipolePminus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[0 0 1],...
                                                                    'LineWidth',1.5)
        for m = 1:size(strdipolePplus)      
            line([strdipolePplus(m, 1) strdipolePminus(m, 1)], [strdipolePplus(m, 2) strdipolePminus(m, 2)], [strdipolePplus(m, 3) strdipolePminus(m, 3)], 'LineWidth', 2, 'Color', 'c');           
        end
    end  
    
    if type == 1        % xy-plane
        scatter(strdipolePminus (:, 1), strdipolePminus (:, 2), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[0 0 1],...
                                                                    'LineWidth',1.5)
        scatter(strdipolePplus (:, 1), strdipolePplus (:, 2), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[1 0 0],...
                                                                    'LineWidth',1.5)       
        for m = 1:size(strdipolePplus, 1)
            line([strdipolePplus(m, 1) strdipolePminus(m, 1)], [strdipolePplus(m, 2) strdipolePminus(m, 2)], 'LineWidth', 2, 'Color', 'm') 
        end
    end  
    
    
    if type == 2        % xz-plane
        scatter(strdipolePplus (:, 1), strdipolePplus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[1 0 0],...
                                                                    'LineWidth',1.5)    
        scatter(strdipolePminus (:, 1), strdipolePminus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[0 0 1],...
                                                                    'LineWidth',1.5)
        for m = 1:size(strdipolePplus, 1)
            line([strdipolePplus(m, 1) strdipolePminus(m, 1)], [strdipolePplus(m, 3) strdipolePminus(m, 3)], 'LineWidth', 2, 'Color', 'm') 
        end
    end
    if type == 3        % yz-plane
        scatter(strdipolePplus (:, 2), strdipolePplus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[1 0 0],...
                                                                    'LineWidth',1.5)    
        scatter(strdipolePminus (:, 2), strdipolePminus (:, 3), sz, 'MarkerEdgeColor',[0 0 0],...
                                                                    'MarkerFaceColor',[0 0 1],...
                                                                    'LineWidth',1.5)
        for m = 1:size(strdipolePplus, 1)
            line([strdipolePplus(m, 2) strdipolePminus(m, 2)], [strdipolePplus(m, 3) strdipolePminus(m, 3)], 'LineWidth', 2, 'Color', 'm') 
        end
    end
end