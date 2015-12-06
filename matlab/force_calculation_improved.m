%Improved Force Calculation, with Verlet Cell Algorithm to Maintain
%Neighbor Lists.
function [f,u]=force_calculation_improved(N,r,Ls,rc2)

            
        
            L_cell = 4.e0;
            %Number of cells per direction
            n_cells = round(Ls/L_cell);
            L_cell = Ls/n_cells;
            %Number of cells per unit length.
           n_per_length = n_cells/Ls;
           %Total number of cells 
           %extended by 2 in every direction.
           n_cells = n_cells + 2;
           n_cells_z_level = n_cells*n_cells;
           M = (n_cells).^3;
           
           x = r(:,1);
           y = r(:,2);
           z = r(:,3);
           
           cell_index_x = floor(x .*n_per_length)+1;
           cell_index_y = floor(y.*n_per_length)+1;
           cell_index_z = floor(z .*n_per_length)+1;
           
           cell_index =  1 + cell_index_x + n_cells.*cell_index_y + ...
                                   n_cells_z_level.*cell_index_z;
                              
             
             
              % Reflect particles across boundaries.
              %should be done more intelligently.
                               
              [pos_x2, ~, cell_index_temp] = find(cell_index_x==n_cells-2);
              cell_index_temp = zeros(length(cell_index_temp),1);
              cell_index_bdry_pts = [1+ cell_index_temp + n_cells.*cell_index_y(pos_x2) + ...
                                                      n_cells_z_level.*cell_index_z(pos_x2)];
                                                  
               isel_x1 = cell_index_x~= n_cells -2;
                                  
              [pos_y2, ~, cell_index_temp] = find(cell_index_y==n_cells-2 & isel_x1 ...
                                                                    );
              cell_index_temp = zeros(length(cell_index_temp),1);
              cell_index_bdry_pts = [cell_index_bdry_pts; ...
                                                      1+ cell_index_x(pos_y2) + n_cells.*cell_index_temp + ...
                                                      n_cells_z_level.*cell_index_z(pos_y2)];
                   
              isel_y = cell_index_y~= n_cells-2 & isel_x1 ;
              
              [pos_z1, ~, cell_index_temp] = find(cell_index_z==1 & isel_y);
              cell_index_temp = (n_cells-1).*ones(length(cell_index_temp),1);
              cell_index_bdry_pts = [cell_index_bdry_pts; ...
                                                     1+ cell_index_x(pos_z1) + n_cells.*cell_index_y(pos_z1) + ...
                                                      n_cells_z_level.*cell_index_temp];
                               
                                                  
             cells = 1:M;
             z_level_cell = floor(cells/n_cells_z_level);
             y_level_cell  = floor((cells - z_level_cell.*n_cells_z_level)/ n_cells);
             x_level_cell = mod(cells - z_level_cell*n_cells_z_level - y_level_cell*n_cells, ...
                                       n_cells);
                           
                                   
             %cells not on boundary
             i_select = x_level_cell >1 &  ...
                               y_level_cell > 0 & ...
                               y_level_cell < n_cells -1 & ...
                               z_level_cell > 0 & ...
                               z_level_cell < n_cells - 1;
                           
              cells = cells(i_select)';
             
             
              bdry_pts = [pos_x2; pos_y2; ...
                                    pos_z1];
              %Construct cell neighbor list
             %14 = round(27/2) neighbors per cell
             %To avoid double counting.
             
             cells_left = cells - 1;
             cells_up = cells + n_cells_z_level;
             cells_front = cells - n_cells;
             cells_up_front = cells_up - n_cells;
             cells_up_back = cells_up + n_cells;
             
             neigh_cells = [cells, cells_left, cells_front, ...
                                    cells_front-1 , cells_front+1, ...
                                    cells_up, cells_up+1, cells_up-1, ...
                                    cells_up_front, cells_up_front+1, ...
                                    cells_up_front-1, ...
                                    cells_up_back, cells_up_back+1, ...
                                    cells_up_back-1];
           
                %Unavoidable for loop  AFAIK
                u = zeros(N,1);
                f = zeros(N,3);
               
                for i = 1:size(cells,1)
                        cell_i = cells(i);
                        [sources,~] = find(cell_index==cell_i);
                        n_sources = length(sources);
                        %Interactions within same cell
                        %Minimum image convention automatically followed.
                        target_cells = neigh_cells(i, :);
                   
                        %Interactions within box.
                        targets = sources;
                        n_targets = n_sources;
                        xij =  triu(repmat(x(sources),1, n_targets)) - ...
                                  triu(repmat(x(targets)',n_sources,1));
                        
                        yij =  triu(repmat(y(sources),1, n_targets)) - ...
                                  triu(repmat(y(targets)',n_sources,1));
                             
                         zij =  triu(repmat(z(sources),1, n_targets)) - ...
                                  triu(repmat(z(targets)',n_sources,1));
                          %Minimum image convention
                         xij = xij - Ls*round(xij/Ls);
                         yij = yij - Ls*round(yij/Ls);
                         zij = zij - Ls*round(zij/Ls);
   
                        rij2 = xij.^2 + yij.^2 + zij.^2;
                        isel = rij2>rc2;
                        rij2(isel) = 0.e0;
                        [k1, k2, rij2] = find(rij2);
                        xij(isel) = 0.e0;
                        yij(isel) = 0.e0;
                        zij(isel) = 0.e0;
   
                        %The 6 term
                        rij6_inv = (1./rij2).*(1./rij2).*(1./rij2);
                        %The 12 term
                        rij12_inv =rij6_inv.*rij6_inv;
                        
                        %Force on i due to j , along rij pointing towards i.
                        magfij = 24.e0.*(2.e0*rij12_inv - rij6_inv)./rij2;
                        magfij = sparse(k1,k2,magfij,n_sources,n_targets);
                        fij_x = magfij.*xij;
                        fij_y = magfij.*yij;
                        fij_z = magfij.*zij;
                        fij_x = fij_x - fij_x';
                        fij_y = fij_y - fij_y';
                        fij_z =  fij_z - fij_z';
                        
                        f(sources,1) = f(sources,1) + sum(fij_x,2);
                        f(sources,2) = f(sources,2) + sum(fij_y,2);
                        f(sources,3) = f(sources,3) + sum(fij_z,2);
                      
                        
                        uij = 4.e0.*(rij12_inv - rij6_inv);
                         uij = sparse(k1,k2,uij,n_sources,n_targets);
                         uij = uij + uij';
                        u(sources) = u(sources) + sum(uij,2);
                     
                        %Interactions with particles in other neighbouring cells.
                        targets = [];
                       for k = 2:size(neigh_cells,2)
                            [tars,~] = find(cell_index==target_cells(k));
                            [moretars,~] = find(cell_index_bdry_pts==target_cells(k));
                            moretars = bdry_pts(moretars);
                            targets = [targets; tars; moretars];
                        end
                      
                        n_targets = length(targets);
                        
                        xij =  repmat(x(sources),1, n_targets) - ...
                                  repmat(x(targets)',n_sources,1);
                        
                        yij =  repmat(y(sources),1, n_targets) - ...
                                  repmat(y(targets)',n_sources,1);
                             
                         zij =  repmat(z(sources),1, n_targets) - ...
                                  repmat(z(targets)',n_sources,1);
                       
                        %Minimum image convention
                         xij = xij - Ls*round(xij/Ls);
                         yij = yij - Ls*round(yij/Ls);
                         zij = zij - Ls*round(zij/Ls);
   
                        rij2 = xij.^2 + yij.^2 + zij.^2;
                        isel = rij2>rc2;
                        rij2(isel) = 0.e0;
                        [k1, k2, rij2] = find(rij2);
                        xij(isel) = 0.e0;
                        yij(isel) = 0.e0;
                        zij(isel) = 0.e0;
   
                        %The 6 term
                        rij6_inv = (1./rij2).*(1./rij2).*(1./rij2);
                        %The 12 term
                        rij12_inv =rij6_inv.*rij6_inv;
    
                        %Force on i due to j , along rij pointing towards i.
                        magfij = 24.e0.*(2.e0*rij12_inv - rij6_inv)./rij2;
                        magfij = sparse(k1,k2,magfij,n_sources,n_targets);
                        fij_x = magfij.*xij;
                        fij_y = magfij.*yij;
                        fij_z = magfij.*zij;
                        
                        f(sources,1) = f(sources,1) + sum(fij_x,2);
                        f(targets, 1) = f(targets,1) - sum(fij_x,1)';
                       
                        f(sources,2) = f(sources,2) + sum(fij_y,2);
                        f(targets, 2) = f(targets,2) - sum(fij_y,1)';
                        
                        f(sources,3) = f(sources,3) + sum(fij_z,2);
                        f(targets, 3) = f(targets,3) - sum(fij_z,1)';
                        
                        uij = 4.e0.*(rij12_inv - rij6_inv);
                         uij = sparse(k1,k2,uij,n_sources,n_targets);
                        u(sources) = u(sources) + sum(uij,2);
                        u(targets)  = u(targets)  + sum(uij,1)';
                end
                % disp(count);
            
end            
            

















%end
