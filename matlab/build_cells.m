%Machinery for neighbour list calculations
%function build_cells(Ls)

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
             %cells next to boundary.
              i_select_x = x_level_cell ==2 ; 
              i_select_y = y_level_cell == 1 ;
              i_select_z = z_level_cell == 1 ;
                          
             
                                 
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
                                
                                