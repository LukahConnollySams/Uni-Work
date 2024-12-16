To Run Files:

Run desired file in IDE of choice

for parts (a) to (c), use game_of_life1.py
for parts (d) to (f), use game_of_life2.py

------------------------------------------------------------------------------------------
(b)For n = 2, the simulation starts with few alive cells that move across the grid
   and eventually they burst into a large number of alive cells, almost exponentially.
   they then procede to stay in a random oscillating state indefinetily (i.e. stablises).

   For n = 3, the simulation again starts with few alive cells, and instantly kills 
   all cells after the first itteration, likely due to the fact that the chance of 
   having three alive neighbours is significatnly low.

Files for part:
    alive_cells_n=2.png
    alive_cells_n=2_steady_state.png
    alive_cells_n=3.png
    alive_cells_n=3_steady_state.png

------------------------------------------------------------------------------------------
(c)The steady state behaviour of the 'box' simulation are four oscillating diagonal pairs
   of alive cells. This is because the only cells that can becom alive after the first
   itteration are the off cells orthogonal to the corners of the alive box, which create a
   diagonal pair that will oscillate with n=2.

Files for part:
    alive_cells_box.png
    alive_cells_box_steady_state.png

------------------------------------------------------------------------------------------
(e)The heatmap of the fraction of average alive cells seems to be split along a diagonal line 
   that is roughly p1 = 0.85p2 + 0.15. where everything above this line is of a clear lower 
   fraction of average alive cells than below.

Files for part:
    heatmap.dat
    heatmap.png

------------------------------------------------------------------------------------------
(f)There doesnt seem to be any clear separation in the variances of the heatmap, and the line
   produced for the maximal values, appears to be random. Equation produced from polyfit: 
   1.12x**2 - 1.18x + 0.64.  

Files for part:
    var_heatmap.dat
    var_heatmap.png
