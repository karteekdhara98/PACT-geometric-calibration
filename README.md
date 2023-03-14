# PACT-geometric-calibration-code

This code is part of the research paper: "A method for the geometric 
calibration of ultrasound transducer arrays with arbitrary geometries".
 
The paper is authored by: Karteekeya Sastry, Yang Zhang, Peng Hu, 
Yilin Luo, Xin Tong, Shuai Na, Lihong V. Wang. The paper is currently 
in preparation.

Author Affiliation: Caltech Optical Imaging Laboratory, Department of 
Electrical Engineering, Department of Medical Engineering, California 
Institute of Technology, 1200 East California Boulevard, Pasadena, 
CA 91125, USA.

For questions or comments about the code, please contact Karteekeya
Sastry at sdharave \at caltech.edu .

Copyright (c) 2023 Karteekeya Sastry, Lihong V. Wang
This code is licensed under the MIT license. 
See the LICENSE file for details.

Author name: Karteekeya Sastry
Caltech Optical Imaging Laboratory
Date: 14 March, 2023

Description of the code:
This is a demo of the geometric calibration procedure described in the 
above-mentioned paper. It estimates the coordinates of a simulated 
arc-shaped array using time-of-arrival (ToA) data from 27 point sources
in a 3 x 3 x 3 arrangement. The code is divided into 5 steps. 
1. Initialize the experiment parameters and load the ToA data.
2. Initialize the point source locations. 
3. Populate the A matrix and b vector as shown in Section 2 of the paper. 
4. Estimate the transducer locations.
5. Visualize the results.
