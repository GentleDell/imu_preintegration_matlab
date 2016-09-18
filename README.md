#PreIntegration Method for the fusion of IMU data with GPS

I am amazed at the optimization based method for sensor fusion. So I do a tiny test to fuse one time stamp GPS data with IMU output. The referrence is [IMU Preintegration on Manifold for Efficient Visual-Inertial Maximum-a-Posteriori Estimation](http://www.roboticsproceedings.org/rss11/p06.pdf). 

The matlab codes here can also differentiate a sequence of SE3 poses into IMU data using interpolation from the referrence: [Spline Fusion: A continuous-time representation for visual-inertial fusion with application to rolling shutter cameras](https://arpg.colorado.edu/wp-content/uploads/Spline-fusion.pdf).

