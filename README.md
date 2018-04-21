## Unscented Kalman Filter

[//]: # (Image References)

[dataset1]: ./etc/dataset1_unscented.png "Original dataset"
[lidar_only]: ./etc/lidar_only_unscented.png "Lidar only"
[radar_only]: ./etc/radar_only_unscented.png "Radar only"
[lidar_nis]: ./etc/lidar_nis.png "Lidar NIS"
[radar_nis]: ./etc/radar_nis.png "Radar NIS"

In this project, we use an unscented kalman filter to predict the location and velocity of a simulated bicycle that is traveling around the vehicle. The measurement data comes from lidar and radar sensors with the main algorithm implemented in C++.

### Example of Predicted location
In the screenshots below, the green triangles represent the predicted location, the red circles are from the laser sensor, and the blue markers are from the radar sensor. We measure the accuracy of the algorithm by calculating the RMSE of the `x, y` positions and the velocity along the `x, y` axis. 

![alt text][dataset1]

* The original dataset starting with lidar measurement 

If we just use one or the other of the sensor measurements to update the algorithm we can start to see that when taken separately the filter performs worse. This makese sense because with both types of sensors, which are each good at a particular type of sensing, we get more information with which to update our understanding.

![alt text][lidar_only]

* The original dataset starting with lidar measurement and only using the lidar measurements to update to algorithm. We can see that compared to using both sources of sensor data the overall algorithm performs worse, especially with regard to the velocity predictions. 

![alt text][radar_only]

* The original dataset starting with lidar measurement and only using the radar measurements to update to algorithm. Compared to using only the lidar data, the radar only updated algorithm is worse at localizing the positon (higher RMSE for `x` and `y`).

### Normalized Innovation Squared (NIS)
To check the consistency of our filter and more specifically to tune our noise parameters, we measured the NIS. The NIS metric gives an approximate idea of whether or not the parameters were initialized in the correct range. Specifically, the NIS follows a Chi-squared distribution and given the number of dimensions we can figure out what the value should be if we expect that only in 5% of the cases the NIS will exceed the value. In our case, for the radar measurement example we have 3 dimensions and the value that we expect to exceed 5% of the time is `7.815`. For the lidar, since there are 2 dimensions (the `x` and `y` position) the value we expect to exceed ~5% of the time is `5.991`.

![alt text][lidar_nis]

* For the lidar NIS, we see that only in about ~5% of the cases does the value of the NIS exceed the `5.991` value.

![alt text][radar_nis]

* Similarly, for the radar NIS, about ~5% of the time the NIS exceeds `7.815`

### Compile and Build
In order to compile and build this project, make sure that the following dependencies are met.

* `cmake`:
  * For Mac make sure that `cmake` is at least version 3.5
* `make`:
  * For Mac make sure that `make` is at least version 4.1
* `gcc/g++`:
  * For Mac make sure that `gcc/g++` is at least version 5.4
* `uWebSocketIO`
  * From the project directory run `install-mac.sh`, which should be linked to the necessary `cmakepatch.txt` file
  * In order to run the above shell script, `homebrew` should also be installed

Once the above dependencies are installed:

1. Clone this repository
2. Create a build directory and navigate into it
  * `mkdir build && cd build`
3. Compile 
  * `cmake .. && make`
4. Run the program
  * Make sure that the [Udacity simulator](https://github.com/udacity/self-driving-car-sim/releases) is installed
  * `./UnscentedKF`


