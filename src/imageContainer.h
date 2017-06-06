#ifndef IMAGECONTAINER_H
#define IMAGECONTAINER_H

#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <vector>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>

static void *startSpinOnceThread(void *arg);
class imageContainer
{
public:
    imageContainer();
    imageContainer(ros::NodeHandle* n){
        nh=*n;
        setup();
    }

public:
    void getCurImage(cv::Mat &imgDest);
    void startImageRev();


private:
    void setup();
    void imgCb(const sensor_msgs::Image &msg);


private:
    cv::Mat curImage;
    ros::NodeHandle nh;
    ros::Subscriber imageSub;
    pthread_t pthreadSpinId;
};

#endif // IMAGECONTAINER_H
