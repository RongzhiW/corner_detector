#include "imageContainer.h"

imageContainer::imageContainer()
{}

void imageContainer::setup(){

    const char* imgTopicName="/ardrone/bottom/image_raw";
    imageSub=nh.subscribe(imgTopicName, 1, &imageContainer::imgCb,this);
}

void imageContainer::imgCb(const sensor_msgs::Image &msg){
    cv_bridge::CvImageConstPtr cv_img
                = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::BGR8);
    cv::Mat imgTmp=cv_img->image;
    imgTmp.copyTo(curImage);

    //cv::imshow("camera view", cv_img->image);
    //cv::waitKey(10);
    //std::cout << "ardrone camera data are read!" << std::endl;
}

void imageContainer::startImageRev(){
    pthread_create(&pthreadSpinId,NULL,&startSpinOnceThread,this);
}

static void* startSpinOnceThread(void* arg){
    ros::Rate loop_rate(20);
    while(ros::ok()){
        ros::spinOnce();
        //std::cout<<"spinOnce"<<std::endl;
        loop_rate.sleep();
    }
}

void imageContainer::getCurImage(cv::Mat &imgDest){

    curImage.copyTo(imgDest);
}






