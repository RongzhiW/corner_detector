#include "iostream"
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "dirent.h"
#include "unistd.h"
#include "detector.h"
#include "imageContainer.h"

using namespace cv;
using namespace std;

int main(int argc, char* argv[])
{
    //DIR *dir;
    //char basePath[100];
    //memset(basePath,'\0',sizeof(basePath));
    //getcwd(basePath,999);
    //cout<<"base_path:"<<basePath<<endl;

    //Load Image
    //string basePath="/home/tianu/onlyCatkin_ws/src/corner_detector/data/";
    //string imagePath=basePath+"10.jpg";
    //cv::Mat image=cv::imread(basePath+"caliPattern_1.jpg");
    //cv::Mat image=cv::imread(imagePath);

    ros::init(argc, argv, "cornerDetectorGroundThuth");
    ros::NodeHandle node;
    cv::Mat img;

    imageContainer *imgCont=new imageContainer(&node);
    imgCont->startImageRev();

    cv::namedWindow("original image");
    while(1){
        imgCont->getCurImage(img);
        if(img.data){
            cv::imshow("original image",img);

            Detector corner_detector(img);
            corner_detector.cornerDetect();
             cv::waitKey(0);

        }


    }


//    if(!image.data){
//        cerr<<"image load failed"<<endl;
//        exit(-1);
//    }

//    cv::namedWindow("original image");
//    cv::imshow("original image",image);
//    cv::waitKey(0);

//    Detector corner_detector(image);
//    corner_detector.cornerDetect();



    return 0;
}
