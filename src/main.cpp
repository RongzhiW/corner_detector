#include "iostream"
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "dirent.h"
#include "unistd.h"
#include "detector.h"

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
    string basePath="/home/tianu/onlyCatkin_ws/src/corner_detector/data/";
    cv::Mat image=cv::imread(basePath+"caliPattern_1.jpg");

    if(!image.data){
        cerr<<"image load failed"<<endl;
        exit(-1);
    }

    cv::namedWindow("original image");
    cv::imshow("original image",image);
    cv::waitKey(0);

    Detector corner_detector(image);
    corner_detector.cornerDetect();



    return 0;
}
