#ifndef DETECTOR_H
#define DETECTOR_H

#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <vector>

class Detector
{
public:
    Detector();
    Detector(cv::Mat img);

    ~Detector();

public:
    void cornerDetect();

private:
    void getImageAngle(cv::Mat imDv,cv::Mat imDu,cv::Mat &result);
    void getImageWeight(cv::Mat imDv,cv::Mat imDu,cv::Mat &result);
    void createCorrelationPatch(const cv::Mat& templateProp);
    float normpdf(float dist,float mu=0,float sigma=1);
    void getMin(const cv::Mat& src1,const cv::Mat& src2,cv::Mat& dst);
    void getMax(const cv::Mat& src1,const cv::Mat& src2,cv::Mat& dst);
    void nonMaximumSuppression(const cv::Mat& img,int n=3,float tau=0.025,int margin=5);


private:
    cv::Mat image;
    cv::Mat imageRaw;
    cv::Mat kernelA;
    cv::Mat kernelB;
    cv::Mat kernelC;
    cv::Mat kernelD;
    cv::Mat sobelKernel;
    cv::Mat sobelKernelTrs;
    cv::Mat templateProps;
    cv::Mat radius;
    cv::Mat templateA1;
    cv::Mat templateA2;
    cv::Mat templateB1;
    cv::Mat templateB2;
    std::vector<cv::Point* > cornerPoints;





};

#endif // DETECTOR_H
