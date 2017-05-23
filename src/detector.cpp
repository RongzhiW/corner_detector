#include "detector.h"


#define pi 3.1415927

Detector::Detector()
{}
Detector::~Detector()
{}

Detector::Detector(cv::Mat img):
    sobelKernel(3,3,CV_32F),sobelKernelTrs(3,3,CV_32F)
{
    if(img.data){
       img.copyTo(image);
       img.copyTo(imageRaw);
    }

    sobelKernel.col(0).setTo(cv::Scalar(-1));
    sobelKernel.col(1).setTo(cv::Scalar(0));
    sobelKernel.col(2).setTo(cv::Scalar(1));

    sobelKernelTrs=sobelKernel.t();

    float* tempRad=new float[3]{4,8,12};
    float* tempTemp=new float[18]{0,pi/2,tempRad[0],
                      pi/4,-pi/4,tempRad[0],
                      0,pi/2,tempRad[1],
                      pi/4,-pi/4,tempRad[1],
                      0,pi/2,tempRad[2],
                      pi/4,-pi/4,tempRad[2]};
    radius=cv::Mat(3,1,CV_32F,tempRad);
    templateProps=cv::Mat(6,3,CV_32F,tempTemp);

    //std::cout<<templateProps<<std::endl;

}

void Detector::getImageAngle(cv::Mat imDv,cv::Mat imDu,cv::Mat &result){
    int rowsLeft=imDv.rows;
    int colsLeft=imDv.cols;
    int rowsRight=imDu.rows;
    int colsRight=imDu.cols;
    if(rowsLeft!=rowsRight||colsLeft!=colsRight)return;

    int channels=imDv.channels();
    result.create(rowsLeft,colsLeft,CV_32F);

    int nr=rowsLeft;
    int nc=colsLeft;
    if(imDv.isContinuous()){
        nc=nc*nr;
        nr=1;
        //std::cout<<"continue"<<std::endl;
    }
    for(int i=0;i<nr;i++){
        const float* dataLeft=imDv.ptr<float>(i);
        const float* dataRight=imDu.ptr<float>(i);
        float* dataResult=result.ptr<float>(i);
        for(int j=0;j<nc*channels;++j){
            if(abs(dataRight[j])>=0.000001){
                dataResult[j]=atan2((float)dataLeft[j],(float)dataRight[j]);
                if(dataResult[j]<0)dataResult[j]+=pi;
                else if(dataResult[j]>pi)dataResult[j]-=pi;
                //std::cout<<(float)dataLeft[j]<<' '<<(float)dat/stl_algobase.h:198:7: error: could not convert 'cv::operator<(caRight[j]<<' '<<dataResult[j]<<std::endl;
            }


        }
    }
}

void Detector::getImageWeight(cv::Mat imDv,cv::Mat imDu,cv::Mat &result){
    int rowsLeft=imDv.rows;
    int colsLeft=imDv.cols;
    int rowsRight=imDu.rows;
    int colsRight=imDu.cols;
    if(rowsLeft!=rowsRight||colsLeft!=colsRight)return;

    int channels=imDv.channels();
    result.create(rowsLeft,colsLeft,CV_32F);

    int nr=rowsLeft;
    int nc=colsLeft;
    if(imDv.isContinuous()){
        nc=nc*nr;
        nr=1;
        //std::cout<<"continue"<<std::endl;
    }
    for(int i=0;i<nr;i++){
        const float* dataLeft=imDv.ptr<float>(i);
        const float* dataRight=imDu.ptr<float>(i);
        float* dataResult=result.ptr<float>(i);
        for(int j=0;j<nc*channels;++j){
            dataResult[j]=sqrt((float)dataLeft[j]*(float)dataLeft[j]+(float)dataRight[j]*(float)dataRight[j]);
            //std::cout<<(float)dataLeft[j]<<' '<<(float)dataRight[j]<<' '<<dataResult[j]<<std::endl;

        }
    }

}

float Detector::normpdf(float dist,float mu,float sigma){
    return exp(-0.5*(dist-mu)*(dist-mu)/(sigma*sigma))/(sqrt(2*pi)*sigma);
}

void Detector::createCorrelationPatch(const cv::Mat& templateProp){
    float angle_1=templateProp.at<float>(0);
    float angle_2=templateProp.at<float>(1);
    float rad=templateProp.at<float>(2);
    int width=(int)rad*2+1;
    int height=(int)rad*2+1;
    templateA1=cv::Mat::zeros(height,width,CV_32F);
    templateA2=cv::Mat::zeros(height,width,CV_32F);
    templateB1=cv::Mat::zeros(height,width,CV_32F);
    templateB2=cv::Mat::zeros(height,width,CV_32F);

    int mu=rad;
    int mv=rad;
    float nor1[]={-sin(angle_1),cos(angle_1)};
    float nor2[]={-sin(angle_2),cos(angle_2)};

    for(int u=0;u<width;++u){
        for(int v=0;v<height;++v){
            float vec[]={u-mu,v-mv};
            float dis=sqrt((u-mu)*(u-mu)+(v-mv)*(v-mv));
            float side1=vec[0]*nor1[0]+vec[1]*nor1[1];
            float side2=vec[0]*nor2[0]+vec[1]*nor2[1];
            if(side1<=-0.1&&side2<=-0.1){
                templateA1.at<float>(v,u)=normpdf(dis,0,rad/2);
            }
            else if(side1>=0.1&&side2>=0.1){
                templateA2.at<float>(v,u)=normpdf(dis,0,rad/2);
            }
            else if(side1<=-0.1&&side2>=0.1){
                templateB1.at<float>(v,u)=normpdf(dis,0,rad/2);
            }
            else if(side1>=0.1&&side2<=-0.1){
                templateB2.at<float>(v,u)=normpdf(dis,0,rad/2);
            }
        }
    }
    templateA1=templateA1/cv::sum(templateA1)[0];
    templateA2=templateA2/cv::sum(templateA2)[0];
    templateB1=templateB1/cv::sum(templateB1)[0];
    templateB2=templateB2/cv::sum(templateB2)[0];

    //std::cout<<templateA1<<' '<<templateA2<<' '<<templateB1<<' '<<templateB2<<std::endl;
    //std::cout<<templateA1<<std::endl;

}

void Detector::getMin(const cv::Mat& src1,const cv::Mat& src2,cv::Mat& dst){
    int rowsLeft=src1.rows;
    int colsLeft=src1.cols;
    int rowsRight=src2.rows;
    int colsRight=src2.cols;
    if(rowsLeft!=rowsRight||colsLeft!=colsRight)return;

    int channels=src1.channels();

    int nr=rowsLeft;
    int nc=colsLeft;
    if(src1.isContinuous()){
        nc=nc*nr;
        nr=1;
        //std::cout<<"continue"<<std::endl;
    }
    for(int i=0;i<nr;i++){
        const float* dataLeft=src1.ptr<float>(i);
        const float* dataRight=src2.ptr<float>(i);
        float* dataResult=dst.ptr<float>(i);
        for(int j=0;j<nc*channels;++j){
            dataResult[j]=(dataLeft[j]<dataRight[j])?dataLeft[j]:dataRight[j];
        }
    }
}

void Detector::getMax(const cv::Mat& src1,const cv::Mat& src2,cv::Mat& dst){
    int rowsLeft=src1.rows;
    int colsLeft=src1.cols;
    int rowsRight=src2.rows;
    int colsRight=src2.cols;
    if(rowsLeft!=rowsRight||colsLeft!=colsRight)return;

    int channels=src1.channels();

    int nr=rowsLeft;
    int nc=colsLeft;
    if(src1.isContinuous()){
        nc=nc*nr;
        nr=1;
        //std::cout<<"continue"<<std::endl;
    }
    for(int i=0;i<nr;i++){
        const float* dataLeft=src1.ptr<float>(i);
        const float* dataRight=src2.ptr<float>(i);
        float* dataResult=dst.ptr<float>(i);
        for(int j=0;j<nc*channels;++j){
            dataResult[j]=(dataLeft[j]>=dataRight[j])?dataLeft[j]:dataRight[j];
        }
    }
}
void Detector::nonMaximumSuppression(const cv::Mat& img,int n,float tau,int margin){
    int width=img.cols;
    int height=img.rows;
    cornerPoints.clear();
    if(width==0||height==0)return;
    int maxi=0;
    int maxj=0;
    float maxVal=img.at<float>(maxj,maxi);
    for(int i=n+margin;i<width-n-margin;i+=(n+1)){
        for(int j=n+margin;j<height-n-margin;j+=(n+1)){
            maxi=i;
            maxj=j;
            maxVal=img.at<float>(j,i);
            for(int i2=i;i2<i+n+1;++i2){
                for(int j2=j;j2<j+n+1;++j2){
                    float val=img.at<float>(j2,i2);
                    if(val>maxVal){
                        maxi=i2;
                        maxj=j2;
                        maxVal=val;
                    }
                }
            }
            bool failed=false;
            for(int i3=maxi-n;i3<std::min(maxi+n,width-margin);++i3){
                for(int j3=maxj-n;j3<std::min(maxj+n,height-margin);++j3){
                    float val=img.at<float>(j3,i3);
                    if(val>maxVal && (i3<i||i3>i+n||j3<j||j3>j+n)){
                        failed=true;
                        break;
                    }
                }
                if(failed)break;
            }
            if(maxVal>=tau&& !failed){
                cv::Point* oneCorner=new cv::Point;
                oneCorner->x=maxi;
                oneCorner->y=maxj;
                cornerPoints.push_back(oneCorner);
            }
        }
    }
}

void Detector::cornerDetect(){
    if(!image.data){
        std::cerr<<"no image data!"<<std::endl;
        return;
    }

    cv::cvtColor(image,image,CV_BGR2GRAY);

    cv::Mat imageDu(image.size(),CV_32F);
    cv::Mat imageDv(image.size(),CV_32F);
    //std::cout<<"image_type:"<<image.type()<<std::endl;

    cv::Mat imageAngle(image.size(),CV_32F);
    cv::Mat imageWeight(image.size(),CV_32F);
    cv::filter2D(image,imageDu,CV_32F,sobelKernel);
    cv::filter2D(image,imageDv,CV_32F,sobelKernelTrs);

    getImageAngle(imageDv,imageDu,imageAngle);
    getImageWeight(imageDv,imageDu,imageWeight);

    cv::Mat imageNorm(image.size(),CV_32F);
    //imageNorm=image;
    cv::normalize(image,imageNorm,0,1,cv::NORM_MINMAX,CV_32F);

    int rowTempPros=templateProps.rows;
    cv::Mat imgCorners(image.size(),CV_32F);
    imgCorners.zeros(image.size(),CV_32F);
    cv::Mat imgCornerA1(image.size(),CV_32F);
    cv::Mat imgCornerA2(image.size(),CV_32F);
    cv::Mat imgCornerB1(image.size(),CV_32F);
    cv::Mat imgCornerB2(image.size(),CV_32F);
    cv::Mat imgCornerMean(image.size(),CV_32F);
    cv::Mat imgCornerA(image.size(),CV_32F);
    cv::Mat imgCornerB(image.size(),CV_32F);
    cv::Mat imgCorner1(image.size(),CV_32F);
    cv::Mat imgCorner2(image.size(),CV_32F);
    for(int i=0;i<rowTempPros;++i){

        createCorrelationPatch(templateProps.row(i));
        cv::filter2D(imageNorm,imgCornerA1,CV_32F,templateA1);
        cv::filter2D(imageNorm,imgCornerA2,CV_32F,templateA2);
        cv::filter2D(imageNorm,imgCornerB1,CV_32F,templateB1);
        cv::filter2D(imageNorm,imgCornerB2,CV_32F,templateB2);

        imgCornerMean=(imgCornerA1+imgCornerA2+imgCornerB1+imgCornerB2)/4;
        getMin(imgCornerA1-imgCornerMean,imgCornerA2-imgCornerMean,imgCornerA);
        getMin(imgCornerMean-imgCornerB1,imgCornerMean-imgCornerB2,imgCornerB);
        getMin(imgCornerA,imgCornerB,imgCorner1);

        getMin(imgCornerMean-imgCornerA1,imgCornerMean-imgCornerA2,imgCornerA);
        getMin(imgCornerB1-imgCornerMean,imgCornerB2-imgCornerMean,imgCornerB);
        getMin(imgCornerA,imgCornerB,imgCorner2);

        getMax(imgCorners,imgCorner1,imgCorners);
        getMax(imgCorners,imgCorner2,imgCorners);

//        nonMaximumSuppression(imgCorners,3,0.05,5);

//        int cornerSize=cornerPoints.size();
//        //std::cout<<"size:"<<cornerSize<<std::endl;
//        for(int i=0;i<cornerSize;++i){
//            cv::Point* p=cornerPoints[i];
//            cv::circle(imageRaw,*p,5,CV_RGB(255,0,0),3,8,0);
//        }

//        cv::namedWindow("corner detector result");
//        cv::imshow("corner detector result",imageRaw);
//        cv::waitKey(0);

    }

    nonMaximumSuppression(imgCorners,6,0.1,10);



    int cornerSize=cornerPoints.size();
    std::cout<<"size:"<<cornerSize<<std::endl;
    for(int i=0;i<cornerSize;++i){
        cv::Point* p=cornerPoints[i];
        cv::circle(imageRaw,*p,5,CV_RGB(255,0,0),3,8,0);
//        std::cout<<"u:"<<t->x<<" v:"<<t->y<<std::endl;
    }




    //std::cout<<"image_type:"<<image_du.type()<<std::endl;
    //std::cout<<"image_mat:"<<imgCorners<<std::endl;
    //std::cout<<templateProps<<' '<<radius<<std::endl;

    cv::namedWindow("corner detector result");
    cv::imshow("corner detector result",imageRaw);
//    cv::namedWindow("corner detector result1");
//    cv::imshow("corner detector result1",imgCorners);

    cv::waitKey(0);

}


