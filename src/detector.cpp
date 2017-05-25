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
            if(std::abs(dataRight[j])>=0.000001){
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
            dataResult[j]=std::sqrt((float)dataLeft[j]*(float)dataLeft[j]+(float)dataRight[j]*(float)dataRight[j]);
            //std::cout<<(float)dataLeft[j]<<' '<<(float)dataRight[j]<<' '<<dataResult[j]<<std::endl;

        }
    }

}

float Detector::normpdf(float dist,float mu,float sigma){
    return exp(-0.5*(dist-mu)*(dist-mu)/(sigma*sigma))/(std::sqrt(2*pi)*sigma);
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
            float dis=std::sqrt((u-mu)*(u-mu)+(v-mv)*(v-mv));
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
    cornersEdge1.clear();
    cornersEdge2.clear();
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
                std::vector<float> e1(2,0.0);
                std::vector<float> e2(2,0.0);
                cornersEdge1.push_back(e1);
                cornersEdge2.push_back(e2);
            }
        }
    }
}

void Detector::refineCorners(const cv::Mat& imgDu,const cv::Mat& imgDv,const cv::Mat& imgAngle,const cv::Mat& imgWeight,int r){
    int width=imgDu.cols;
    int height=imgDv.rows;
    int cornersSize=cornerPoints.size();
    for(int i=0;i<cornersSize;++i){
        cv::Point* pTmp=cornerPoints[i];
        int cu=pTmp->x;
        int cv=pTmp->y;
        int left=std::max(cu-r,0);
        int right=std::min(cu+r,width-1);
        int top=std::max(cv-r,0);
        int down=std::min(cv+r,height-1);
        edgeOrientations(imgAngle,imgWeight,left,right,top,down,i);

        if((cornersEdge1[i][0]==0.0&&cornersEdge1[i][1]==0.0) || (cornersEdge2[i][0]==0.0&&cornersEdge2[i][1]==0.0))continue;
        //static int counter=0;
        //std::cout<<"count:"<<++counter<<" ang1:"<<cornersEdge1[i][0]<<' '<<cornersEdge1[i][1\
                       ]<<" ang2:"<<cornersEdge2[i][0]<<' '<<cornersEdge2[i][1]<<std::endl;


        cv::Mat A1=cv::Mat::zeros(2,2,CV_32F);
        cv::Mat A2=cv::Mat::zeros(2,2,CV_32F);
        //std::cout<<"A1:"<<A1<<std::endl;
        //std::cout<<"A2:"<<A2<<std::endl;

        cv::Mat G=cv::Mat::zeros(2,2,CV_32F);
        cv::Mat b=cv::Mat::zeros(2,1,CV_32F);
        //std::cout<<"G:"<<G<<std::endl;
        //std::cout<<"b:"<<b<<std::endl;

        std::vector<float> e1(cornersEdge1[i]);
        std::vector<float> e2(cornersEdge2[i]);
        //std::cout<<"e1:"<<e1[0]<<' '<<e1[1]<<" e2:"<<e2[0]<<' '<<e2[1]<<std::endl;

        for(int u=left;u<=right;++u){
            for(int v=top;v<=down;v++){
                std::vector<float> gp;
                gp.push_back(imgDu.at<float>(v,u));
                gp.push_back(imgDv.at<float>(v,u));
                float normLen=std::sqrt(gp[0]*gp[0]+gp[1]*gp[1]);
                //std::cout<<"normLen:"<<normLen<<std::endl;
                if(normLen<0.1)continue;
                gp[0]=gp[0]/normLen;
                gp[1]=gp[1]/normLen;


//                float dot_gp_e1=gp[0]*e1[0]+gp[1]*e1[1];
//                dot_gp_e1=(dot_gp_e1>0)?dot_gp_e1:-dot_gp_e1;
                float dot_gp_e1=std::abs(gp[0]*e1[0]+gp[1]*e1[1]);
                if(dot_gp_e1<0.25){
                    A1.at<float>(0,0)+=gp[0]*gp[0];
                    A1.at<float>(0,1)+=gp[0]*gp[1];
                    A1.at<float>(1,0)+=gp[1]*gp[0];
                    A1.at<float>(1,1)+=gp[1]*gp[1];
                }

//                float dot_gp_e2=gp[0]*e2[0]+gp[1]*e2[1];
//                dot_gp_e2=(dot_gp_e2>0)?dot_gp_e2:-dot_gp_e2;
                float dot_gp_e2=std::abs(gp[0]*e2[0]+gp[1]*e2[1]);
                if(dot_gp_e2<0.25){
                    A2.at<float>(0,0)+=gp[0]*gp[0];
                    A2.at<float>(0,1)+=gp[0]*gp[1];
                    A2.at<float>(1,0)+=gp[1]*gp[0];
                    A2.at<float>(1,1)+=gp[1]*gp[1];
                }
                //std::cout<<"A1:"<<A1<<" A2:"<<A2<<std::endl;
                //std::cout<<"dot_gp_e1:"<<dot_gp_e1<<" dot_gp_e2:"<<dot_gp_e2<<std::endl;
                std::vector<float> vec_cp;
                vec_cp.push_back(u-cu);
                vec_cp.push_back(v-cv);
                if(u!=cu||v!=cv){
                    float dot=vec_cp[0]*e1[0]+vec_cp[1]*e1[1];
                    float d1u=vec_cp[0]-dot*e1[0];
                    float d1v=vec_cp[1]-dot*e1[1];
                    float d1=std::sqrt(d1u*d1u+d1v*d1v);

                    float dot2=vec_cp[0]*e2[0]+vec_cp[1]*e2[1];
                    float d2u=vec_cp[0]-dot2*e2[0];
                    float d2v=vec_cp[1]-dot2*e2[1];
                    float d2=std::sqrt(d2u*d2u+d2v*d2v);

                    if((d1<3&&dot_gp_e1<0.25) || (d2<3&&dot_gp_e2<0.25)){

                        G.at<float>(0,0)+=gp[0]*gp[0];
                        G.at<float>(0,1)+=gp[0]*gp[1];
                        G.at<float>(1,0)+=gp[1]*gp[0];
                        G.at<float>(1,1)+=gp[1]*gp[1];
                        b.at<float>(0,0)+=gp[0]*gp[0]*u+gp[0]*gp[1]*v;
                        b.at<float>(1,0)+=gp[1]*gp[0]*u+gp[1]*gp[1]*v;
                        //std::cout<<"G:"<<G<<std::endl;
                        //std::cout<<"b:"<<b<<std::endl;
                        //std::cout<<gp[0]*gp[0]<<' '<<gp[0]*gp[1]<<' '<<gp[1]*gp[0]<<' '\
                                   <<gp[1]*gp[1]<<' '<<gp[0]*gp[0]*u+gp[0]*gp[1]*v<<' '<<gp[1]*gp[0]*u+gp[1]*gp[1]*v<<std::endl;

                    }
                    //std::cout<<"G:"<<G<<" b:"<<b<<std::endl;
                }

            }
        }
        cv::Mat u,w,vt;
        cv::SVDecomp(A1,w,u,vt,cv::SVD::MODIFY_A | cv::SVD::FULL_UV);
        cornersEdge1[i][0]=vt.at<float>(1,0);
        cornersEdge1[i][1]=vt.at<float>(1,1);

        cv::Mat u2,w2,vt2;
        cv::SVDecomp(A2,w2,u2,vt2,cv::SVD::MODIFY_A | cv::SVD::FULL_UV);
        cornersEdge2[i][0]=vt2.at<float>(1,0);
        cornersEdge2[i][1]=vt2.at<float>(1,1);

        //std::cout<<"newang1:"<<cornersEdge1[i][0]<<' '<<cornersEdge1[i][1\
                       ]<<" newang2:"<<cornersEdge2[i][0]<<' '<<cornersEdge2[i][1]<<std::endl;
//        float absDeterG=cv::determinant(G);
//        absDeterG=(absDeterG>0)?absDeterG:-absDeterG;
        //std::cout<<"G:"<<G<<std::endl;
        //std::cout<<"det(G):"<<std::abs(cv::determinant(G))<<std::endl;
        if(std::abs(cv::determinant(G))>0.001){
            cv::Mat cornerNew(2,1,CV_32F);
            cornerNew=G.inv()*b;

            cv::Point* pTmp=cornerPoints[i];
            float new_u=cornerNew.at<float>(0,0);
            float new_v=cornerNew.at<float>(1,0);
            float cur_u=pTmp->x;
            float cur_v=pTmp->y;
            //std::cout<<"new_u:"<<new_u<<" cur_u:"<<cur_u<<" new_v:"<<new_v<<" cur_v:"<<cur_v<<std::endl;
//            if(sqrt((new_u-cur_u)*(new_u-cur_u)+(new_v-cur_v)*(new_v-cur_v))>=4){
//                cornersEdge1[i][0]=0.0;
//                cornersEdge1[i][1]=0.0;
//                cornersEdge2[i][0]=0.0;
//                cornersEdge2[i][1]=0.0;

//            }
//            pTmp->x=new_u;
//            pTmp->x=new_v;

        }
        else{
            cornersEdge1[i][0]=0.0;
            cornersEdge1[i][1]=0.0;
            cornersEdge2[i][0]=0.0;
            cornersEdge2[i][1]=0.0;
        }


    }
}

void Detector::edgeOrientations(const cv::Mat& imgAngle,const cv::Mat& imgWeight,int left,int right,int top,int down,int index){
    std::vector<float> vecAngle;
    std::vector<float> vecWeight;
    for(int u=left;u<=right;++u){
        for(int v=top;v<=down;++v){
            float ang=imgAngle.at<float>(v,u)+pi/2;
            ang=ang>pi?(ang-pi):ang;
            vecAngle.push_back(ang);
            vecWeight.push_back(imgWeight.at<float>(v,u));
        }
    }
    int binNum=32;
    int sizeVec=vecAngle.size();
    std::vector<float> angleHist(binNum,0);
    for(int i=0;i<sizeVec;++i){
        int bin=std::max(std::min((int)floor(vecAngle[i]/(pi/binNum)),binNum-1),0);
        angleHist[bin]+=vecWeight[i];
    }
    std::vector<float> angleHistSmoothed(angleHist);
    std::vector<std::pair<float,int> > modes;
    findModesMeanShift(angleHist,angleHistSmoothed,modes,1.0);
    int sizeModes=modes.size();
    if(sizeModes<2)return;
    std::pair<float,int> most1=modes[sizeModes-1];
    std::pair<float,int> most2=modes[sizeModes-2];
    float most1Angle=most1.second*pi/binNum;
    float most2Angle=most2.second*pi/binNum;
    float tmp=most1Angle;
    most1Angle=(most1Angle>most2Angle)?most1Angle:most2Angle;
    most2Angle=(tmp>most2Angle)?most2Angle:tmp;
    float deltaAngle=std::min(most1Angle-most2Angle,most2Angle+(float)pi-most1Angle);
    if(deltaAngle<=1.2)return;
    cornersEdge1[index][0]=cos(most1Angle);
    cornersEdge1[index][1]=sin(most1Angle);
    cornersEdge2[index][0]=cos(most2Angle);
    cornersEdge2[index][1]=sin(most2Angle);

//    static int counter=0;
//    static float most1AngleMin=4;
//    static float most1AngleMax=0;
//    if(most1Angle-most2Angle>most1AngleMax)most1AngleMax=most1Angle-most2Angle;
//    if(most1Angle-most2Angle<most1AngleMin)most1AngleMin=most1Angle-most2Angle;
//    std::cout<<"count:"<<++counter<<" ang1:"<<most1Angle<<" ang1-ang2:"<<most1Angle-most2Angle<<" min:"<<most1AngleMin<<" max:"<<most1AngleMax<<std::endl;
    //std::cout<<"count:"<<++counter<<" ang1:"<<cornersEdge1[i][0]<<' '<<cornersEdge1[i][1\
               ]<<" ang2:"<<cornersEdge2[i][0]<<' '<<cornersEdge2[i][1]<<std::endl;




}

void Detector::findModesMeanShift(std::vector<float> &angleHist,std::vector<float> &angleHistSmoothed,std::vector<std::pair<float,int> > &modes,float sigma){
    int sizeHist=angleHist.size();
    bool allZeros=true;
    for(int i=0;i<sizeHist;++i){
        float angSmoothed=0;
        for(int j=-(int)ceil(2*sigma);j<=(int)ceil(2*sigma);++j){
            int modIndex=(i+j)<0?i+j+sizeHist:(i+j);
            angSmoothed+=angleHist[modIndex]*normpdf(j,0,sigma);
        }
        angleHistSmoothed[i]=angSmoothed;
        if(std::abs(angSmoothed-angleHistSmoothed[0])>0.0001)allZeros=false;
    }
    if(allZeros)return;
    for(int i=0;i<sizeHist;++i){
        int j=i;
        int curLeft=(j-1)<0?j-1+sizeHist:j-1;
        int curRight=(j+1)>sizeHist-1?j+1-sizeHist:j+1;
        if(angleHistSmoothed[curLeft]<angleHistSmoothed[i]&&angleHistSmoothed[curRight]<angleHistSmoothed[i]){
            modes.push_back(std::make_pair(angleHistSmoothed[i],i));
        }
    }
    std::sort(modes.begin(),modes.end());
    //std::cout<<"modes size:"<<modes.size()<<std::endl;


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
    cv::Mat imgCorners=cv::Mat::zeros(image.size(),CV_32F);
    //imgCorners.zeros(image.size(),CV_32F);
    //std::cout<<"imgCorners:"<<imgCorners<<std::endl;
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

    refineCorners(imageDu,imageDv,imageAngle,imageWeight,10);


    int cornerSize=cornerPoints.size();
    std::cout<<"size:"<<cornerSize<<std::endl;
    for(int i=0;i<cornerSize;++i){
        cv::Point* p=cornerPoints[i];
        std::vector<float> e1=cornersEdge1[i];
        std::vector<float> e2=cornersEdge2[i];
        if((e1[0]==0.0&&e1[1]==0.0)||(e2[0]==0.0&&e2[1]==0.0))continue;
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


