#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}



MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_startBtn_clicked()
{
    setWidth(640);
    setLength(480);
    setFPS(60);
    dev = ctx.get_device(0);
    // Configure depth to run at VGA resolution at 30 frames per second
    dev->enable_stream(rs::stream::depth, getWidth(), getLength(), rs::format::z16, getFPS());
    dev->start();
    ui->label->setGeometry(0,0,getWidth(),getLength());
    runBool = true;
    rs_apply_ivcam_preset((rs_device *)dev,RS_IVCAM_PRESET_BACKGROUND_SEGMENTATION);
    updateSliders();
    while(runBool){
        dev->wait_for_frames();
        depth_intrin = dev->get_stream_intrinsics(rs::stream::depth);
        depth_image = (uint16_t *)dev->get_frame_data(rs::stream::depth);
        // Create depth image
        cv::Mat depth16( depth_intrin.height,
                         depth_intrin.width,
                         CV_16U,
                         depth_image);
        cv::Mat depth8u = depth16;
        QImage depthimage((uchar*)depth8u.data, depth8u.cols, depth8u.rows, QImage::Format_RGB16);
        ui->label->setPixmap(QPixmap::fromImage(depthimage));
        ui->label->show();
        cvWaitKey(1);
    }
}

void MainWindow::setWidth(int widthIn){
    Width = widthIn;
}

void MainWindow::setLength(int lengthIn){
    Length = lengthIn;
}

void MainWindow::setFPS( int FPSin ){
    FPS = FPSin;
}

int MainWindow::getWidth(){
    return Width;
}

int MainWindow::getLength(){
    return Length;
}

int MainWindow::getFPS(){
    return FPS;
}

void MainWindow::on_StopBtn_clicked()
{
    runBool = false;
}


void MainWindow::on_ExitBtn_clicked()
{
    runBool = false;
    close();
}

// Sliders

void MainWindow::on_laserPwrSlider_sliderMoved(int position)
{
    dev->set_option(rs::option::f200_laser_power, position);
}


void MainWindow::on_accuracySlider_sliderMoved(int position)
{
    dev->set_option(rs::option::f200_accuracy, position);
}


void MainWindow::on_filterOptionSlider_sliderMoved(int position)
{
    dev->set_option(rs::option::f200_filter_option, position);
}



void MainWindow::on_confidenceSlider_sliderMoved(int position)
{
    dev->set_option(rs::option::f200_confidence_threshold, position);
}



void MainWindow::on_motionrangeSlider_sliderMoved(int position)
{
    dev->set_option(rs::option::f200_motion_range, position);
}


void MainWindow::updateSliders()
{
    ui->accuracySlider->setValue(dev->get_option(rs::option::f200_accuracy));
    ui->laserPwrSlider->setValue(dev->get_option(rs::option::f200_laser_power));
    ui->filterOptionSlider->setValue(dev->get_option(rs::option::f200_laser_power));
    ui->confidenceSlider->setValue(dev->get_option(rs::option::f200_confidence_threshold));
    ui->motionrangeSlider->setValue(dev->get_option(rs::option::f200_motion_range));
}

void MainWindow::on_regBtn_clicked()
{
  writeWRL();
}


void MainWindow::writeWRL()
{
    //Prepare wrl file
    std::ostringstream pathString;
    index = getFaceIndex();
    pathString << "data/face_" << index << ".wrl";
    std::ofstream outputWRL;
    outputWRL.open(pathString.str());
    outputWRL << "#VRML V2.0 utf8 " << std::endl;
    outputWRL << "Transform {" << std::endl;
    outputWRL << "scale 1 1 1 " << std::endl;
    outputWRL << "translation 0 0 0 " << std::endl;
    outputWRL << "children[ Shape{ " << std::endl;
    outputWRL << "geometry IndexedFaceSet{" << std::endl;
    outputWRL << "coord Coordinate{" << std::endl;
    outputWRL << "point[" << std::endl;

    float scale = dev->get_depth_scale()*400;
    //
    for(int dy=0; dy<depth_intrin.height; ++dy)
    //for(int dy=depth_intrin.height; dy>0; dy--)
    {
        for(int dx=0; dx<depth_intrin.width; ++dx)
        //for(int dx=depth_intrin.width; dx>0; dx--)
        {

            // Retrieve the 16-bit depth value and map it into a depth in meters
            uint16_t depth_value = depth_image[dy * depth_intrin.width + dx];
            float depth_in_meters = depth_value * scale;
            rs::float2 depth_pixel = {(float)(dx), (float)(depth_intrin.height - 1 - dy)};
            rs::float3 depth_point = depth_intrin.deproject(depth_pixel, depth_in_meters);
            outputWRL << depth_point.x << " " << depth_point.y << " " <<  (-depth_point.z) << ',' << std::endl;
        }
     }
    outputWRL << "]}}}]}" << std::endl;
    outputWRL.close();
}


int MainWindow::getFaceIndex(){
   int index;
   std::fstream indexFile;
   indexFile.open("/home/s1594907/Desktop/index.txt" , std::ios::in);
   if (indexFile.good()){
       indexFile >> index;
       index += 1;
   }
   else{
       index = 1;
   }
   indexFile.close();
   indexFile.open("/home/s1594907/Desktop/index.txt", std::ios::out | std::ios::trunc );
   indexFile << index;
   indexFile.close();
   return index;
}

void MainWindow::on_matchBtn_clicked()
{
//Gallery
    std::fstream probegalleryStream;
    std::ostringstream queryString;
    probegalleryStream.open("gallery.list", std::ios::out | std::ios::trunc);
    probegalleryStream << "data/face_" << (index-1) << ".wrl" << std::endl;
    probegalleryStream.close();
//Probe
    probegalleryStream.open("probe.list", std::ios::out | std::ios::trunc);
    probegalleryStream << "data/face_" << index << ".wrl" << std::endl;
    probegalleryStream.close();
    FILE *fp;
    char path[1035];
    //queryString <<
    /* Open the command for reading. */
    fp = popen("/home/s1594907/FaceRecLuuk/FaceUT3D/faceut3dbi -P probe.list -G gallery.list", "r");
    if (fp == NULL) {
      printf("Failed to run command\n" );
      exit(1);
    }

    /* Read the output a line at a time - output it. */
    while (fgets(path, sizeof(path)-1, fp) != NULL) {
      printf("%s", path);
      ui->textBrowser->setText(path);
      //qDebug(path);
    }

    /* close */
    pclose(fp);

}
