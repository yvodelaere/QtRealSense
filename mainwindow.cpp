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
        // Create depth image
        cv::Mat depth16( depth_intrin.height,
                         depth_intrin.width,
                         CV_16U,
                         (void*)dev->get_frame_data(rs::stream::depth));
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
