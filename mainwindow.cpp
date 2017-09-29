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


void MainWindow::on_startButton_clicked()
{
//testcommand
    rs::context ctx;
    rs::device * dev = ctx.get_device(0);
    // Configure depth to run at VGA resolution at 30 frames per second
    dev->enable_stream(rs::stream::depth, 640, 480, rs::format::z16, 60);
    dev->start();

    while(true){
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
        //panda
    }
}
