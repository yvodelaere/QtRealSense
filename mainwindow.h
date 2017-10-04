#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <librealsense/rs.hpp>
#include <cstdio>
#include <opencv2/opencv.hpp>
#include <QDebug>
#include <string>
#include <sstream>
#include <QString>
#include <fstream>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    rs_error * e = 0;
    rs::intrinsics   depth_intrin;
    void setWidth(int Width);
    void setLength(int Length);
    void setFPS(int FPS);
    int getWidth();
    int getLength();
    int getFPS();
    void updateSliders();
    void writeWRL();
    int getFaceIndex();
private:
    Ui::MainWindow *ui;
    int Width;
    int Length;
    int FPS;
    uint16_t * depth_image;
    bool runBool;
    rs::context ctx;
    rs::device * dev;
    int index;
private slots:
    void on_startBtn_clicked();
    void on_StopBtn_clicked();
    void on_ExitBtn_clicked();
    void on_laserPwrSlider_sliderMoved(int position);
    void on_accuracySlider_sliderMoved(int position);
    void on_filterOptionSlider_sliderMoved(int position);
    void on_confidenceSlider_sliderMoved(int position);
    void on_motionrangeSlider_sliderMoved(int position);
    void on_regBtn_clicked();
    void on_matchBtn_clicked();
};

#endif // MAINWINDOW_H
