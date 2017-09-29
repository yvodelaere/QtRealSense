#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <librealsense/rs.hpp>
#include <cstdio>
#include <opencv2/opencv.hpp>

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

private:
    Ui::MainWindow *ui;

private slots:
    void on_startButton_clicked();
};

#endif // MAINWINDOW_H
