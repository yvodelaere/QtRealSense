#-------------------------------------------------
#
# Project created by QtCreator 2017-09-27T15:42:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ShowCam
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

HEADERS += \
        $$PWD/FaceUT3D/builtinclassifier.h \
        $$PWD/FaceUT3D/faceut3d.h \
        $$PWD/lda/asciimatrixio.h \
        $$PWD/lda/beematrixio.h \
        $$PWD/lda/lda.h \
        $$PWD/lda/matrixio.h \
        $$PWD/registration_v3/coordlist.h \
        $$PWD/registration_v3/filter.h \
        $$PWD/registration_v3/holes.h \
        $$PWD/registration_v3/nose.h \
        $$PWD/registration_v3/opt.h \
        $$PWD/registration_v3/orderedpointset.h \
        $$PWD/registration_v3/point3d.h \
        $$PWD/registration_v3/rangeimage.h \
        $$PWD/registration_v3/register.h \
        $$PWD/registration_v3/roi.h \
        $$PWD/registration_v3/spline.h \
        $$PWD/registration_v3/symmetry.h \
        $$PWD/registration_v3/unorderedpointset.h \
        mainwindow.h

SOURCES += \
        $$PWD/FaceUT3D/builtinclassifier.cc \
        $$PWD/FaceUT3D/faceut3d.cc \
        #$$PWD/FaceUT3D/main.cc \
        $$PWD/lda/asciimatrixio.cc \
        $$PWD/lda/beematrixio.cc \
        $$PWD/lda/lda.cc \
        $$PWD/lda/matrixio.cc \
        $$PWD/registration_v3/coordlist.cc \
        $$PWD/registration_v3/filter.cc \
        $$PWD/registration_v3/holes.cc \
        $$PWD/registration_v3/nose.cc \
        $$PWD/registration_v3/opt.cc \
        $$PWD/registration_v3/orderedpointset.cc \
        $$PWD/registration_v3/point3d.cc \
        $$PWD/registration_v3/rangeimage.cc \
        $$PWD/registration_v3/register.cc \
        $$PWD/registration_v3/roi.cc \
        $$PWD/registration_v3/spline.cc \
        $$PWD/registration_v3/symmetry.cc \
        $$PWD/registration_v3/unorderedpointset.cc \
        main.cpp \
        mainwindow.cpp



FORMS += \
        mainwindow.ui

INCLUDEPATH += /usr/local/include/opencv
INCLUDEPATH += /usr/local/include
INCLUDEPATH += /usr/include/x86_64-linux-gnu
INCLUDEPATH +=    $$PWD/FaceUT3D
INCLUDEPATH +=    $$PWD/lda
INCLUDEPATH +=    $$PWD/registration_v3

LIBS += -L/usr/local/lib -lopencv_dnn -lopencv_ml -lopencv_objdetect -lopencv_shape -lopencv_stitching -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_imgproc -lopencv_flann -lopencv_core -lm -ljpeg -lpthread -lrt -lrealsense -lopencv_core -lopencv_imgproc -lopencv_highgui
