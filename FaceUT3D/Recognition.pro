# Created by and for Qt Creator. This file was created for editing the project sources only.
# You may attempt to use it for building too, by modifying this file here.

#TARGET = Recognition

HEADERS = \
   $$PWD/../FaceUT3D/builtinclassifier.h \
   $$PWD/../FaceUT3D/faceut3d.h \
   $$PWD/../lda/asciimatrixio.h \
   $$PWD/../lda/beematrixio.h \
   $$PWD/../lda/lda.h \
   $$PWD/../lda/matrixio.h \
   #$$PWD/../lda/TDFAuthent1_utw.h \
   #$$PWD/../lda/TDFExtract_utw.h \
   #$$PWD/../registration_v3/compat.h \
   $$PWD/../registration_v3/coordlist.h \
   $$PWD/../registration_v3/filter.h \
   $$PWD/../registration_v3/holes.h \
   $$PWD/../registration_v3/nose.h \
   $$PWD/../registration_v3/opt.h \
   $$PWD/../registration_v3/orderedpointset.h \
   #$$PWD/../registration_v3/plane3d.h \
   $$PWD/../registration_v3/point3d.h \
   $$PWD/../registration_v3/rangeimage.h \
   $$PWD/../registration_v3/register.h \
   $$PWD/../registration_v3/roi.h \
   $$PWD/../registration_v3/spline.h \
   $$PWD/../registration_v3/symmetry.h \
   $$PWD/../registration_v3/unorderedpointset.h

SOURCES = \
   $$PWD/../FaceUT3D/builtinclassifier.cc \
   $$PWD/../FaceUT3D/faceut3d.cc \
   $$PWD/../FaceUT3D/main.cc \
   #$$PWD/../FaceUT3D/nobuiltinclassifier.cc \
   $$PWD/../lda/asciimatrixio.cc \
   $$PWD/../lda/beematrixio.cc \
   $$PWD/../lda/lda.cc \
   #$$PWD/../lda/main.cc \
   $$PWD/../lda/matrixio.cc \
   #$$PWD/../lda/reconstruct.cc \
   #$$PWD/../lda/standardLDA.cc \
   #$$PWD/../lda/TDFAuthent1_utw.cc \
  #$$PWD/../lda/TDFExtract_utw.cc \
   #$$PWD/../lda/tstlda.cc \
   #$$PWD/../registration_v3/abs2wrl.cc \
   #$$PWD/../registration_v3/absscale.cc \
   $$PWD/../registration_v3/coordlist.cc \
   $$PWD/../registration_v3/filter.cc \
   $$PWD/../registration_v3/holes.cc \
   #$$PWD/../registration_v3/main.cc \
   $$PWD/../registration_v3/nose.cc \
   $$PWD/../registration_v3/opt.cc \
   $$PWD/../registration_v3/orderedpointset.cc \
   #$$PWD/../registration_v3/pnm2sfi.cc \
   $$PWD/../registration_v3/point3d.cc \
   $$PWD/../registration_v3/rangeimage.cc \
   $$PWD/../registration_v3/register.cc \
   $$PWD/../registration_v3/roi.cc \
   #$$PWD/../registration_v3/sfi2pgm.cc \
   #$$PWD/../registration_v3/sficurv.cc \
   #$$PWD/../registration_v3/sfidiff.cc \
   $$PWD/../registration_v3/spline.cc \
   $$PWD/../registration_v3/symmetry.cc \
   #$$PWD/../registration_v3/transformabs.cc \
   #$$PWD/../registration_v3/tst.cc \
   $$PWD/../registration_v3/unorderedpointset.cc

INCLUDEPATH = \
    $$PWD/../FaceUT3D \
    $$PWD/../lda \
    $$PWD/../registration_v3
INCLUDEPATH += /usr/local/include/opencv
INCLUDEPATH += /usr/local/include
INCLUDEPATH += /usr/include/x86_64-linux-gnu






LIBS += -L/usr/local/lib -lopencv_dnn -lopencv_ml -lopencv_objdetect -lopencv_shape -lopencv_stitching -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_imgproc -lopencv_flann -lopencv_core -lm -ljpeg -lpthread -lrt



#DEFINES =

