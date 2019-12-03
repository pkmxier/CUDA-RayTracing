TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -Ofast

SOURCES += \
        main.cpp \
    figure.cpp \
    vector3d.cpp \
    color.cpp \
    camera.cpp \
    light.cpp \
    scene.cpp

HEADERS += \
    figure.h \
    vector3d.h \
    color.h \
    camera.h \
    light.h \
    scene.h

DISTFILES += \
    out.ppm \
    images/out0.000000.ppm \
    images/out120.000000.ppm \
    images/out135.000000.ppm \
    images/out150.000000.ppm \
    images/out180.000000.ppm \
    images/out210.000000.ppm \
    images/out225.000000.ppm \
    images/out240.000000.ppm \
    images/out270.000000.ppm \
    images/out30.000000.ppm \
    images/out300.000000.ppm \
    images/out315.000000.ppm \
    images/out330.000000.ppm \
    images/out360.000000.ppm \
    images/out45.000000.ppm \
    images/out60.000000.ppm \
    images/out90.000000.ppm
