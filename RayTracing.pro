TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    figure.cpp \
    vector3d.cpp \
    color.cpp \
    camera.cpp \
    light.cpp

HEADERS += \
    figure.h \
    vector3d.h \
    color.h \
    camera.h \
    light.h

DISTFILES += \
    out.ppm
