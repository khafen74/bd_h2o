#-------------------------------------------------
#
# Project created by QtCreator 2015-10-03T09:02:39
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = bdStorage2
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    raster.cpp \
    geometry.cpp \
    random.cpp \
    statistics.cpp \
    storagemodel.cpp \
    dampoints.cpp \
    dampolygons.cpp \
    vectorops.cpp \
    reachlines.cpp \
    raster_beaverpond.cpp

win32: LIBS += -L$$PWD/../../../../../../../MinGW/msys/1.0/local/lib/ -llibgdal

INCLUDEPATH += $$PWD/../../../../../../../MinGW/msys/1.0/local/include
DEPENDPATH += $$PWD/../../../../../../../MinGW/msys/1.0/local/include

HEADERS += \
    raster.h \
    geometry.h \
    random.h \
    statistics.h \
    storagemodel.h \
    dampoints.h \
    dampolygons.h \
    vectorops.h \
    reachlines.h \
    raster_beaverpond.h
