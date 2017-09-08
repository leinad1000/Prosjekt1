TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    project1b.cpp \
    project1c.cpp

HEADERS += \
    project1.h

LIBS += -larmadillo -llapack -lblas
